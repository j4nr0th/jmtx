//
// Created by jan on 3.1.2024.
//
#include <assert.h>
#include "dense_row_major_internal.h"


/**
 * Initializes a new Dense Row-Major matrix
 * @param p_mtx address that receives the pointer to the matrix
 * @param rows number of rows of the matrix
 * @param cols number of columns of the matrix
 * @param set_value pointer to value with which to initialize all entries. If NULL, then matrix is left uninitialized
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxd_matrix_drm_new(
        jmtxd_matrix_drm** p_mtx, uint32_t cols, uint32_t rows, const double* set_value,
        const jmtx_allocator_callbacks* allocator_callbacks)
{
    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }


    jmtxd_matrix_drm* mtx = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*mtx));
    if (!mtx)
    {
        return JMTX_RESULT_BAD_ALLOC;
    }


    memset(mtx, 0xCC, sizeof*mtx);
    mtx->base.cols = cols;
    mtx->base.type = JMTX_TYPE_DRM;
    mtx->base.rows = rows;
    mtx->base.allocator_callbacks = *allocator_callbacks;
    const uint64_t entry_count = rows * cols;

    double* const values = allocator_callbacks->alloc(allocator_callbacks->state, entry_count * sizeof(*values));
    if (!values)
    {
        allocator_callbacks->free(allocator_callbacks->state, mtx);
        return JMTX_RESULT_BAD_ALLOC;
    }

    if (set_value)
    {
        const double v = *set_value;
        if (v == 0.0)
        {
            memset(values, 0, sizeof(*values) * entry_count);
        }
        else
        {
            for (uint64_t i = 0; i < entry_count; ++i)
            {
                values[i] = v;
            }
        }
    }
    mtx->values = values;
    mtx->permutations = NULL;
    mtx->rperm = NULL;

    *p_mtx = mtx;

    return JMTX_RESULT_SUCCESS;
}


/**
 * Cleans up the DRM matrix and frees all of its memory
 * @param mtx pointer to memory where the matrix is stored
 */
void jmtxd_matrix_drm_destroy(jmtxd_matrix_drm* mtx)
{
    jmtx_allocator_callbacks callbacks = mtx->base.allocator_callbacks;
    if (mtx->permutations)
    {
        callbacks.free(callbacks.state, mtx->permutations);
    }
    if (mtx->rperm)
    {
        callbacks.free(callbacks.state, mtx->rperm);
    }
    callbacks.free(callbacks.state, mtx->values);
    callbacks.free(callbacks.state, mtx);
}

/**
 * Sets the row of the matrix. More efficient than setting it element by element
 * @param mtx pointer to the memory where the matrix is stored
 * @param row index of the row to set
 * @param values values of entries
 */
void jmtxd_matrix_drm_set_row(const jmtxd_matrix_drm* mtx, uint32_t row, const double values[])
{
    double* ptr = mtx->values + mtx->base.cols * (mtx->permutations ? mtx->permutations[row] : row);
    for (uint32_t i = 0; i < mtx->base.cols; ++i)
    {
        ptr[i] = values[i];
    }
}

/**
 * Sets the column of the matrix. More efficient than setting it element by element
 * @param mtx pointer to the memory where the matrix is stored
 * @param col index of the column to set
 * @param values values of entries
 */
void jmtxd_matrix_drm_set_col(const jmtxd_matrix_drm* mtx, uint32_t col, const double values[])
{
    if (mtx->permutations)
    {
        double* const ptr = mtx->values + col;
        for (uint32_t i = 0; i < mtx->base.rows; ++i)
        {
            ptr[mtx->base.cols * mtx->permutations[i]] = values[i];
        }
    }
    else
    {
        double* const ptr = mtx->values + col;
        for (uint32_t i = 0; i < mtx->base.rows; ++i)
        {
            ptr[mtx->base.cols * i] = values[i];
        }
    }
}

/**
 * Returns the pointers to arrays of column indices and element values for that row
 * @param mtx pointer to the memory where the matrix is stored
 * @param row index of the row to get
 * @param p_elements pointer to array of values
 * @return number of elements in the row, which is the number of valid elements in arrays given to p_indices and
 * p_elements
 */
uint_fast32_t jmtxd_matrix_drm_get_row(const jmtxd_matrix_drm* mtx, uint32_t row, double* p_elements[1])
{
    if (mtx->permutations)
    {
        *p_elements = mtx->values + mtx->base.cols * mtx->permutations[row];
    }
    else
    {
        *p_elements = mtx->values + mtx->base.cols * row;
    }
    return mtx->base.cols;
}

/**
 * Multiplies a dense column vector x by the sparse matrix and stores the result at y
 * @param mtx pointer to the memory where the matrix is stored
 * @param x pointer to vector to be multiplied
 * @param y pointer to vector where the result of multiplication is to be stored
 */
void jmtxd_matrix_drm_vector_multiply(const jmtxd_matrix_drm* mtx, const double* restrict x, double* restrict y)
{
    if (mtx->permutations)
    {
        for (uint32_t i = 0; i < mtx->base.rows; ++i)
        {
            const double* const ptr = mtx->values + mtx->base.cols * mtx->permutations[i];
            double v = 0;
            for (uint32_t j = 0; j < mtx->base.cols; ++j)
            {
                v += x[j] * ptr[j];
            }
            y[i] = v;
        }
    }
    else
    {

        for (uint32_t i = 0; i < mtx->base.rows; ++i)
        {
            const double* const ptr = mtx->values + mtx->base.cols * i;
            double v = 0;
            for (uint32_t j = 0; j < mtx->base.cols; ++j)
            {
                v += x[j] * ptr[j];
            }
            y[i] = v;
        }
    }
}

/**
 * Sets a single entry in the matrix.
 * @param mtx pointer to the memory where the matrix is stored
 * @param i row index
 * @param j column index
 * @param value value to which the value is set
 */
void jmtxd_matrix_drm_set_entry(const jmtxd_matrix_drm* mtx, uint32_t i, uint32_t j, double value)
{
    if (mtx->permutations)
    {
        double* const ptr = mtx->values + mtx->permutations[i] * mtx->base.cols;
        ptr[j] = value;
    }
    else
    {
        double* const ptr = mtx->values + i * mtx->base.cols;
        ptr[j] = value;
    }
}

/**
 * Returns a single entry from the matrix.
 * @param mtx pointer to the memory where the matrix is stored
 * @param i row index
 * @param j column index
 * @return value of the entry
 */
double jmtxd_matrix_drm_get_entry(const jmtxd_matrix_drm* mtx, uint32_t i, uint32_t j)
{
    if (mtx->permutations)
    {
        const double* const ptr = mtx->values + mtx->permutations[i] * mtx->base.cols;
        return ptr[j];
    }
    const double* const ptr = mtx->values + i * mtx->base.cols;
    return ptr[j];
}

/**
 * Returns a pointer to an entry in the matrix.
 * @param mtx pointer to the memory where the matrix is stored
 * @param i row index
 * @param j column index
 * @return pointer to an entry
 */
double* jmtxd_matrix_drm_entry_ptr(const jmtxd_matrix_drm* mtx, uint32_t i, uint32_t j)
{
    if (mtx->permutations)
    {
        double* const ptr = mtx->values + mtx->permutations[i] * mtx->base.cols;
        return ptr + j;
    }
    double* const ptr = mtx->values + i * mtx->base.cols;
    return ptr + j;
}


/**
 * Zeros all entries within a matrix, but does not remove them in case they need to be reused
 * @param mtx matrix to zero
 */
void jmtxd_matrix_drm_zero_all_entries(const jmtxd_matrix_drm* mtx)
{
    memset(mtx->values, 0, sizeof(*mtx->values) * mtx->base.cols * mtx->base.rows);
}

/**
 * Similar to jmtxd_matrix_drm_set_all_entries, but slower, since it can not use memset. On the other hand, it allows for
 * the value to be other than 0
 * @param mtx matrix to set
 * @param x value to which to set all entries to
 */
void jmtxd_matrix_drm_set_all_entries(const jmtxd_matrix_drm* mtx, double x)
{
    for (uint32_t i = 0; i < mtx->base.cols * mtx->base.rows; ++i)
    {
        mtx->values[i] = x;
    }
}

/**
 * Returns the values of entries in the matrix, along with what row of the matrix they were located in
 * @param mtx pointer to the memory where the matrix is stored
 * @param col column index of the matrix to look at
 * @param p_values a buffer of at least n values which receives the values of the column
 * @return number of entries that were extracted from the column (may be less than are really in the column if n was too
 * small)
 */
uint32_t
jmtxd_matrix_drm_get_col(const jmtxd_matrix_drm* mtx, uint32_t col, double values[])
{
    if (mtx->permutations)
    {
        for (uint32_t i = 0; i < mtx->base.rows; ++i)
        {
            values[i] = mtx->values[mtx->permutations[i] * mtx->base.cols + col];
        }
    }
    else
    {
        for (uint32_t i = 0; i < mtx->base.rows; ++i)
        {
            values[i] = mtx->values[i * mtx->base.cols + col];
        }
    }
    return mtx->base.rows;
}

/**
 * Creates a transpose of a matrix
 * @param mtx pointer to the memory where the input matrix is stored
 * @param p_out address where the pointer to the output matrix will be returned
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxd_matrix_drm_transpose(
        const jmtxd_matrix_drm* mtx, jmtxd_matrix_drm** p_out, const jmtx_allocator_callbacks* allocator_callbacks)
{
    jmtxd_matrix_drm* new;
    const uint32_t new_cols = mtx->base.rows, new_rows = mtx->base.cols;
    jmtx_result res = jmtxd_matrix_drm_new(&new, new_cols, new_rows, NULL, allocator_callbacks);
    if (res != JMTX_RESULT_SUCCESS)
    {
        return res;
    }

    if (mtx->permutations)
    {
        for (uint32_t i = 0; i < mtx->base.rows; ++i)
        {
            for (uint32_t j = 0; j < mtx->base.cols; ++j)
            {
                new->values[j * new->base.cols + i] = mtx->values[mtx->permutations[i] * mtx->base.cols + j];
            }
        }
    }
    else
    {
        for (uint32_t i = 0; i < mtx->base.rows; ++i)
        {
            for (uint32_t j = 0; j < mtx->base.cols; ++j)
            {
                new->values[j * new->base.cols + i] = mtx->values[i * mtx->base.cols + j];
            }
        }
    }

    *p_out = new;
    return JMTX_RESULT_SUCCESS;
}

/**
 * Creates a transpose of a matrix
 * @param mtx pointer to the memory where the input matrix is stored
 * @param p_out address where the pointer to the output matrix will be returned
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxd_matrix_drm_transpose_inplace(jmtxd_matrix_drm* mtx, double* aux_row)
{
    const uint32_t new_cols = mtx->base.rows, new_rows = mtx->base.cols;

    if (mtx->permutations)
    {
        //  Commit to permutations now
        jmtxd_matrix_drm_commit_permutations2(mtx, aux_row);
    }
    assert(mtx->permutations == NULL);
    assert(mtx->rperm == NULL);
    for (uint32_t i = 0; i < (mtx->base.rows + 1) / 2; ++i)
    {
        for (uint32_t j = 0; j < mtx->base.cols; ++j)
        {
            const double tmp = mtx->values[j * mtx->base.cols + i];
            mtx->values[j * mtx->base.cols + i] = mtx->values[i * mtx->base.cols + j];
            mtx->values[i * mtx->base.cols + j] = tmp;
        }
    }


    mtx->base.rows = new_rows;
    mtx->base.cols = new_cols;
    return JMTX_RESULT_SUCCESS;
}

/**
 * Creates a copy of the matrix
 * @param mtx pointer to the memory where the input matrix is stored
 * @param p_out address where the pointer to the output matrix will be returned
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxd_matrix_drm_copy(const jmtxd_matrix_drm* mtx, jmtxd_matrix_drm** p_out, const jmtx_allocator_callbacks* allocator_callbacks)
{
    jmtxd_matrix_drm* new;
    jmtx_result res = jmtxd_matrix_drm_new(&new, mtx->base.cols, mtx->base.rows, NULL, allocator_callbacks);
    if (res != JMTX_RESULT_SUCCESS)
    {
        return res;
    }

    if (mtx->permutations)
    {
        for (uint32_t i = 0; i < mtx->base.rows; ++i)
        {
            for (uint32_t j = 0; j < mtx->base.cols; ++j)
            {
                new->values[j * new->base.cols + i] = mtx->values[mtx->permutations[i] * mtx->base.cols + j];
            }
        }
    }
    else
    {
        for (uint32_t i = 0; i < mtx->base.rows; ++i)
        {
            for (uint32_t j = 0; j < mtx->base.cols; ++j)
            {
                new->values[j * new->base.cols + i] = mtx->values[i * mtx->base.cols + j];
            }
        }
    }

    *p_out = new;
    return JMTX_RESULT_SUCCESS;
}


/**
 * Computes one entry of Ax. This function only computes the i-th entry to make it possible to compute it in parallel.
 * @param mtx pointer to the memory where the matrix A is stored.
 * @param x pointer to the memory where the vector x is stored
 * @param i what entry of the residual to compute
 * @return result of the multiplication
 */
double jmtxd_matrix_drm_vector_multiply_row(const jmtxd_matrix_drm* mtx, const double* x, uint32_t i)
{
    const double* ptr;
    if (mtx->permutations)
    {
        ptr = mtx->values + mtx->base.cols * mtx->permutations[i];
    }
    else
    {
        ptr = mtx->values + mtx->base.cols * i;
    }
    double v = 0;
    for (uint32_t j = 0; j < mtx->base.cols; ++j)
    {
        v += ptr[i] * x[i];
    }
    return v;
}

/**
 * Returns the upper bandwidth and the lower bandwidth of the BRM matrix
 * @param mtx pointer to the memory where the matrix is stored
 * @param ubw pointer which receives the upper bandwidth of the matrix
 * @param lbw pointer which receives the lower bandwidth of the matrix
 */
void jmtxd_matrix_drm_get_bandwidths(const jmtxd_matrix_drm* mtx, uint32_t* ubw, uint32_t* lbw)
{
    uint32_t l = 0, u = 0;
    //  Quick check for extreme case
    if (mtx->values[mtx->base.cols - 1] != 0 && mtx->values[mtx->base.cols * (mtx->base.rows - 1)] != 0)
    {
        *ubw = mtx->base.cols - 1;
        *lbw = mtx->base.rows - 1;
        return;
    }

    for (uint32_t i = 0; i < mtx->base.rows; ++i)
    {
        const double* const row_ptr = mtx->values + mtx->base.cols * (mtx->permutations ? mtx->permutations[i] : i);
        //  Check if ubw is greater
        for (uint32_t j = i + u; j < mtx->base.cols; ++j)
        {
            if (row_ptr[j] != 0)
            {
                u = j - i;
            }
        }
        //  Check if lbw is greater
        for (uint32_t j = i - l; j > 0; --j)
        {
            if (row_ptr[j - 1] != 0)
            {
                l = i - (j - 1);
            }
        }
    }
    *lbw = l;
    *ubw = u;
}

/**
 * Exchanges row of matrix through the use of a permutation matrix. This does not cause memory copying and is thus fast
 * @param mtx pointer to the memory where the matrix is stored
 * @param row1 index of a row to exchange with row2
 * @param row2 index of a row to exchange with row1
 */
jmtx_result jmtxd_matrix_drm_swap_rows(jmtxd_matrix_drm* mtx, uint32_t row1, uint32_t row2)
{
    if (!mtx->permutations)
    {
        mtx->permutations = mtx->base.allocator_callbacks.alloc(mtx->base.allocator_callbacks.state, sizeof(*mtx->permutations) * mtx->base.rows);
        if (!mtx->permutations)
        {
            return JMTX_RESULT_BAD_ALLOC;
        }

        mtx->rperm = mtx->base.allocator_callbacks.alloc(mtx->base.allocator_callbacks.state, sizeof(*mtx->rperm) * mtx->base.rows);
        if (!mtx->rperm)
        {
            mtx->base.allocator_callbacks.free(mtx->base.allocator_callbacks.state, mtx->permutations);
            mtx->permutations = NULL;
            return JMTX_RESULT_BAD_ALLOC;
        }
        for (uint32_t i = 0; i < mtx->base.rows; ++i)
        {
            mtx->permutations[i] = i;
            mtx->rperm[i] = i;
        }
    }

    //  Swap two rows in permutation list
    const uint32_t tmp = mtx->permutations[row1];
    mtx->permutations[row1] = mtx->permutations[row2];
    mtx->permutations[row2] = tmp;
    //  Correct reverse permutation lists
    const uint32_t rtmp = mtx->rperm[mtx->permutations[row1]];
    mtx->rperm[mtx->permutations[row1]] = mtx->rperm[mtx->permutations[row2]];
    mtx->rperm[mtx->permutations[row2]] = rtmp;
    return JMTX_RESULT_SUCCESS;
}


/**
 * Permutes all rows at once. The provided list must contain all entries in the set [0, n) exactly once, where n is the
 * number or rows of the matrix.
 * @param mtx pointer to the memory where the matrix is stored
 * @param perm list of permutation indices
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on allocation failure
 */
jmtx_result jmtxd_matrix_drm_set_permutation(jmtxd_matrix_drm* mtx, const uint32_t* perm)
{
    if (!mtx->permutations)
    {
        mtx->permutations = mtx->base.allocator_callbacks.alloc(mtx->base.allocator_callbacks.state, sizeof(*mtx->permutations) * mtx->base.rows);
        if (!mtx->permutations)
        {
            return JMTX_RESULT_BAD_ALLOC;
        }

        mtx->rperm = mtx->base.allocator_callbacks.alloc(mtx->base.allocator_callbacks.state, sizeof(*mtx->rperm) * mtx->base.rows);
        if (!mtx->rperm)
        {
            mtx->base.allocator_callbacks.free(mtx->base.allocator_callbacks.state, mtx->permutations);
            mtx->permutations = NULL;
            return JMTX_RESULT_BAD_ALLOC;
        }
    }
    for (uint32_t i = 0; i < mtx->base.rows; ++i)
    {
        mtx->permutations[i] = perm[i];
    }
    //  Deal with reverse permutations
    uint32_t first = 0, last = mtx->base.rows - 1;
    for (uint32_t i = 0; i < mtx->base.rows; ++i)
    {
        //  Find j such that perm[j] == i
        for (uint32_t j = first; j <= last; ++j)
        {
            if (perm[j] == i)
            {
                if (j == first)
                {
                    first += 1;
                }
                else if (j == last)
                {
                    last -= 1;
                }
                mtx->rperm[i] = j;
                break;
            }
        }
    }
    return JMTX_RESULT_SUCCESS;
}

/**
 * Executes permutations of matrix rows and reorders rows in memory. This may allow for better memory access on a
 * permuted matrix. Should not be done on decompositions. This function requires additional memory, but that allows
 * it to swap whole rows at once, which makes it more cache friendly and potentially faster for large matrices.
 * @param mtx pointer to the memory where the matrix is stored
 * @param aux_row memory that can be used by the function to store an intermediate matrix row
 */
void jmtxd_matrix_drm_commit_permutations(jmtxd_matrix_drm* mtx)
{
    uint32_t pos;
    const uint32_t n = mtx->base.cols;
    for (pos = 0; pos < mtx->base.rows; ++pos)
    {
        uint32_t src = pos;
        uint32_t dst = mtx->rperm[pos];
        if (src == dst)
        {
            //  Already there
            continue;
        }
        if (mtx->rperm[dst] == src)
        {
            assert(mtx->permutations[src] == dst);
            //  Simple swap
            double* const p1 = mtx->values + n * dst;
            double* const p2 = mtx->values + n * src;
            for (uint32_t i = 0; i < n; ++i)
            {
                const double tmp = p1[i];
                p1[i] = p2[i];
                p2[i] = tmp;
            }
            mtx->permutations[src] = src;
            mtx->permutations[dst] = dst;
            mtx->rperm[src] = src;
            mtx->rperm[dst] = dst;
            continue;
        }
        //  This has to be done the hard way
        uint32_t v = 0;
        while (mtx->rperm[src] != pos)
        {
            uint32_t old_dst = dst;
            dst = mtx->rperm[dst];
            src = old_dst;
            v += 1;
        }

        for (uint32_t i = 0; i < n; ++i)
        {
            dst = mtx->rperm[pos];
            src = pos;
            double aux = mtx->values[n * src + i];
            while (mtx->rperm[src] != pos)
            {
                const double tmp1 = aux;
                aux = mtx->values[n * dst + i];
                mtx->values[n * dst + i] = tmp1;
                uint32_t old_dst = dst;
                dst = mtx->rperm[dst];
                src = old_dst;
            }
            mtx->values[n * pos + i] = aux;
        }
        dst = mtx->rperm[pos];
        for (uint32_t i = 0; i < v; ++i)
        {
            const uint32_t old_dst = dst;
            dst = mtx->rperm[dst];
            mtx->rperm[old_dst] = old_dst;
        }
        mtx->rperm[pos] = pos;
    }
#ifndef NDEBUG
    for (uint32_t i = 0; i < mtx->base.rows; ++i)
    {
        assert(mtx->rperm[i] == i);
    }
#endif
    mtx->base.allocator_callbacks.free(mtx->base.allocator_callbacks.state, mtx->permutations);
    mtx->base.allocator_callbacks.free(mtx->base.allocator_callbacks.state, mtx->rperm);
    mtx->permutations = NULL;
    mtx->rperm = NULL;
}

/**
 * Executes permutations of matrix rows and reorders rows in memory. This may allow for better memory access on a
 * permuted matrix. Should not be done on decompositions. This function requires additional memory, but that allows
 * it to swap whole rows at once, which makes it more cache friendly and potentially faster for large matrices.
 * @param mtx pointer to the memory where the matrix is stored
 * @param aux_row memory that can be used by the function to store an intermediate matrix row
 */
void jmtxd_matrix_drm_commit_permutations2(jmtxd_matrix_drm* mtx, double* aux_row)
{
    uint32_t pos;
    const uint32_t n = mtx->base.cols;
    for (pos = 0; pos < mtx->base.rows; ++pos)
    {
        uint32_t src = pos;
        uint32_t dst = mtx->rperm[pos];
        if (src == dst)
        {
            //  Already there
            continue;
        }
        if (mtx->rperm[dst] == src)
        {
            assert(mtx->permutations[src] == dst);
            //  Simple swap
            double* const p1 = mtx->values + n * dst;
            double* const p2 = mtx->values + n * src;
            for (uint32_t i = 0; i < n; ++i)
            {
                const double tmp = p1[i];
                p1[i] = p2[i];
                p2[i] = tmp;
            }
            mtx->permutations[src] = src;
            mtx->permutations[dst] = dst;
            mtx->rperm[src] = src;
            mtx->rperm[dst] = dst;
            continue;
        }
        //  This has to be done the hard way
        memcpy(aux_row, mtx->values + n * src, sizeof(*aux_row) * n);
        uint32_t v = 0;
        while (mtx->rperm[src] != pos)
        {
            for (uint32_t i = 0; i < n; ++i)
            {
                const double tmp = aux_row[i];
                aux_row[i] = mtx->values[n * dst + i];
                mtx->values[n * dst + i] = tmp;
            }
            uint32_t old_dst = dst;
            dst = mtx->rperm[dst];
            src = old_dst;
            v += 1;
        }
        memcpy(mtx->values + n * pos, aux_row, sizeof(*aux_row) * n);
        dst = mtx->rperm[pos];
        for (uint32_t i = 0; i < v; ++i)
        {
            const uint32_t old_dst = dst;
            dst = mtx->rperm[dst];
            mtx->rperm[old_dst] = old_dst;
        }
        mtx->rperm[pos] = pos;
    }
#ifndef NDEBUG
    for (uint32_t i = 0; i < mtx->base.rows; ++i)
    {
        assert(mtx->rperm[i] == i);
    }
#endif
    mtx->base.allocator_callbacks.free(mtx->base.allocator_callbacks.state, mtx->permutations);
    mtx->base.allocator_callbacks.free(mtx->base.allocator_callbacks.state, mtx->rperm);
    mtx->permutations = NULL;
    mtx->rperm = NULL;
}

