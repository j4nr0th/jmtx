//
// Created by jan on 3.1.2024.
//
#include "dense_row_major_internal.h"
#include <assert.h>

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
jmtx_result jmtxc_matrix_drm_new(jmtxc_matrix_drm **p_mtx, uint32_t rows, uint32_t cols,
                                 const _Complex float *set_value, const jmtx_allocator_callbacks *allocator_callbacks)
{
    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }

    jmtxc_matrix_drm *mtx = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*mtx));
    if (!mtx)
    {
        return JMTX_RESULT_BAD_ALLOC;
    }

    memset(mtx, 0xCC, sizeof *mtx);
    mtx->base.cols = cols;
    mtx->base.type = JMTX_TYPE_DRM;
    mtx->base.rows = rows;
    mtx->base.allocator_callbacks = *allocator_callbacks;
    const uint64_t entry_count = rows * cols;

    _Complex float *const values =
        allocator_callbacks->alloc(allocator_callbacks->state, entry_count * sizeof(*values));
    if (!values)
    {
        allocator_callbacks->free(allocator_callbacks->state, mtx);
        return JMTX_RESULT_BAD_ALLOC;
    }

    if (set_value)
    {
        const _Complex float v = *set_value;
        if (v == 0.0f)
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
void jmtxc_matrix_drm_destroy(jmtxc_matrix_drm *mtx)
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
void jmtxc_matrix_drm_set_row(const jmtxc_matrix_drm *mtx, uint32_t row, const _Complex float values[])
{
    _Complex float *ptr = mtx->values + mtx->base.cols * (mtx->permutations ? mtx->permutations[row] : row);
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
void jmtxc_matrix_drm_set_col(const jmtxc_matrix_drm *mtx, uint32_t col, const _Complex float values[])
{
    if (mtx->permutations)
    {
        _Complex float *const ptr = mtx->values + col;
        for (uint32_t i = 0; i < mtx->base.rows; ++i)
        {
            ptr[mtx->base.cols * mtx->permutations[i]] = values[i];
        }
    }
    else
    {
        _Complex float *const ptr = mtx->values + col;
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
uint_fast32_t jmtxc_matrix_drm_get_row(const jmtxc_matrix_drm *mtx, uint32_t row, _Complex float *p_elements[1])
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
void jmtxc_matrix_drm_vector_multiply(const jmtxc_matrix_drm *mtx, const _Complex float *restrict x,
                                      _Complex float *restrict y)
{
    if (mtx->permutations)
    {
        for (uint32_t i = 0; i < mtx->base.rows; ++i)
        {
            const _Complex float *const ptr = mtx->values + mtx->base.cols * mtx->permutations[i];
            _Complex float v = 0;
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
            const _Complex float *const ptr = mtx->values + mtx->base.cols * i;
            _Complex float v = 0;
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
void jmtxc_matrix_drm_set_entry(const jmtxc_matrix_drm *mtx, uint32_t i, uint32_t j, _Complex float value)
{
    if (mtx->permutations)
    {
        _Complex float *const ptr = mtx->values + mtx->permutations[i] * mtx->base.cols;
        ptr[j] = value;
    }
    else
    {
        _Complex float *const ptr = mtx->values + i * mtx->base.cols;
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
_Complex float jmtxc_matrix_drm_get_entry(const jmtxc_matrix_drm *mtx, uint32_t i, uint32_t j)
{
    if (mtx->permutations)
    {
        const _Complex float *const ptr = mtx->values + mtx->permutations[i] * mtx->base.cols;
        return ptr[j];
    }
    const _Complex float *const ptr = mtx->values + i * mtx->base.cols;
    return ptr[j];
}

/**
 * Returns a pointer to an entry in the matrix.
 * @param mtx pointer to the memory where the matrix is stored
 * @param i row index
 * @param j column index
 * @return pointer to an entry
 */
_Complex float *jmtxc_matrix_drm_entry_ptr(const jmtxc_matrix_drm *mtx, uint32_t i, uint32_t j)
{
    if (mtx->permutations)
    {
        _Complex float *const ptr = mtx->values + mtx->permutations[i] * mtx->base.cols;
        return ptr + j;
    }
    _Complex float *const ptr = mtx->values + i * mtx->base.cols;
    return ptr + j;
}

/**
 * Zeros all entries within a matrix, but does not remove them in case they need to be reused
 * @param mtx matrix to zero
 */
void jmtxc_matrix_drm_zero_all_entries(const jmtxc_matrix_drm *mtx)
{
    memset(mtx->values, 0, sizeof(*mtx->values) * mtx->base.cols * mtx->base.rows);
}

/**
 * Similar to jmtxc_matrix_drm_set_all_entries, but slower, since it can not use memset. On the other hand, it allows
 * for the value to be other than 0
 * @param mtx matrix to set
 * @param x value to which to set all entries to
 */
void jmtxc_matrix_drm_set_all_entries(const jmtxc_matrix_drm *mtx, _Complex float x)
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
 * @param values a buffer of at least n values which receives the values of the column
 * @return number of entries that were extracted from the column (may be less than are really in the column if n was too
 * small)
 */
uint32_t jmtxc_matrix_drm_get_col(const jmtxc_matrix_drm *mtx, uint32_t col, _Complex float values[])
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
jmtx_result jmtxc_matrix_drm_transpose(const jmtxc_matrix_drm *mtx, jmtxc_matrix_drm **p_out,
                                       const jmtx_allocator_callbacks *allocator_callbacks)
{
    jmtxc_matrix_drm *new;
    const uint32_t new_cols = mtx->base.rows, new_rows = mtx->base.cols;
    jmtx_result res = jmtxc_matrix_drm_new(&new, new_rows, new_cols, NULL, allocator_callbacks);
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
jmtx_result jmtxc_matrix_drm_transpose_inplace(jmtxc_matrix_drm *mtx, _Complex float *aux_row)
{
    const uint32_t new_cols = mtx->base.rows, new_rows = mtx->base.cols;

    if (mtx->permutations)
    {
        //  Commit to permutations now
        jmtxc_matrix_drm_commit_permutations2(mtx, aux_row);
    }
    assert(mtx->permutations == NULL);
    assert(mtx->rperm == NULL);
    for (uint32_t i = 0; i < (mtx->base.rows + 1) / 2; ++i)
    {
        for (uint32_t j = 0; j < mtx->base.cols; ++j)
        {
            const _Complex float tmp = mtx->values[j * mtx->base.cols + i];
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
jmtx_result jmtxc_matrix_drm_copy(const jmtxc_matrix_drm *mtx, jmtxc_matrix_drm **p_out,
                                  const jmtx_allocator_callbacks *allocator_callbacks)
{
    jmtxc_matrix_drm *new;
    jmtx_result res = jmtxc_matrix_drm_new(&new, mtx->base.rows, mtx->base.cols, NULL, allocator_callbacks);
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
_Complex float jmtxc_matrix_drm_vector_multiply_row(const jmtxc_matrix_drm *mtx, const _Complex float *x, uint32_t i)
{
    const _Complex float *ptr;
    if (mtx->permutations)
    {
        ptr = mtx->values + mtx->base.cols * mtx->permutations[i];
    }
    else
    {
        ptr = mtx->values + mtx->base.cols * i;
    }
    _Complex float v = 0;
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
void jmtxc_matrix_drm_get_bandwidths(const jmtxc_matrix_drm *mtx, uint32_t *ubw, uint32_t *lbw)
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
        const _Complex float *const row_ptr =
            mtx->values + mtx->base.cols * (mtx->permutations ? mtx->permutations[i] : i);
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
jmtx_result jmtxc_matrix_drm_swap_rows(jmtxc_matrix_drm *mtx, uint32_t row1, uint32_t row2)
{
    if (!mtx->permutations)
    {
        mtx->permutations = mtx->base.allocator_callbacks.alloc(mtx->base.allocator_callbacks.state,
                                                                sizeof(*mtx->permutations) * mtx->base.rows);
        if (!mtx->permutations)
        {
            return JMTX_RESULT_BAD_ALLOC;
        }

        mtx->rperm = mtx->base.allocator_callbacks.alloc(mtx->base.allocator_callbacks.state,
                                                         sizeof(*mtx->rperm) * mtx->base.rows);
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
jmtx_result jmtxc_matrix_drm_set_permutation(jmtxc_matrix_drm *mtx, const uint32_t *perm)
{
    if (!mtx->permutations)
    {
        mtx->permutations = mtx->base.allocator_callbacks.alloc(mtx->base.allocator_callbacks.state,
                                                                sizeof(*mtx->permutations) * mtx->base.rows);
        if (!mtx->permutations)
        {
            return JMTX_RESULT_BAD_ALLOC;
        }

        mtx->rperm = mtx->base.allocator_callbacks.alloc(mtx->base.allocator_callbacks.state,
                                                         sizeof(*mtx->rperm) * mtx->base.rows);
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
void jmtxc_matrix_drm_commit_permutations(jmtxc_matrix_drm *mtx)
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
            _Complex float *const p1 = mtx->values + n * dst;
            _Complex float *const p2 = mtx->values + n * src;
            for (uint32_t i = 0; i < n; ++i)
            {
                const _Complex float tmp = p1[i];
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
            _Complex float aux = mtx->values[n * src + i];
            while (mtx->rperm[src] != pos)
            {
                const _Complex float tmp1 = aux;
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
void jmtxc_matrix_drm_commit_permutations2(jmtxc_matrix_drm *mtx, _Complex float *aux_row)
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
            _Complex float *const p1 = mtx->values + n * dst;
            _Complex float *const p2 = mtx->values + n * src;
            for (uint32_t i = 0; i < n; ++i)
            {
                const _Complex float tmp = p1[i];
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
                const _Complex float tmp = aux_row[i];
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
/**
 * Computes the result of a Givens rotation being applied on a (square) matrix as multiplication on the left.
 * The rotation is characterized by a rotation angle theta and two indices, which indicate which two rows the rotation
 * is applied to.
 *
 * This function takes cos(theta) and sin(theta) instead of just the angle directly, because in some cases sine and
 * cosine may be computed directly without computing the angle. In that case it would be redundant to convert those into
 * an angle, then convert them back to sine and cosine.
 *
 * @param mtx input matrix, to which the Givens rotation is applied to
 * @param r1 first integer characterizing the rotation
 * @param r2 second integer characterizing the rotation
 * @param ct value of cos(theta), which is used for the rotation
 * @param st value of sin(theta), which is used for the rotation
 */
void jmtxc_matrix_drm_givens_rotation_left(jmtxc_matrix_drm *mtx, unsigned r1, unsigned r2, const _Complex float ct,
                                           const _Complex float st)
{
    const unsigned n = mtx->base.cols;
    if (mtx->permutations)
    {
        r1 = mtx->permutations[r1];
        r2 = mtx->permutations[r2];
    }
    _Complex float *const row1 = mtx->values + n * r1;
    _Complex float *const row2 = mtx->values + n * r2;

    for (unsigned j = 0; j < n; ++j)
    {
        const _Complex float a = row1[j];
        const _Complex float b = row2[j];
        row1[j] = ct * a - st * b;
        row2[j] = st * a + ct * b;
    }
}

/**
 * Computes the result of a Givens rotation being applied on a (square) matrix as multiplication on the left.
 * The rotation is characterized by a rotation angle theta and two indices, which indicate which two rows the rotation
 * is applied to.
 *
 * This function takes cos(theta) and sin(theta) instead of just the angle directly, because in some cases sine and
 * cosine may be computed directly without computing the angle. In that case it would be redundant to convert those into
 * an angle, then convert them back to sine and cosine.
 *
 * @param mtx input matrix, to which the Givens rotation is applied to
 * @param r1 first integer characterizing the rotation
 * @param r2 second integer characterizing the rotation
 * @param ct value of cos(theta), which is used for the rotation
 * @param st value of sin(theta), which is used for the rotation
 *
 * @return JMTX_RESULT_INDEX_OUT_OF_BOUNDS if either r1 or r2 are out of bounds for the matrix.
 */
jmtx_result jmtxcs_matrix_drm_givens_rotation_left(jmtxc_matrix_drm *mtx, const unsigned r1, const unsigned r2,
                                                   _Complex float ct, _Complex float st)
{
    if (r1 > mtx->base.rows || r2 > mtx->base.rows)
    {
        return JMTX_RESULT_INDEX_OUT_OF_BOUNDS;
    }
    jmtxc_matrix_drm_givens_rotation_left(mtx, r1, r2, ct, st);
    return JMTX_RESULT_SUCCESS;
}

/**
 * Computes the result of a Givens rotation being applied on a (square) matrix as multiplication on the right.
 * The rotation is characterized by a rotation angle theta and two indices, which indicate which two columns the
 * rotation is applied to.
 *
 * This function takes cos(theta) and sin(theta) instead of just the angle directly, because in some cases sine and
 * cosine may be computed directly without computing the angle. In that case it would be redundant to convert those into
 * an angle, then convert them back to sine and cosine.
 *
 * @param mtx input matrix, to which the Givens rotation is applied to
 * @param c1 first integer characterizing the rotation
 * @param c2 second integer characterizing the rotation
 * @param ct value of cos(theta), which is used for the rotation
 * @param st value of sin(theta), which is used for the rotation
 */
void jmtxc_matrix_drm_givens_rotation_right(jmtxc_matrix_drm *mtx, unsigned c1, unsigned c2, _Complex float ct,
                                            _Complex float st)
{
    const unsigned n = mtx->base.rows;

    _Complex float *const col1 = mtx->values + c1;
    _Complex float *const col2 = mtx->values + c2;

    if (!mtx->permutations)
    {
        for (unsigned i = 0; i < n; ++i)
        {
            const _Complex float a = col1[i * mtx->base.cols];
            const _Complex float b = col2[i * mtx->base.cols];
            col1[i * mtx->base.cols] = ct * a - st * b;
            col2[i * mtx->base.cols] = st * a + ct * b;
        }
    }
    else
    {
        for (unsigned i = 0; i < n; ++i)
        {
            const _Complex float a = col1[mtx->permutations[i] * mtx->base.cols];
            const _Complex float b = col2[mtx->permutations[i] * mtx->base.cols];
            col1[mtx->rperm[i] * mtx->base.cols] = ct * a - st * b;
            col2[mtx->rperm[i] * mtx->base.cols] = st * a + ct * b;
        }
    }
}

/**
 * Computes the result of a Givens rotation being applied on a (square) matrix as multiplication on the right.
 * The rotation is characterized by a rotation angle theta and two indices, which indicate which two columns the
 * rotation is applied to.
 *
 * This function takes cos(theta) and sin(theta) instead of just the angle directly, because in some cases sine and
 * cosine may be computed directly without computing the angle. In that case it would be redundant to convert those into
 * an angle, then convert them back to sine and cosine.
 *
 * @param mtx input matrix, to which the Givens rotation is applied to
 * @param c1 first integer characterizing the rotation
 * @param c2 second integer characterizing the rotation
 * @param ct value of cos(theta), which is used for the rotation
 * @param st value of sin(theta), which is used for the rotation
 *
 * @return JMTX_RESULT_INDEX_OUT_OF_BOUNDS if either r1 or r2 are out of bounds for the matrix.
 */
jmtx_result jmtxcs_matrix_drm_givens_rotation_right(jmtxc_matrix_drm *mtx, unsigned c1, unsigned c2, _Complex float ct,
                                                    _Complex float st)
{
    if (c1 > mtx->base.cols || c2 > mtx->base.cols)
    {
        return JMTX_RESULT_INDEX_OUT_OF_BOUNDS;
    }
    jmtxc_matrix_drm_givens_rotation_right(mtx, c1, c2, ct, st);
    return JMTX_RESULT_SUCCESS;
}

/**
 * Computes the matrix product of two matrices A and B as C = A B. This is can be written as:
 * $$
 *  C_{i,j} = \sum\limits_{k=0}^{N-1} A_{i,k} B_{k,j}
 * $$
 *
 * This version of the function uses a pre-exsiting matrix as output.
 *
 * @param a matrix A
 * @param b matrix B
 * @param out where the output should be returned
 * @return JMTX_RESULT_SUCCES if successful, JMTX_RESULT_DIMS_MISMATCH if the dimensions of a and b don't allow
 * multiplication, or if c does not have correct dimensions
 */
jmtx_result jmtxc_matrix_drm_multiply_matrix(jmtxc_matrix_drm *a, jmtxc_matrix_drm *b, jmtxc_matrix_drm *out)
{
    if (out->permutations)
    {
        return JMTX_RESULT_BAD_PARAM;
    }
    const unsigned n = a->base.cols;
    if (n != b->base.rows)
    {
        return JMTX_RESULT_DIMS_MISMATCH;
    }

    const unsigned new_rows = out->base.rows;
    const unsigned new_cols = out->base.cols;
    if (a->base.rows != new_rows || b->base.cols != new_cols)
    {
        return JMTX_RESULT_DIMS_MISMATCH;
    }
    if (!b->permutations && !a->permutations)
    {
        for (unsigned i = 0; i < new_rows; ++i)
        {
            const _Complex float *restrict p_row = a->values + i * n;
            for (unsigned j = 0; j < new_cols; ++j)
            {
                _Complex float v = 0;
                const _Complex float *restrict p_col = b->values + j;
                for (unsigned k = 0; k < n; ++k)
                {
                    v += p_row[k] * p_col[b->base.cols * k];
                    out->values[j + i * new_cols] = v;
                }
            }
        }
    }
    else if (a->permutations && b->permutations)
    {
        for (unsigned i = 0; i < new_rows; ++i)
        {
            const _Complex float *restrict p_row = a->values + n * a->permutations[i];
            for (unsigned j = 0; j < new_cols; ++j)
            {
                _Complex float v = 0;
                const _Complex float *restrict p_col = b->values + j;
                for (unsigned k = 0; k < n; ++k)
                {
                    v += p_row[k] * p_col[b->base.cols * b->permutations[k]];
                    out->values[j + i * new_cols] = v;
                }
            }
        }
    }
    else if (a->permutations)
    {
        for (unsigned i = 0; i < new_rows; ++i)
        {
            const _Complex float *restrict p_row = a->values + n * a->permutations[i];
            for (unsigned j = 0; j < new_cols; ++j)
            {
                _Complex float v = 0;
                const _Complex float *restrict p_col = b->values + j;
                for (unsigned k = 0; k < n; ++k)
                {
                    v += p_row[k] * p_col[b->base.cols * k];
                    out->values[j + i * new_cols] = v;
                }
            }
        }
    }
    else // if (b->permutations)
    {
        for (unsigned i = 0; i < new_rows; ++i)
        {
            const _Complex float *restrict p_row = a->values + n * i;
            for (unsigned j = 0; j < new_cols; ++j)
            {
                _Complex float v = 0;
                const _Complex float *restrict p_col = b->values + j;
                for (unsigned k = 0; k < n; ++k)
                {
                    v += p_row[k] * p_col[b->base.cols * k];
                    out->values[j + i * new_cols] = v;
                }
            }
        }
    }

    return JMTX_RESULT_SUCCESS;
}

/**
 * Shifts the diagonal of the matrix mtx by the value v, so that: A_{i,i} = A_{i,i} + v for all i
 * @param mtx matrix which should have its diagonal shifted
 * @param v value by which to shift the diagonal
 */
void jmtxc_matrix_drm_shift_diagonal(jmtxc_matrix_drm *mtx, _Complex float v)
{
    const unsigned n = mtx->base.rows > mtx->base.cols ? mtx->base.cols : mtx->base.rows;
    if (!mtx->permutations)
    {
        for (unsigned i = 0; i < n; ++i)
        {
            mtx->values[i * (mtx->base.cols + 1)] += v;
        }
    }
    else
    {
        for (unsigned i = 0; i < n; ++i)
        {
            mtx->values[mtx->rperm[i] * (mtx->base.cols) + i] += v;
        }
    }
}

jmtxc_matrix_drm jmtxc_matrix_drm_from_data(const unsigned rows, const unsigned cols,
                                            _Complex float values[JMTX_ARRAY_ATTRIB(static rows * cols)])
{
    return (jmtxc_matrix_drm){.base = {.type = JMTX_TYPE_DRM,
                                       .rows = rows,
                                       .cols = cols,
                                       .allocator_callbacks = {.state = NULL, .alloc = NULL, .free = NULL}},
                              .permutations = NULL,
                              .rperm = NULL,
                              .values = values};
}
