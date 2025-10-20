#include <assert.h>
#include <math.h>

#include "dense_row_major.h"

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
jmtx_result JMTX_NAME_TYPED(matrix_drm_new)(JMTX_NAME_TYPED(matrix_drm) * *p_mtx, JMTX_INDEX_T rows, JMTX_INDEX_T cols,
                                            const JMTX_SCALAR_T *set_value,
                                            const jmtx_allocator_callbacks *allocator_callbacks)
{
    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }

    JMTX_NAME_TYPED(matrix_drm) *mtx = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*mtx));
    if (!mtx)
    {
        return JMTX_RESULT_BAD_ALLOC;
    }

    memset(mtx, 0xCC, sizeof *mtx);
    mtx->base.cols = cols;
    mtx->base.type = JMTXF_TYPE_DRM;
    mtx->base.rows = rows;
    mtx->base.allocator_callbacks = *allocator_callbacks;
    const uint64_t entry_count = rows * cols;

    JMTX_SCALAR_T *const values = allocator_callbacks->alloc(allocator_callbacks->state, entry_count * sizeof(*values));
    if (!values)
    {
        allocator_callbacks->free(allocator_callbacks->state, mtx);
        return JMTX_RESULT_BAD_ALLOC;
    }

    if (set_value)
    {
        const JMTX_SCALAR_T v = *set_value;
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
void JMTX_NAME_TYPED(matrix_drm_destroy)(JMTX_NAME_TYPED(matrix_drm) * mtx)
{
    const jmtx_allocator_callbacks callbacks = mtx->base.allocator_callbacks;
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
void JMTX_NAME_TYPED(matrix_drm_set_row)(const JMTX_NAME_TYPED(matrix_drm) * mtx, JMTX_INDEX_T row,
                                         const JMTX_SCALAR_T values[])
{
    JMTX_SCALAR_T *ptr = mtx->values + mtx->base.cols * (mtx->permutations ? mtx->permutations[row] : row);
    for (JMTX_INDEX_T i = 0; i < mtx->base.cols; ++i)
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
void JMTX_NAME_TYPED(matrix_drm_set_col)(const JMTX_NAME_TYPED(matrix_drm) * mtx, JMTX_INDEX_T col,
                                         const JMTX_SCALAR_T values[])
{
    if (mtx->permutations)
    {
        JMTX_SCALAR_T *const ptr = mtx->values + col;
        for (JMTX_INDEX_T i = 0; i < mtx->base.rows; ++i)
        {
            ptr[mtx->base.cols * mtx->permutations[i]] = values[i];
        }
    }
    else
    {
        JMTX_SCALAR_T *const ptr = mtx->values + col;
        for (JMTX_INDEX_T i = 0; i < mtx->base.rows; ++i)
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
JMTX_FAST_INT_T JMTX_NAME_TYPED(matrix_drm_get_row)(const JMTX_NAME_TYPED(matrix_drm) * mtx, JMTX_INDEX_T row,
                                                    JMTX_SCALAR_T *p_elements[1])
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
void JMTX_NAME_TYPED(matrix_drm_vector_multiply)(const JMTX_NAME_TYPED(matrix_drm) * mtx,
                                                 const JMTX_SCALAR_T *restrict x, JMTX_SCALAR_T *restrict y)
{
    if (mtx->permutations)
    {
        for (JMTX_INDEX_T i = 0; i < mtx->base.rows; ++i)
        {
            const JMTX_SCALAR_T *const ptr = mtx->values + mtx->base.cols * mtx->permutations[i];
            JMTX_SCALAR_T v = 0;
            for (JMTX_INDEX_T j = 0; j < mtx->base.cols; ++j)
            {
                v += x[j] * ptr[j];
            }
            y[i] = v;
        }
    }
    else
    {

        for (JMTX_INDEX_T i = 0; i < mtx->base.rows; ++i)
        {
            const JMTX_SCALAR_T *const ptr = mtx->values + mtx->base.cols * i;
            JMTX_SCALAR_T v = 0;
            for (JMTX_INDEX_T j = 0; j < mtx->base.cols; ++j)
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
void JMTX_NAME_TYPED(matrix_drm_set_entry)(const JMTX_NAME_TYPED(matrix_drm) * mtx, JMTX_INDEX_T i, JMTX_INDEX_T j,
                                           JMTX_SCALAR_T value)
{
    if (mtx->permutations)
    {
        JMTX_SCALAR_T *const ptr = mtx->values + mtx->permutations[i] * mtx->base.cols;
        ptr[j] = value;
    }
    else
    {
        JMTX_SCALAR_T *const ptr = mtx->values + i * mtx->base.cols;
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
JMTX_SCALAR_T JMTX_NAME_TYPED(matrix_drm_get_entry)(const JMTX_NAME_TYPED(matrix_drm) * mtx, JMTX_INDEX_T i,
                                                    JMTX_INDEX_T j)
{
    if (mtx->permutations)
    {
        const JMTX_SCALAR_T *const ptr = mtx->values + mtx->permutations[i] * mtx->base.cols;
        return ptr[j];
    }
    const JMTX_SCALAR_T *const ptr = mtx->values + i * mtx->base.cols;
    return ptr[j];
}

/**
 * Returns a pointer to an entry in the matrix.
 * @param mtx pointer to the memory where the matrix is stored
 * @param i row index
 * @param j column index
 * @return pointer to an entry
 */
JMTX_SCALAR_T *JMTX_NAME_TYPED(matrix_drm_entry_ptr)(const JMTX_NAME_TYPED(matrix_drm) * mtx, JMTX_INDEX_T i,
                                                     JMTX_INDEX_T j)
{
    if (mtx->permutations)
    {
        JMTX_SCALAR_T *const ptr = mtx->values + mtx->permutations[i] * mtx->base.cols;
        return ptr + j;
    }
    JMTX_SCALAR_T *const ptr = mtx->values + i * mtx->base.cols;
    return ptr + j;
}

/**
 * Zeros all entries within a matrix, but does not remove them in case they need to be reused
 * @param mtx matrix to zero
 */
void JMTX_NAME_TYPED(matrix_drm_zero_all_entries)(const JMTX_NAME_TYPED(matrix_drm) * mtx)
{
    memset(mtx->values, 0, sizeof(*mtx->values) * mtx->base.cols * mtx->base.rows);
}

/**
 * Similar to JMTX_NAME_TYPED(matrix_drm_set_all_entries, but slower, since it can not use memset. On the other hand, it
 * allows for the value to be other than 0
 * @param mtx matrix to set
 * @param x value to which to set all entries to
 */
void JMTX_NAME_TYPED(matrix_drm_set_all_entries)(const JMTX_NAME_TYPED(matrix_drm) * mtx, JMTX_SCALAR_T x)
{
    for (JMTX_INDEX_T i = 0; i < mtx->base.cols * mtx->base.rows; ++i)
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
JMTX_INDEX_T JMTX_NAME_TYPED(matrix_drm_get_col)(const JMTX_NAME_TYPED(matrix_drm) * mtx, JMTX_INDEX_T col,
                                                 JMTX_SCALAR_T values[])
{
    if (mtx->permutations)
    {
        for (JMTX_INDEX_T i = 0; i < mtx->base.rows; ++i)
        {
            values[i] = mtx->values[mtx->permutations[i] * mtx->base.cols + col];
        }
    }
    else
    {
        for (JMTX_INDEX_T i = 0; i < mtx->base.rows; ++i)
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
jmtx_result JMTX_NAME_TYPED(matrix_drm_transpose)(const JMTX_NAME_TYPED(matrix_drm) * mtx,
                                                  JMTX_NAME_TYPED(matrix_drm) * *p_out,
                                                  const jmtx_allocator_callbacks *allocator_callbacks)
{
    JMTX_NAME_TYPED(matrix_drm) * new;
    const JMTX_INDEX_T new_cols = mtx->base.rows, new_rows = mtx->base.cols;
    jmtx_result res = JMTX_NAME_TYPED(matrix_drm_new)(&new, new_rows, new_cols, NULL, allocator_callbacks);
    if (res != JMTX_RESULT_SUCCESS)
    {
        return res;
    }

    if (mtx->permutations)
    {
        for (JMTX_INDEX_T i = 0; i < mtx->base.rows; ++i)
        {
            for (JMTX_INDEX_T j = 0; j < mtx->base.cols; ++j)
            {
                new->values[j * new->base.cols + i] = mtx->values[mtx->permutations[i] * mtx->base.cols + j];
            }
        }
    }
    else
    {
        for (JMTX_INDEX_T i = 0; i < mtx->base.rows; ++i)
        {
            for (JMTX_INDEX_T j = 0; j < mtx->base.cols; ++j)
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
jmtx_result JMTX_NAME_TYPED(matrix_drm_transpose_inplace)(JMTX_NAME_TYPED(matrix_drm) * mtx, JMTX_SCALAR_T *aux_row)
{
    const JMTX_INDEX_T new_cols = mtx->base.rows, new_rows = mtx->base.cols;

    if (mtx->permutations)
    {
        //  Commit to permutations now
        JMTX_NAME_TYPED(matrix_drm_commit_permutations2)(mtx, aux_row);
    }
    assert(mtx->permutations == NULL);
    assert(mtx->rperm == NULL);
    for (JMTX_INDEX_T i = 0; i < (mtx->base.rows + 1) / 2; ++i)
    {
        for (JMTX_INDEX_T j = 0; j < mtx->base.cols; ++j)
        {
            const JMTX_SCALAR_T tmp = mtx->values[j * mtx->base.cols + i];
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
jmtx_result JMTX_NAME_TYPED(matrix_drm_copy)(const JMTX_NAME_TYPED(matrix_drm) * mtx,
                                             JMTX_NAME_TYPED(matrix_drm) * *p_out,
                                             const jmtx_allocator_callbacks *allocator_callbacks)
{
    JMTX_NAME_TYPED(matrix_drm) * new;
    jmtx_result res = JMTX_NAME_TYPED(matrix_drm_new)(&new, mtx->base.rows, mtx->base.cols, NULL, allocator_callbacks);
    if (res != JMTX_RESULT_SUCCESS)
    {
        return res;
    }

    if (mtx->permutations)
    {
        for (JMTX_INDEX_T i = 0; i < mtx->base.rows; ++i)
        {
            for (JMTX_INDEX_T j = 0; j < mtx->base.cols; ++j)
            {
                new->values[j * new->base.cols + i] = mtx->values[mtx->permutations[i] * mtx->base.cols + j];
            }
        }
    }
    else
    {
        for (JMTX_INDEX_T i = 0; i < mtx->base.rows; ++i)
        {
            for (JMTX_INDEX_T j = 0; j < mtx->base.cols; ++j)
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
JMTX_SCALAR_T JMTX_NAME_TYPED(matrix_drm_vector_multiply_row)(const JMTX_NAME_TYPED(matrix_drm) * mtx,
                                                              const JMTX_SCALAR_T *x, JMTX_INDEX_T i)
{
    const JMTX_SCALAR_T *ptr;
    if (mtx->permutations)
    {
        ptr = mtx->values + mtx->base.cols * mtx->permutations[i];
    }
    else
    {
        ptr = mtx->values + mtx->base.cols * i;
    }
    JMTX_SCALAR_T v = 0;
    for (JMTX_INDEX_T j = 0; j < mtx->base.cols; ++j)
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
void JMTX_NAME_TYPED(matrix_drm_get_bandwidths)(const JMTX_NAME_TYPED(matrix_drm) * mtx, JMTX_INDEX_T *ubw,
                                                JMTX_INDEX_T *lbw)
{
    JMTX_INDEX_T l = 0, u = 0;
    //  Quick check for extreme case
    if (mtx->values[mtx->base.cols - 1] != 0 && mtx->values[mtx->base.cols * (mtx->base.rows - 1)] != 0)
    {
        *ubw = mtx->base.cols - 1;
        *lbw = mtx->base.rows - 1;
        return;
    }

    for (JMTX_INDEX_T i = 0; i < mtx->base.rows; ++i)
    {
        const JMTX_SCALAR_T *const row_ptr =
            mtx->values + mtx->base.cols * (mtx->permutations ? mtx->permutations[i] : i);
        //  Check if ubw is greater
        for (JMTX_INDEX_T j = i + u; j < mtx->base.cols; ++j)
        {
            if (row_ptr[j] != 0)
            {
                u = j - i;
            }
        }
        //  Check if lbw is greater
        for (JMTX_INDEX_T j = i - l; j > 0; --j)
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
jmtx_result JMTX_NAME_TYPED(matrix_drm_swap_rows)(JMTX_NAME_TYPED(matrix_drm) * mtx, JMTX_INDEX_T row1,
                                                  JMTX_INDEX_T row2)
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
        for (JMTX_INDEX_T i = 0; i < mtx->base.rows; ++i)
        {
            mtx->permutations[i] = i;
            mtx->rperm[i] = i;
        }
    }

    //  Swap two rows in permutation list
    const JMTX_INDEX_T tmp = mtx->permutations[row1];
    mtx->permutations[row1] = mtx->permutations[row2];
    mtx->permutations[row2] = tmp;
    //  Correct reverse permutation lists
    const JMTX_INDEX_T rtmp = mtx->rperm[mtx->permutations[row1]];
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
jmtx_result JMTX_NAME_TYPED(matrix_drm_set_permutation)(JMTX_NAME_TYPED(matrix_drm) * mtx, const JMTX_INDEX_T *perm)
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
    for (JMTX_INDEX_T i = 0; i < mtx->base.rows; ++i)
    {
        mtx->permutations[i] = perm[i];
    }
    //  Deal with reverse permutations
    JMTX_INDEX_T first = 0, last = mtx->base.rows - 1;
    for (JMTX_INDEX_T i = 0; i < mtx->base.rows; ++i)
    {
        //  Find j such that perm[j] == i
        for (JMTX_INDEX_T j = first; j <= last; ++j)
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
void JMTX_NAME_TYPED(matrix_drm_commit_permutations)(JMTX_NAME_TYPED(matrix_drm) * mtx)
{
    JMTX_INDEX_T pos;
    const JMTX_INDEX_T n = mtx->base.cols;
    for (pos = 0; pos < mtx->base.rows; ++pos)
    {
        JMTX_INDEX_T src = pos;
        JMTX_INDEX_T dst = mtx->rperm[pos];
        if (src == dst)
        {
            //  Already there
            continue;
        }
        if (mtx->rperm[dst] == src)
        {
            assert(mtx->permutations[src] == dst);
            //  Simple swap
            JMTX_SCALAR_T *const p1 = mtx->values + n * dst;
            JMTX_SCALAR_T *const p2 = mtx->values + n * src;
            for (JMTX_INDEX_T i = 0; i < n; ++i)
            {
                const JMTX_SCALAR_T tmp = p1[i];
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
        JMTX_INDEX_T v = 0;
        while (mtx->rperm[src] != pos)
        {
            JMTX_INDEX_T old_dst = dst;
            dst = mtx->rperm[dst];
            src = old_dst;
            v += 1;
        }

        for (JMTX_INDEX_T i = 0; i < n; ++i)
        {
            dst = mtx->rperm[pos];
            src = pos;
            JMTX_SCALAR_T aux = mtx->values[n * src + i];
            while (mtx->rperm[src] != pos)
            {
                const JMTX_SCALAR_T tmp1 = aux;
                aux = mtx->values[n * dst + i];
                mtx->values[n * dst + i] = tmp1;
                JMTX_INDEX_T old_dst = dst;
                dst = mtx->rperm[dst];
                src = old_dst;
            }
            mtx->values[n * pos + i] = aux;
        }
        dst = mtx->rperm[pos];
        for (JMTX_INDEX_T i = 0; i < v; ++i)
        {
            const JMTX_INDEX_T old_dst = dst;
            dst = mtx->rperm[dst];
            mtx->rperm[old_dst] = old_dst;
        }
        mtx->rperm[pos] = pos;
    }
#ifndef NDEBUG
    for (JMTX_INDEX_T i = 0; i < mtx->base.rows; ++i)
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
void JMTX_NAME_TYPED(matrix_drm_commit_permutations2)(JMTX_NAME_TYPED(matrix_drm) * mtx, JMTX_SCALAR_T *aux_row)
{
    JMTX_INDEX_T pos;
    const JMTX_INDEX_T n = mtx->base.cols;
    for (pos = 0; pos < mtx->base.rows; ++pos)
    {
        JMTX_INDEX_T src = pos;
        JMTX_INDEX_T dst = mtx->rperm[pos];
        if (src == dst)
        {
            //  Already there
            continue;
        }
        if (mtx->rperm[dst] == src)
        {
            assert(mtx->permutations[src] == dst);
            //  Simple swap
            JMTX_SCALAR_T *const p1 = mtx->values + n * dst;
            JMTX_SCALAR_T *const p2 = mtx->values + n * src;
            for (JMTX_INDEX_T i = 0; i < n; ++i)
            {
                const JMTX_SCALAR_T tmp = p1[i];
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
        JMTX_INDEX_T v = 0;
        while (mtx->rperm[src] != pos)
        {
            for (JMTX_INDEX_T i = 0; i < n; ++i)
            {
                const JMTX_SCALAR_T tmp = aux_row[i];
                aux_row[i] = mtx->values[n * dst + i];
                mtx->values[n * dst + i] = tmp;
            }
            JMTX_INDEX_T old_dst = dst;
            dst = mtx->rperm[dst];
            src = old_dst;
            v += 1;
        }
        memcpy(mtx->values + n * pos, aux_row, sizeof(*aux_row) * n);
        dst = mtx->rperm[pos];
        for (JMTX_INDEX_T i = 0; i < v; ++i)
        {
            const JMTX_INDEX_T old_dst = dst;
            dst = mtx->rperm[dst];
            mtx->rperm[old_dst] = old_dst;
        }
        mtx->rperm[pos] = pos;
    }
#ifndef NDEBUG
    for (JMTX_INDEX_T i = 0; i < mtx->base.rows; ++i)
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
void JMTX_NAME_TYPED(matrix_drm_givens_rotation_left)(JMTX_NAME_TYPED(matrix_drm) * mtx, unsigned r1, unsigned r2,
                                                      const JMTX_SCALAR_T ct, const JMTX_SCALAR_T st)
{
    const unsigned n = mtx->base.cols;
    if (mtx->permutations)
    {
        r1 = mtx->permutations[r1];
        r2 = mtx->permutations[r2];
    }
    JMTX_SCALAR_T *const row1 = mtx->values + n * r1;
    JMTX_SCALAR_T *const row2 = mtx->values + n * r2;

    for (unsigned j = 0; j < n; ++j)
    {
        const JMTX_SCALAR_T a = row1[j];
        const JMTX_SCALAR_T b = row2[j];
        row1[j] = ct * a - st * b;
        row2[j] = st * a + ct * b;
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
 */
void JMTX_NAME_TYPED(matrix_drm_givens_rotation_right)(JMTX_NAME_TYPED(matrix_drm) * mtx, unsigned c1, unsigned c2,
                                                       JMTX_SCALAR_T ct, JMTX_SCALAR_T st)
{
    const unsigned n = mtx->base.rows;

    JMTX_SCALAR_T *const col1 = mtx->values + c1;
    JMTX_SCALAR_T *const col2 = mtx->values + c2;

    if (!mtx->permutations)
    {
        for (unsigned i = 0; i < n; ++i)
        {
            const JMTX_SCALAR_T a = col1[i * mtx->base.cols];
            const JMTX_SCALAR_T b = col2[i * mtx->base.cols];
            col1[i * mtx->base.cols] = ct * a - st * b;
            col2[i * mtx->base.cols] = st * a + ct * b;
        }
    }
    else
    {
        for (unsigned i = 0; i < n; ++i)
        {
            const JMTX_SCALAR_T a = col1[mtx->permutations[i] * mtx->base.cols];
            const JMTX_SCALAR_T b = col2[mtx->permutations[i] * mtx->base.cols];
            col1[mtx->rperm[i] * mtx->base.cols] = ct * a - st * b;
            col2[mtx->rperm[i] * mtx->base.cols] = st * a + ct * b;
        }
    }
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
jmtx_result JMTX_NAME_TYPED(matrix_drm_multiply_matrix)(JMTX_NAME_TYPED(matrix_drm) * a,
                                                        JMTX_NAME_TYPED(matrix_drm) * b,
                                                        JMTX_NAME_TYPED(matrix_drm) * out)
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
            const JMTX_SCALAR_T *restrict p_row = a->values + i * n;
            for (unsigned j = 0; j < new_cols; ++j)
            {
                JMTX_SCALAR_T v = 0;
                const JMTX_SCALAR_T *restrict p_col = b->values + j;
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
            const JMTX_SCALAR_T *restrict p_row = a->values + n * a->permutations[i];
            for (unsigned j = 0; j < new_cols; ++j)
            {
                JMTX_SCALAR_T v = 0;
                const JMTX_SCALAR_T *restrict p_col = b->values + j;
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
            const JMTX_SCALAR_T *restrict p_row = a->values + n * a->permutations[i];
            for (unsigned j = 0; j < new_cols; ++j)
            {
                JMTX_SCALAR_T v = 0;
                const JMTX_SCALAR_T *restrict p_col = b->values + j;
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
            const JMTX_SCALAR_T *restrict p_row = a->values + n * i;
            for (unsigned j = 0; j < new_cols; ++j)
            {
                JMTX_SCALAR_T v = 0;
                const JMTX_SCALAR_T *restrict p_col = b->values + j;
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
void JMTX_NAME_TYPED(matrix_drm_shift_diagonal)(JMTX_NAME_TYPED(matrix_drm) * mtx, JMTX_SCALAR_T v)
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

JMTX_NAME_TYPED(matrix_drm)
JMTX_NAME_TYPED(matrix_drm_from_data)(const unsigned rows, const unsigned cols,
                                      JMTX_SCALAR_T values[JMTX_ARRAY_ATTRIB(static rows * cols)])
{
    return (JMTX_NAME_TYPED(matrix_drm)){.base = {.type = JMTXD_TYPE_DRM,
                                                  .rows = rows,
                                                  .cols = cols,
                                                  .allocator_callbacks = {.state = NULL, .alloc = NULL, .free = NULL}},
                                         .permutations = NULL,
                                         .rperm = NULL,
                                         .values = values};
}
