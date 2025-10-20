// Automatically generated from source/float/matrices/sparse_diagonal_compressed.c on Fri Dec  1 06:43:01 2023
//
// Created by jan on 27.11.2023.
//

#include "sparse_diagonal_compressed.h"
#include <assert.h>
#include <math.h>

enum
{
    MINIMUM_RESERVED_DIAGONALS = 8
};

static inline jmtx_result JMTX_NAME_TYPED(matrix_cds_diagonal_array_init)(
    const jmtx_allocator_callbacks *allocator_callbacks, JMTX_NAME_TYPED(matrix_cds_diagonal_array) * ptr,
    JMTX_INDEX_T capacity)
{
    JMTX_INDEX_T *const idx = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*idx) * capacity);
    if (!idx)
    {
        return JMTX_RESULT_BAD_ALLOC;
    }
    JMTX_SCALAR_T **const diags = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*diags) * capacity);
    if (!diags)
    {
        allocator_callbacks->free(allocator_callbacks->state, idx);
        return JMTX_RESULT_BAD_ALLOC;
    }
    ptr->capacity = capacity;
    ptr->count = 0;
    ptr->indices = idx;
    ptr->diagonals = diags;
    return JMTX_RESULT_SUCCESS;
}

static inline jmtx_result JMTX_NAME_TYPED(matrix_cds_diagonal_array_insert)(
    const jmtx_allocator_callbacks *allocator_callbacks, JMTX_NAME_TYPED(matrix_cds_diagonal_array) * this,
    JMTX_INDEX_T idx, JMTX_SCALAR_T *dia)
{
    assert(this->capacity != 0);
    if (this->capacity == this->count)
    {
        const JMTX_INDEX_T new_capacity = this->capacity << 1;
        JMTX_INDEX_T *const new_idx =
            allocator_callbacks->realloc(allocator_callbacks->state, this->indices, sizeof(*new_idx) * new_capacity);
        if (!idx)
        {
            return JMTX_RESULT_BAD_ALLOC;
        }
        this->indices = new_idx;
        JMTX_SCALAR_T **const diags =
            allocator_callbacks->realloc(allocator_callbacks->state, this->diagonals, sizeof(*diags) * new_capacity);
        if (!diags)
        {
            return JMTX_RESULT_BAD_ALLOC;
        }
        this->diagonals = diags;
        this->capacity = new_capacity;
    }

    JMTX_INDEX_T *const indices = this->indices;
    JMTX_SCALAR_T **const diagonals = this->diagonals;

    JMTX_INDEX_T pos;
    for (pos = 0; pos < this->count; ++pos)
    {
        if (indices[pos] > idx)
        {
            break;
        }
    }
    if (pos != this->count)
    {
        memmove(indices + pos + 1, indices + pos, sizeof(*indices) * (this->count - pos));
        memmove(diagonals + pos + 1, diagonals + pos, sizeof(*diagonals) * (this->count - pos));
    }
    diagonals[pos] = dia;
    indices[pos] = idx;

    this->count += 1;
    return JMTX_RESULT_SUCCESS;
}

static inline JMTX_SCALAR_T *JMTX_NAME_TYPED(matrix_cds_diagonal_array_get_ptr)(
    const JMTX_NAME_TYPED(matrix_cds_diagonal_array) * this, JMTX_INDEX_T idx)
{
    for (JMTX_INDEX_T i = 0; i < this->count; ++i)
    {
        if (this->indices[i] == idx)
        {
            return this->diagonals[i];
        }
    }
    return NULL;
}

static inline void JMTX_NAME_TYPED(matrix_cds_diagonal_array_free)(const jmtx_allocator_callbacks *allocator_callbacks,
                                                                   JMTX_NAME_TYPED(matrix_cds_diagonal_array) * ptr)
{
    for (JMTX_INDEX_T i = 0; i < ptr->count; ++i)
    {
        allocator_callbacks->free(allocator_callbacks->state, ptr->diagonals[i]);
    }
    allocator_callbacks->free(allocator_callbacks->state, ptr->diagonals);
    allocator_callbacks->free(allocator_callbacks->state, ptr->indices);
#ifndef NDEBUG
    ptr->diagonals = (void *)0xCCCCCCCCCCCCCCCC;
    ptr->indices = (void *)0xCCCCCCCCCCCCCCCC;
    ptr->count = 0xCCCCCCCC;
    ptr->capacity = 0xCCCCCCCC;
#endif
}

static inline JMTX_FAST_INT_T cds_subdiagonal_length(const JMTX_FAST_INT_T cols, const JMTX_FAST_INT_T rows,
                                                     const JMTX_FAST_INT_T dia)
{
    assert(rows > dia);
    if (rows - dia > cols)
    {
        return cols;
    }
    return rows - dia;
}
static inline JMTX_FAST_INT_T cds_superdiagonal_length(const JMTX_FAST_INT_T cols, const JMTX_FAST_INT_T rows,
                                                       const JMTX_FAST_INT_T dia)
{
    assert(cols > dia);
    if (cols - dia > rows)
    {
        return rows;
    }
    return cols - dia;
}

/**
 * Initializes a new Compressed Diagonal Sparse matrix
 * @param p_mtx address that receives the pointer to the matrix
 * @param rows number of rows of the sparse matrix
 * @param cols number of columns of the sparse matrix
 * @param reserved_entries how many entries should the space be reserved for in the matrix initially
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result JMTX_NAME_TYPED(matrix_cds_new)(JMTX_NAME_TYPED(matrix_cds) * *p_mtx, JMTX_INDEX_T rows, JMTX_INDEX_T cols,
                                            JMTX_INDEX_T n_diagonals,
                                            const int32_t p_dia_idx[JMTX_ARRAY_ATTRIB(static n_diagonals)],
                                            const jmtx_allocator_callbacks *allocator_callbacks)
{
    const int32_t backup_diag[1] = {0};
    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }
    if (n_diagonals == 0)
    {
        //  Always prepare the main diagonal
        p_dia_idx = backup_diag;
        n_diagonals = 1;
    }

    jmtx_result mtx_res;

    JMTX_NAME_TYPED(matrix_cds) *const mtx = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*mtx));
    if (!mtx)
    {
        return JMTX_RESULT_BAD_ALLOC;
    }

    JMTX_INDEX_T n_sub = 0, n_sup = 0;
    for (JMTX_INDEX_T i = 0; i < n_diagonals; ++i)
    {
        if (p_dia_idx[i] < 0)
        {
            n_sub += 1;
        }
        else if (p_dia_idx[i] > 0)
        {
            n_sup += 1;
        }
    }

    mtx_res = JMTX_NAME_TYPED(matrix_cds_diagonal_array_init)(
        allocator_callbacks, &mtx->sub_diagonals,
        n_sub < MINIMUM_RESERVED_DIAGONALS ? MINIMUM_RESERVED_DIAGONALS : n_sub);
    if (mtx_res != JMTX_RESULT_SUCCESS)
    {
        goto failed_sub;
    }

    mtx_res = JMTX_NAME_TYPED(matrix_cds_diagonal_array_init)(
        allocator_callbacks, &mtx->super_diagonals,
        n_sup < MINIMUM_RESERVED_DIAGONALS ? MINIMUM_RESERVED_DIAGONALS : n_sup);
    if (mtx_res != JMTX_RESULT_SUCCESS)
    {
        goto failed_sup;
    }

    for (JMTX_INDEX_T i = 0; i < n_diagonals; ++i)
    {
        if (p_dia_idx[i] < 0)
        {
            JMTX_SCALAR_T *const ptr = allocator_callbacks->alloc(
                allocator_callbacks->state, sizeof(*ptr) * cds_subdiagonal_length(cols, rows, -p_dia_idx[i]));
            if (ptr == NULL)
            {
                mtx_res = JMTX_RESULT_BAD_ALLOC;
                goto failed_adding;
            }
            JMTX_NAME_TYPED(matrix_cds_diagonal_array_insert)(allocator_callbacks, &mtx->sub_diagonals, -p_dia_idx[i],
                                                              ptr);
        }
        else if (p_dia_idx[i] > 0)
        {
            JMTX_SCALAR_T *const ptr = allocator_callbacks->alloc(
                allocator_callbacks->state, sizeof(*ptr) * cds_superdiagonal_length(cols, rows, p_dia_idx[i]));
            if (ptr == NULL)
            {
                mtx_res = JMTX_RESULT_BAD_ALLOC;
                goto failed_adding;
            }
            JMTX_NAME_TYPED(matrix_cds_diagonal_array_insert)(allocator_callbacks, &mtx->super_diagonals, p_dia_idx[i],
                                                              ptr);
        }
    }
    if (n_diagonals > n_sub + n_sup)
    {
        assert(n_sup + n_sub + 1 == n_diagonals);
        const JMTX_FAST_INT_T len_main = cds_subdiagonal_length(cols, rows, 0);
        assert(len_main == cds_superdiagonal_length(cols, rows, 0));
        JMTX_SCALAR_T *const main_ptr =
            allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*main_ptr) * len_main);
        if (!main_ptr)
        {
            mtx_res = JMTX_RESULT_BAD_ALLOC;
            goto failed_adding;
        }
        mtx->main_diagonal = main_ptr;
    }
    else
    {
        mtx->main_diagonal = 0;
    }

    mtx->base.type = JMTXD_TYPE_CDS;
    mtx->base.cols = cols;
    mtx->base.rows = rows;
    mtx->base.allocator_callbacks = *allocator_callbacks;

    *p_mtx = mtx;

    return JMTX_RESULT_SUCCESS;
failed_adding:
    JMTX_NAME_TYPED(matrix_cds_diagonal_array_free)(allocator_callbacks, &mtx->sub_diagonals);
failed_sup:
    JMTX_NAME_TYPED(matrix_cds_diagonal_array_free)(allocator_callbacks, &mtx->super_diagonals);
failed_sub:
    allocator_callbacks->free(allocator_callbacks->state, mtx);

    return mtx_res;
}

/**
 * Cleans up the cds matrix and frees all of its memory
 * @param mtx pointer to memory where the matrix is stored
 */
void JMTX_NAME_TYPED(matrix_cds_destroy)(JMTX_NAME_TYPED(matrix_cds) * mtx)
{
    const jmtx_allocator_callbacks allocator_callbacks = mtx->base.allocator_callbacks;
    allocator_callbacks.free(allocator_callbacks.state, mtx->main_diagonal);
    JMTX_NAME_TYPED(matrix_cds_diagonal_array_free)(&allocator_callbacks, &mtx->sub_diagonals);
    JMTX_NAME_TYPED(matrix_cds_diagonal_array_free)(&allocator_callbacks, &mtx->super_diagonals);
    allocator_callbacks.free(allocator_callbacks.state, mtx);
}

/**
 * Version of JMTX_NAME_TYPED(matrix_cds_set_row which does not touch the count of entries after the current row. This
 * is useful when building a new matrix, as it avoids unnecessary setting and resetting of these entries. Must be called
 * for each row in order to ensure that the matrix is properly built. Makes no checks on the input parameters
 * @param mtx pointer to the memory where the matrix is stored
 * @param dia_idx index of the diagonal to set
 * @param values values of non-zero values
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result JMTX_NAME_TYPED(matrix_cds_set_diagonal_full)(JMTX_NAME_TYPED(matrix_cds) * mtx, int32_t dia_idx,
                                                          const JMTX_SCALAR_T values[JMTX_ARRAY_ATTRIB(const)])
{
    JMTX_INDEX_T len;
    JMTX_SCALAR_T *const ptr = JMTX_NAME_TYPED(matrix_cds_allocate_diagonal)(mtx, dia_idx, &len);
    if (!ptr)
    {
        return JMTX_RESULT_BAD_ALLOC;
    }
    memcpy(ptr, values, sizeof(*ptr) * len);
    return JMTX_RESULT_SUCCESS;
}

/**
 * Sets the specified number of entries in the diagonal, starting at the given offset for up to n entries at most
 * @param mtx pointer to the memory where the matrix is stored
 * @param dia_idx index of the diagonal to set
 * @param offset offset from the start of the diagonal to start at
 * @param n maximum number of entries to set
 * @param p_count (optional) pointer which receives the number of entries that were actually set (may be less than n if
 * it would go out of bounds)
 * @param values values of non-zero values
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result JMTX_NAME_TYPED(matrix_cds_set_diagonal_part)(JMTX_NAME_TYPED(matrix_cds) * mtx, int32_t dia_idx,
                                                          JMTX_INDEX_T offset, JMTX_INDEX_T n, JMTX_INDEX_T *p_count,
                                                          const JMTX_SCALAR_T values[JMTX_ARRAY_ATTRIB(static n)])
{
    JMTX_INDEX_T len;
    JMTX_SCALAR_T *const ptr = JMTX_NAME_TYPED(matrix_cds_allocate_diagonal)(mtx, dia_idx, &len);
    if (!ptr)
    {
        return JMTX_RESULT_BAD_ALLOC;
    }
    if (len > offset + n)
    {
        len = offset + n;
    }
    else if (offset > len)
    {
        if (p_count)
        {
            *p_count = 0;
        }
        return JMTX_RESULT_SUCCESS;
    }
    JMTX_SCALAR_T *const p = ptr + offset;
    for (JMTX_FAST_INT_T i = 0; i < len - offset; ++i)
    {
        p[i] = values[i];
    }
    if (p_count)
    {
        *p_count = len - offset;
    }
    return JMTX_RESULT_SUCCESS;
}

/**
 * Allocates a new diagonal if one does not already exits, otherwise it returns the pointer to the existing one.
 * @param mtx matrix which to allocate the diagonal for
 * @param dia offset of the diagonal from the main diagonal
 * @param p_size pointer which receives the number of elements in the allocated diagonal. May be NULL
 * @return pointer to the newly allocated diagonal, or NULL in case it failed to do so
 */
JMTX_SCALAR_T *JMTX_NAME_TYPED(matrix_cds_allocate_diagonal)(JMTX_NAME_TYPED(matrix_cds) * mtx, int32_t dia,
                                                             JMTX_INDEX_T *p_size)
{
    JMTX_SCALAR_T *ptr;
    JMTX_FAST_INT_T len;
    JMTX_FAST_INT_T idx;
    if (dia == 0)
    {
        len = cds_superdiagonal_length(mtx->base.cols, mtx->base.rows, 0);
        //  Main diagonal
        if (mtx->main_diagonal == NULL)
        {
            mtx->main_diagonal = mtx->base.allocator_callbacks.alloc(mtx->base.allocator_callbacks.state,
                                                                     sizeof(*mtx->main_diagonal) * len);
            if (!mtx->main_diagonal)
            {
                return NULL;
            }
        }
        ptr = mtx->main_diagonal;
    }
    else if (dia < 0)
    {
        idx = -dia;
        len = cds_subdiagonal_length(mtx->base.cols, mtx->base.rows, idx);
        ptr = JMTX_NAME_TYPED(matrix_cds_diagonal_array_get_ptr)(&mtx->sub_diagonals, idx);
        if (!ptr)
        {
            ptr = mtx->base.allocator_callbacks.alloc(mtx->base.allocator_callbacks.state, sizeof(*ptr) * len);
            if (!ptr)
            {
                return NULL;
            }
            const jmtx_result res = JMTX_NAME_TYPED(matrix_cds_diagonal_array_insert)(&mtx->base.allocator_callbacks,
                                                                                      &mtx->sub_diagonals, idx, ptr);
            if (res != JMTX_RESULT_SUCCESS)
            {
                mtx->base.allocator_callbacks.free(mtx->base.allocator_callbacks.state, ptr);
                return NULL;
            }
        }
    }
    else // if (dia_idx > 0)
    {
        idx = dia;
        len = cds_superdiagonal_length(mtx->base.cols, mtx->base.rows, idx);
        ptr = JMTX_NAME_TYPED(matrix_cds_diagonal_array_get_ptr)(&mtx->super_diagonals, idx);
        if (!ptr)
        {
            ptr = mtx->base.allocator_callbacks.alloc(mtx->base.allocator_callbacks.state, sizeof(*ptr) * len);
            if (!ptr)
            {
                return NULL;
            }
            const jmtx_result res = JMTX_NAME_TYPED(matrix_cds_diagonal_array_insert)(&mtx->base.allocator_callbacks,
                                                                                      &mtx->super_diagonals, idx, ptr);
            if (res != JMTX_RESULT_SUCCESS)
            {
                mtx->base.allocator_callbacks.free(mtx->base.allocator_callbacks.state, ptr);
                return NULL;
            }
        }
    }
    if (p_size)
    {
        *p_size = len;
    }
    return ptr;
}

/**
 * Allocates a new diagonal if one does not already exits, otherwise it returns the pointer to the existing one.
 * Sets the diagonal to zero if it was not allocated before.
 * @param mtx matrix which to allocate the diagonal for
 * @param dia offset of the diagonal from the main diagonal
 * @param p_size pointer which receives the number of elements in the allocated diagonal. May be NULL
 * @return pointer to the newly allocated diagonal, or NULL in case it failed to do so
 */
JMTX_SCALAR_T *JMTX_NAME_TYPED(matrix_cds_allocate_zero_diagonal)(JMTX_NAME_TYPED(matrix_cds) * mtx, int32_t dia,
                                                                  JMTX_INDEX_T *p_size)
{
    JMTX_SCALAR_T *ptr;
    JMTX_FAST_INT_T len;
    JMTX_FAST_INT_T idx;
    if (dia == 0)
    {
        len = cds_superdiagonal_length(mtx->base.cols, mtx->base.rows, 0);
        //  Main diagonal
        if (mtx->main_diagonal == NULL)
        {
            mtx->main_diagonal = mtx->base.allocator_callbacks.alloc(mtx->base.allocator_callbacks.state,
                                                                     sizeof(*mtx->main_diagonal) * len);
            if (!mtx->main_diagonal)
            {
                return NULL;
            }
            memset(mtx->main_diagonal, 0, sizeof(*mtx->main_diagonal) * len);
        }
        ptr = mtx->main_diagonal;
    }
    else if (dia < 0)
    {
        idx = -dia;
        len = cds_subdiagonal_length(mtx->base.cols, mtx->base.rows, idx);
        ptr = JMTX_NAME_TYPED(matrix_cds_diagonal_array_get_ptr)(&mtx->sub_diagonals, idx);
        if (!ptr)
        {
            ptr = mtx->base.allocator_callbacks.alloc(mtx->base.allocator_callbacks.state, sizeof(*ptr) * len);
            if (!ptr)
            {
                return NULL;
            }
            const jmtx_result res = JMTX_NAME_TYPED(matrix_cds_diagonal_array_insert)(&mtx->base.allocator_callbacks,
                                                                                      &mtx->sub_diagonals, idx, ptr);
            if (res != JMTX_RESULT_SUCCESS)
            {
                mtx->base.allocator_callbacks.free(mtx->base.allocator_callbacks.state, ptr);
                return NULL;
            }
            memset(ptr, 0, sizeof(*ptr) * len);
        }
    }
    else // if (dia_idx > 0)
    {
        idx = dia;
        len = cds_superdiagonal_length(mtx->base.cols, mtx->base.rows, idx);
        ptr = JMTX_NAME_TYPED(matrix_cds_diagonal_array_get_ptr)(&mtx->super_diagonals, idx);
        if (!ptr)
        {
            ptr = mtx->base.allocator_callbacks.alloc(mtx->base.allocator_callbacks.state, sizeof(*ptr) * len);
            if (!ptr)
            {
                return NULL;
            }
            const jmtx_result res = JMTX_NAME_TYPED(matrix_cds_diagonal_array_insert)(&mtx->base.allocator_callbacks,
                                                                                      &mtx->super_diagonals, idx, ptr);
            if (res != JMTX_RESULT_SUCCESS)
            {
                mtx->base.allocator_callbacks.free(mtx->base.allocator_callbacks.state, ptr);
                return NULL;
            }
            memset(ptr, 0, sizeof(*ptr) * len);
        }
    }
    if (p_size)
    {
        *p_size = len;
    }
    return ptr;
}

/**
 * Returns the pointer to the diagonal if one already exits, otherwise it returns NULL
 * @param mtx matrix which to allocate the diagonal for
 * @param dia offset of the diagonal from the main diagonal
 * @return pointer to the diagonal, or NULL in case it does not exist
 */
JMTX_SCALAR_T *JMTX_NAME_TYPED(matrix_cds_get_diagonal)(const JMTX_NAME_TYPED(matrix_cds) * mtx, int32_t dia)
{
    if (dia == 0)
    {
        return mtx->main_diagonal;
    }
    else if (dia < 0)
    {
        return JMTX_NAME_TYPED(matrix_cds_diagonal_array_get_ptr)(&mtx->sub_diagonals, -dia);
    }
    //  else if (dia > 0)
    return JMTX_NAME_TYPED(matrix_cds_diagonal_array_get_ptr)(&mtx->super_diagonals, dia);
}

/**
 * Returns the number of entries in the row of the matrix
 * @param mtx pointer to the memory where the matrix is stored
 * @param row row index of the matrix to look at
 * @return number of entries in the row
 */
JMTX_INDEX_T JMTX_NAME_TYPED(matrix_cds_entries_in_row)(const JMTX_NAME_TYPED(matrix_cds) * mtx, JMTX_INDEX_T row)
{
    JMTX_FAST_INT_T k = 0;
    for (JMTX_FAST_INT_T i = 0; i < mtx->sub_diagonals.count; ++i)
    {
        //  Check the subdiagonals
        const JMTX_FAST_INT_T d = mtx->sub_diagonals.indices[mtx->sub_diagonals.count - 1 - i];
        if (row < d || row >= d + cds_subdiagonal_length(mtx->base.cols, mtx->base.rows, d))
        {
            continue;
        }

        k += 1;
    }
    if (mtx->main_diagonal && row < mtx->base.cols)
    {
        k += 1;
    }
    for (JMTX_FAST_INT_T i = 0; i < mtx->super_diagonals.count; ++i)
    {
        //  Check the superdiagonals
        const JMTX_FAST_INT_T idx = mtx->super_diagonals.indices[i];
        if (row >= mtx->base.cols - idx)
        {
            break;
        }

        k += 1;
    }

    return k;
}

/**
 * Returns the values of entries in the matrix, along with what column of the matrix they were located in
 * @param mtx pointer to the memory where the matrix is stored
 * @param row row index of the matrix to look at
 * @param n number of values in the row to be extracted at most
 * @param p_values a buffer of at least n values which receives the values of the row
 * @param p_cols a buffer of at least n values which receives the column indices of the row
 * @return number of entries that were extracted from the column (may be less than are really in the column if n was too
 * small)
 */
JMTX_INDEX_T JMTX_NAME_TYPED(matrix_cds_get_row)(const JMTX_NAME_TYPED(matrix_cds) * mtx, JMTX_INDEX_T row,
                                                 JMTX_INDEX_T n, JMTX_SCALAR_T p_values[JMTX_ARRAY_ATTRIB(restrict n)],
                                                 JMTX_INDEX_T p_cols[JMTX_ARRAY_ATTRIB(restrict n)])
{
    JMTX_FAST_INT_T k = 0;
    for (JMTX_FAST_INT_T i = 0; i < mtx->sub_diagonals.count && k < n; ++i)
    {
        //  Check the subdiagonals
        const JMTX_FAST_INT_T d = mtx->sub_diagonals.indices[mtx->sub_diagonals.count - 1 - i];
        if (row < d || row >= d + cds_subdiagonal_length(mtx->base.cols, mtx->base.rows, d))
        {
            continue;
        }

        const JMTX_FAST_INT_T pos = row - d;
        p_cols[k] = pos;
        const JMTX_SCALAR_T *const diagonal = mtx->sub_diagonals.diagonals[mtx->sub_diagonals.count - 1 - i];
        p_values[k] = diagonal[pos];
        k += 1;
    }
    if (mtx->main_diagonal && row < mtx->base.cols && k < n)
    {
        p_cols[k] = row;
        p_values[k] = mtx->main_diagonal[row];
        k += 1;
    }
    for (JMTX_FAST_INT_T i = 0; i < mtx->super_diagonals.count && k < n; ++i)
    {
        //  Check the superdiagonals
        const JMTX_FAST_INT_T idx = mtx->super_diagonals.indices[i];
        if (row >= mtx->base.cols - idx)
        {
            break;
        }

        p_cols[k] = idx + row;
        const JMTX_SCALAR_T *const diagonal = mtx->super_diagonals.diagonals[i];
        p_values[k] = diagonal[row];
        k += 1;
    }

    return k;
}

/**
 * Returns the values of entries in the matrix, along with what column of the matrix they were located in
 * @param mtx pointer to the memory where the matrix is stored
 * @param row row index of the matrix to look at
 * @param n number of values to be set in the row
 * @param p_values a buffer of at least n values which contain the values to be set to
 * @param p_cols a buffer of at least n values which receives the column indices of the row
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result JMTX_NAME_TYPED(matrix_cds_set_row)(JMTX_NAME_TYPED(matrix_cds) * mtx, JMTX_INDEX_T row, JMTX_INDEX_T n,
                                                const JMTX_SCALAR_T p_values[JMTX_ARRAY_ATTRIB(restrict static n)],
                                                const JMTX_INDEX_T p_cols[JMTX_ARRAY_ATTRIB(restrict static n)])
{
    JMTX_NAME_TYPED(matrix_cds_zero_row)(mtx, row);
    for (JMTX_FAST_INT_T m = 0; m < n; ++m)
    {
        const int32_t dia = (int32_t)p_cols[m] - (int32_t)row;
        JMTX_SCALAR_T *ptr;
        JMTX_SCALAR_T v;
        ptr = JMTX_NAME_TYPED(matrix_cds_allocate_zero_diagonal)(mtx, dia, NULL);
        if (!ptr)
        {
            return JMTX_RESULT_BAD_ALLOC;
        }
        v = p_values[m];
        if (dia < 0)
        {
            ptr[p_cols[m]] = v;
        }
        else
        {
            ptr[row] = v;
        }
    }

    return JMTX_RESULT_SUCCESS;
}

/**
 * Sets all the entries in the column of the matrix, zeroing non-specified entries and allocating new diagonals as
 * needed
 * @param mtx pointer to the memory where the matrix is stored
 * @param col column index of the matrix to look at
 * @param n number of values in the column to be set
 * @param p_values a buffer of n values of column entries
 * @param p_rows a buffer of n indices of column indices of columns to set
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on allocation failure
 */
jmtx_result JMTX_NAME_TYPED(matrix_cds_set_col)(JMTX_NAME_TYPED(matrix_cds) * mtx, JMTX_INDEX_T col, JMTX_INDEX_T n,
                                                const JMTX_SCALAR_T p_values[JMTX_ARRAY_ATTRIB(restrict static n)],
                                                const JMTX_INDEX_T p_rows[JMTX_ARRAY_ATTRIB(restrict static n)])
{
    JMTX_NAME_TYPED(matrix_cds_zero_col)(mtx, col);
    for (JMTX_FAST_INT_T m = 0; m < n; ++m)
    {
        const int32_t dia = (int32_t)col - (int32_t)p_rows[m];
        JMTX_SCALAR_T *ptr;
        JMTX_SCALAR_T v;
        ptr = JMTX_NAME_TYPED(matrix_cds_allocate_zero_diagonal)(mtx, dia, NULL);
        if (!ptr)
        {
            return JMTX_RESULT_BAD_ALLOC;
        }
        v = p_values[m];
        if (dia <= 0)
        {
            ptr[col] = v;
        }
        else
        {
            ptr[p_rows[m]] = v;
        }
    }

    return JMTX_RESULT_SUCCESS;
}

/**
 * Returns the number of entries in the column of the matrix
 * @param mtx pointer to the memory where the matrix is stored
 * @param col column index of the matrix to look at
 * @return number of entries in the column
 */
JMTX_INDEX_T JMTX_NAME_TYPED(matrix_cds_entries_in_col)(const JMTX_NAME_TYPED(matrix_cds) * mtx, JMTX_INDEX_T col)
{
    JMTX_FAST_INT_T k = 0;
    for (JMTX_FAST_INT_T i = 0; i < mtx->super_diagonals.count; ++i)
    {
        //  Check the superdiagonals
        const JMTX_FAST_INT_T d = mtx->super_diagonals.indices[mtx->super_diagonals.count - 1 - i];
        if (col < d || col >= d + cds_superdiagonal_length(mtx->base.cols, mtx->base.rows, d))
        {
            continue;
        }

        k += 1;
    }
    if (mtx->main_diagonal && col < mtx->base.rows)
    {
        k += 1;
    }
    for (JMTX_FAST_INT_T i = 0; i < mtx->sub_diagonals.count; ++i)
    {
        //  Check the subdiagonals
        const JMTX_FAST_INT_T idx = mtx->sub_diagonals.indices[i];
        if (col >= mtx->base.rows - idx)
        {
            break;
        }

        k += 1;
    }
    return k;
}

/**
 * Returns the values of entries in the matrix, along with what row of the matrix they were located in
 * @param mtx pointer to the memory where the matrix is stored
 * @param col column index of the matrix to look at
 * @param n number of values in the column to be extracted at most
 * @param p_values a buffer of at least n values which receives the values of the column
 * @param p_rows a buffer of at least n values which receives the row indices of the column
 * @return number of entries that were extracted from the column (may be less than are really in the column if n was too
 * small)
 */
JMTX_INDEX_T JMTX_NAME_TYPED(matrix_cds_get_col)(const JMTX_NAME_TYPED(matrix_cds) * mtx, JMTX_INDEX_T col,
                                                 JMTX_INDEX_T n, JMTX_SCALAR_T p_values[JMTX_ARRAY_ATTRIB(restrict n)],
                                                 JMTX_INDEX_T p_rows[JMTX_ARRAY_ATTRIB(restrict n)])
{
    JMTX_FAST_INT_T k = 0;
    for (JMTX_FAST_INT_T i = 0; i < mtx->super_diagonals.count && k < n; ++i)
    {
        //  Check the superdiagonals
        const JMTX_FAST_INT_T d = mtx->super_diagonals.indices[mtx->super_diagonals.count - 1 - i];
        if (col < d || col >= d + cds_superdiagonal_length(mtx->base.cols, mtx->base.rows, d))
        {
            continue;
        }
        const JMTX_FAST_INT_T pos = col - d;
        p_rows[k] = pos;
        const JMTX_SCALAR_T *const diagonal = mtx->super_diagonals.diagonals[mtx->super_diagonals.count - 1 - i];
        p_values[k] = diagonal[pos];
        k += 1;
    }
    if (mtx->main_diagonal && col < mtx->base.rows && k < n)
    {
        p_rows[k] = col;
        p_values[k] = mtx->main_diagonal[col];
        k += 1;
    }
    for (JMTX_FAST_INT_T i = 0; i < mtx->sub_diagonals.count && k < n; ++i)
    {
        //  Check the subdiagonals
        const JMTX_FAST_INT_T idx = mtx->sub_diagonals.indices[i];
        if (col >= mtx->base.rows - idx)
        {
            break;
        }

        p_rows[k] = idx + col;
        const JMTX_SCALAR_T *const diagonal = mtx->sub_diagonals.diagonals[i];
        p_values[k] = diagonal[col];
        k += 1;
    }
    return k;
}

/**
 * Returns the number of entries in the diagonal of the matrix
 * @param mtx pointer to the memory where the matrix is stored
 * @param dia index of the diagonal of the matrix to look at
 * @return number of entries in the diagonal
 */
JMTX_FAST_INT_T JMTX_NAME_TYPED(matrix_cds_entries_in_dia)(JMTX_NAME_TYPED(matrix_cds) * mtx, int32_t dia)
{
    if (dia >= 0)
    {
        return cds_superdiagonal_length(mtx->base.cols, mtx->base.rows, dia);
    }
    return cds_subdiagonal_length(mtx->base.cols, mtx->base.rows, -dia);
}

/**
 * Returns the total number of diagonals in the matrix
 * @param mtx matrix for which to return this value for
 * @return the number of non-zero diagonals
 */
JMTX_FAST_INT_T JMTX_NAME_TYPED(matrix_cds_diagonal_count)(const JMTX_NAME_TYPED(matrix_cds) * mtx)
{
    return (mtx->main_diagonal != NULL) + mtx->super_diagonals.count + mtx->sub_diagonals.count;
}

/**
 * Multiplies a dense column vector x by the sparse matrix and stores the result at y
 * @param mtx pointer to the memory where the matrix is stored
 * @param x pointer to vector to be multiplied
 * @param y pointer to vector where the result of multiplication is to be stored
 */
void JMTX_NAME_TYPED(matrix_cds_vector_multiply)(const JMTX_NAME_TYPED(matrix_cds) * mtx,
                                                 const JMTX_SCALAR_T *restrict x, JMTX_SCALAR_T *restrict y)
{
    if (mtx->base.cols > mtx->base.rows)
    {
        memset(y + mtx->base.rows, 0, sizeof(*y) * (mtx->base.cols - mtx->base.rows));
    }
    for (JMTX_FAST_INT_T i = 0; i < cds_superdiagonal_length(mtx->base.cols, mtx->base.rows, 0); ++i)
    {
        y[i] = mtx->main_diagonal[i] * x[i];
    }
    for (JMTX_FAST_INT_T i = 0; i < mtx->super_diagonals.count; ++i)
    {
        const JMTX_FAST_INT_T offset = mtx->super_diagonals.indices[i];
        const JMTX_SCALAR_T *const vals = mtx->super_diagonals.diagonals[i];
        for (JMTX_FAST_INT_T j = 0; j < cds_superdiagonal_length(mtx->base.cols, mtx->base.rows, offset); ++j)
        {
            y[j] += vals[j] * x[offset + j];
        }
    }
    for (JMTX_FAST_INT_T i = 0; i < mtx->sub_diagonals.count; ++i)
    {
        const JMTX_FAST_INT_T offset = mtx->sub_diagonals.indices[i];
        const JMTX_SCALAR_T *const vals = mtx->sub_diagonals.diagonals[i];
        for (JMTX_FAST_INT_T j = 0; j < cds_subdiagonal_length(mtx->base.cols, mtx->base.rows, offset); ++j)
        {
            y[offset + j] += vals[j] * x[j];
        }
    }
}

/**
 * Sets a single entry in the matrix. This is about as fast as setting the entire row of the matrix at once, if the
 * element was previously zero
 * @param mtx pointer to the memory where the matrix is stored
 * @param i row index
 * @param j column index
 * @param value value to which the value is set
 */
void JMTX_NAME_TYPED(matrix_cds_set_entry)(const JMTX_NAME_TYPED(matrix_cds) * mtx, JMTX_INDEX_T i, JMTX_INDEX_T j,
                                           JMTX_SCALAR_T value)
{
    if (i == j)
    {
        mtx->main_diagonal[i] = value;
    }
    else if (i > j)
    {
        //  Row idx larger than col idx -> subdiagonal
        JMTX_SCALAR_T *const ptr = JMTX_NAME_TYPED(matrix_cds_diagonal_array_get_ptr)(&mtx->sub_diagonals, i - j);
        ptr[j] = value;
    }
    else // if (i < j)
    {
        //  Row idx smaller than col idx -> superdiagonal
        JMTX_SCALAR_T *const ptr = JMTX_NAME_TYPED(matrix_cds_diagonal_array_get_ptr)(&mtx->super_diagonals, j - i);
        ptr[i] = value;
    }
}

/**
 * Returns a single entry from the matrix.
 * @param mtx pointer to the memory where the matrix is stored
 * @param i row index
 * @param j column index
 * @return value of the entry (0 if the entry was not manually set to anything else)
 */
JMTX_SCALAR_T JMTX_NAME_TYPED(matrix_cds_get_entry)(const JMTX_NAME_TYPED(matrix_cds) * mtx, JMTX_INDEX_T i,
                                                    JMTX_INDEX_T j)
{
    if (i == j)
    {
        if (mtx->main_diagonal)
        {
            return mtx->main_diagonal[i];
        }
    }
    else if (i > j)
    {
        //  Row idx larger than col idx -> subdiagonal
        JMTX_SCALAR_T *const ptr = JMTX_NAME_TYPED(matrix_cds_diagonal_array_get_ptr)(&mtx->sub_diagonals, i - j);
        if (ptr)
        {
            return ptr[j];
        }
    }
    else // if (i < j)
    {
        //  Row idx smaller than col idx -> superdiagonal
        JMTX_SCALAR_T *const ptr = JMTX_NAME_TYPED(matrix_cds_diagonal_array_get_ptr)(&mtx->super_diagonals, j - i);
        if (ptr)
        {
            return ptr[i];
        }
    }
    return 0.0f;
}

/**
 * Inserts and entry into the matrix, even if diagonal was not present before in the matrix
 * @param mtx pointer to the memory where the matrix is stored
 * @param i row index
 * @param j column index
 * @param value value to which the value is to be added
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result JMTX_NAME_TYPED(matrix_cds_insert_entry)(JMTX_NAME_TYPED(matrix_cds) * mtx, JMTX_INDEX_T i, JMTX_INDEX_T j,
                                                     JMTX_SCALAR_T value)
{
    int32_t dia = (int32_t)j - (int32_t)i;
    JMTX_SCALAR_T *const ptr = JMTX_NAME_TYPED(matrix_cds_allocate_diagonal)(mtx, dia, NULL);
    if (!ptr)
    {
        return JMTX_RESULT_BAD_ALLOC;
    }

    if (dia < 0)
    {
        //  Row idx larger than col idx -> subdiagonal
        ptr[j] = value;
    }
    else // if (dia >= 0)
    {
        //  Row idx smaller than or equal to col idx -> superdiagonal ro main diagonal
        ptr[i] = value;
    }
    return JMTX_RESULT_SUCCESS;
}

/**
 * Adds a value to an entry in the matrix when it exists or sets it to that value if it does not.
 * @param mtx pointer to the memory where the matrix is stored
 * @param i row index
 * @param j column index
 * @param value value to which the value is to be added
 */
jmtx_result JMTX_NAME_TYPED(matrix_cds_add_to_entry)(JMTX_NAME_TYPED(matrix_cds) * mtx, JMTX_INDEX_T i, JMTX_INDEX_T j,
                                                     JMTX_SCALAR_T value)
{
    int32_t dia = (int32_t)j - (int32_t)i;
    JMTX_SCALAR_T *const ptr = JMTX_NAME_TYPED(matrix_cds_allocate_diagonal)(mtx, dia, NULL);
    if (!ptr)
    {
        return JMTX_RESULT_BAD_ALLOC;
    }

    if (dia < 0)
    {
        //  Row idx larger than col idx -> subdiagonal
        ptr[j] += value;
    }
    else // if (dia >= 0)
    {
        //  Row idx smaller than or equal to col idx -> superdiagonal ro main diagonal
        ptr[i] += value;
    }
    return JMTX_RESULT_SUCCESS;
}

/**
 * Counts the number of times a specific value occurs in the matrix
 * @param mtx pointer to the memory where the matrix is stored
 * @param v value which to search for
 * @return number of times the value appeared in the matrix
 */
JMTX_INDEX_T JMTX_NAME_TYPED(matrix_cds_count_values)(const JMTX_NAME_TYPED(matrix_cds) * mtx, JMTX_SCALAR_T v)
{
    JMTX_FAST_INT_T k = 0;
    for (JMTX_FAST_INT_T i = 0; i < mtx->super_diagonals.count; ++i)
    {
        const JMTX_SCALAR_T *const ptr = mtx->super_diagonals.diagonals[i];
        JMTX_FAST_INT_T const idx = mtx->super_diagonals.indices[i];
        for (JMTX_FAST_INT_T j = 0; j < cds_superdiagonal_length(mtx->base.cols, mtx->base.rows, idx); ++j)
        {
            if (v == ptr[j])
            {
                k += 1;
            }
        }
    }
    for (JMTX_FAST_INT_T i = 0; i < mtx->sub_diagonals.count; ++i)
    {
        const JMTX_SCALAR_T *const ptr = mtx->sub_diagonals.diagonals[i];
        JMTX_FAST_INT_T const idx = mtx->sub_diagonals.indices[i];
        for (JMTX_FAST_INT_T j = 0; j < cds_subdiagonal_length(mtx->base.cols, mtx->base.rows, idx); ++j)
        {
            if (v == ptr[j])
            {
                k += 1;
            }
        }
    }
    if (mtx->main_diagonal)
    {
        for (JMTX_FAST_INT_T j = 0; j < cds_subdiagonal_length(mtx->base.cols, mtx->base.rows, 0); ++j)
        {
            if (v == mtx->main_diagonal[j])
            {
                k += 1;
            }
        }
    }
    return k;
}

/**
 * Zeros all entries within a matrix, but does not remove them in case they need to be reused
 * @param mtx matrix to zero
 */
void JMTX_NAME_TYPED(matrix_cds_zero_all_entries)(const JMTX_NAME_TYPED(matrix_cds) * mtx)
{
    for (JMTX_FAST_INT_T i = 0; i < mtx->super_diagonals.count; ++i)
    {
        JMTX_SCALAR_T *const ptr = mtx->super_diagonals.diagonals[i];
        JMTX_FAST_INT_T const idx = mtx->super_diagonals.indices[i];
        memset(ptr, 0, sizeof(*ptr) * cds_superdiagonal_length(mtx->base.cols, mtx->base.rows, idx));
    }
    for (JMTX_FAST_INT_T i = 0; i < mtx->sub_diagonals.count; ++i)
    {
        JMTX_SCALAR_T *const ptr = mtx->sub_diagonals.diagonals[i];
        JMTX_FAST_INT_T const idx = mtx->sub_diagonals.indices[i];
        memset(ptr, 0, sizeof(*ptr) * cds_subdiagonal_length(mtx->base.cols, mtx->base.rows, idx));
    }
    if (mtx->main_diagonal)
    {
        memset(mtx->main_diagonal, 0,
               sizeof(*mtx->main_diagonal) * cds_subdiagonal_length(mtx->base.cols, mtx->base.rows, 0));
    }
}

/**
 * Similar to JMTX_NAME_TYPED(matrix_cds_zero_all_entries, but slower, since it can not use memset. On the other hand,
 * it allows for the value to be other than 0
 * @param mtx matrix to set
 * @param x value to which to set all entries to
 */
void JMTX_NAME_TYPED(matrix_cds_set_all_entries)(const JMTX_NAME_TYPED(matrix_cds) * mtx, JMTX_SCALAR_T x)
{
    for (JMTX_FAST_INT_T i = 0; i < mtx->super_diagonals.count; ++i)
    {
        JMTX_SCALAR_T *const ptr = mtx->super_diagonals.diagonals[i];
        JMTX_FAST_INT_T const idx = mtx->super_diagonals.indices[i];
        for (JMTX_FAST_INT_T j = 0; j < cds_superdiagonal_length(mtx->base.cols, mtx->base.rows, idx); ++j)
        {
            ptr[j] = x;
        }
    }
    for (JMTX_FAST_INT_T i = 0; i < mtx->sub_diagonals.count; ++i)
    {
        JMTX_SCALAR_T *const ptr = mtx->sub_diagonals.diagonals[i];
        JMTX_FAST_INT_T const idx = mtx->sub_diagonals.indices[i];
        for (JMTX_FAST_INT_T j = 0; j < cds_subdiagonal_length(mtx->base.cols, mtx->base.rows, idx); ++j)
        {
            ptr[j] = x;
        }
    }
    if (mtx->main_diagonal)
    {
        for (JMTX_FAST_INT_T j = 0; j < cds_subdiagonal_length(mtx->base.cols, mtx->base.rows, 0); ++j)
        {
            mtx->main_diagonal[j] = x;
        }
    }
}

/**
 * Keeps the memory for the matrix, but sets the entry count to 0, so that matrix can be rebuilt.
 * @param mtx matrix to clear
 */
void JMTX_NAME_TYPED(matrix_cds_clear)(JMTX_NAME_TYPED(matrix_cds) * mtx)
{
    const jmtx_allocator_callbacks allocator_callbacks = mtx->base.allocator_callbacks;
    JMTX_FAST_INT_T count = mtx->super_diagonals.count;
    for (JMTX_FAST_INT_T i = 0; i < count; ++i)
    {
        allocator_callbacks.free(allocator_callbacks.state, mtx->super_diagonals.diagonals[i]);
    }
    mtx->super_diagonals.count = 0;
    count = mtx->sub_diagonals.count;
    for (JMTX_FAST_INT_T i = 0; i < count; ++i)
    {
        allocator_callbacks.free(allocator_callbacks.state, mtx->sub_diagonals.diagonals[i]);
    }
    mtx->sub_diagonals.count = 0;

    if (mtx->main_diagonal)
    {
        allocator_callbacks.free(allocator_callbacks.state, mtx->main_diagonal);
        mtx->main_diagonal = NULL;
    }
}

/**
 * Creates a transpose of a matrix
 * @param mtx pointer to the memory where the input matrix is stored
 * @param p_out address where the pointer to the output matrix will be returned
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result JMTX_NAME_TYPED(matrix_cds_transpose)(const JMTX_NAME_TYPED(matrix_cds) * mtx,
                                                  JMTX_NAME_TYPED(matrix_cds) * *p_out,
                                                  const jmtx_allocator_callbacks *allocator_callbacks)
{
    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }
    JMTX_NAME_TYPED(matrix_cds) *tps = NULL;
    int32_t dummy = 0;
    const JMTX_INDEX_T new_cols = mtx->base.rows;
    const JMTX_INDEX_T new_rows = mtx->base.cols;
    const JMTX_FAST_INT_T dia_count = JMTX_NAME_TYPED(matrix_cds_diagonal_count)(mtx);
    int32_t *p_dia = &dummy;
    if (dia_count != 0)
    {
        p_dia = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*p_dia) * dia_count);
        if (!p_dia)
        {
            return JMTX_RESULT_BAD_ALLOC;
        }
        JMTX_FAST_INT_T k = 0;
        for (JMTX_FAST_INT_T i = 0; i < mtx->super_diagonals.count; ++i)
        {
            const JMTX_INDEX_T idx = mtx->super_diagonals.indices[i];
            const int32_t new_idx = -(int32_t)idx;
            p_dia[k++] = new_idx;
        }
        if (mtx->main_diagonal)
        {
            p_dia[k++] = 0;
        }
        for (JMTX_FAST_INT_T i = 0; i < mtx->sub_diagonals.count; ++i)
        {
            p_dia[k++] = (int32_t)mtx->sub_diagonals.indices[i];
        }
        assert(k == dia_count);
    }

    jmtx_result res = JMTX_NAME_TYPED(matrix_cds_new)(&tps, new_rows, new_cols, dia_count, p_dia, allocator_callbacks);
    allocator_callbacks->free(allocator_callbacks->state, p_dia);
    if (res != JMTX_RESULT_SUCCESS)
    {
        return res;
    }
    JMTX_NAME_TYPED(matrix_cds_zero_all_entries)(tps);
    for (JMTX_FAST_INT_T i = 0; i < mtx->super_diagonals.count; ++i)
    {
        const JMTX_FAST_INT_T idx = mtx->super_diagonals.indices[i];
        const int32_t new_idx = -(int32_t)idx;
        const JMTX_FAST_INT_T len = cds_superdiagonal_length(mtx->base.cols, mtx->base.rows, idx);
        JMTX_SCALAR_T const *dia_ptr = mtx->super_diagonals.diagonals[i];
        res = JMTX_NAME_TYPED(matrix_cds_set_diagonal_part)(tps, new_idx, 0, len, NULL, dia_ptr);
        if (res != JMTX_RESULT_SUCCESS)
        {
            JMTX_NAME_TYPED(matrix_cds_destroy)(tps);
            return res;
        }
    }
    if (mtx->main_diagonal)
    {
        const JMTX_FAST_INT_T len = cds_subdiagonal_length(mtx->base.cols, mtx->base.rows, 0);
        res = JMTX_NAME_TYPED(matrix_cds_set_diagonal_part)(tps, 0, 0, len, NULL, mtx->main_diagonal);
        if (res != JMTX_RESULT_SUCCESS)
        {
            JMTX_NAME_TYPED(matrix_cds_destroy)(tps);
            return res;
        }
    }
    for (JMTX_FAST_INT_T i = 0; i < mtx->sub_diagonals.count; ++i)
    {
        const JMTX_FAST_INT_T idx = mtx->sub_diagonals.indices[i];
        const int32_t new_idx = (int32_t)idx;
        const JMTX_FAST_INT_T len = cds_subdiagonal_length(mtx->base.cols, mtx->base.rows, idx);
        JMTX_SCALAR_T const *dia_ptr = mtx->sub_diagonals.diagonals[i];
        res = JMTX_NAME_TYPED(matrix_cds_set_diagonal_part)(tps, new_idx, 0, len, NULL, dia_ptr);
        if (res != JMTX_RESULT_SUCCESS)
        {
            JMTX_NAME_TYPED(matrix_cds_destroy)(tps);
            return res;
        }
    }
    *p_out = tps;

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
jmtx_result JMTX_NAME_TYPED(matrix_cds_copy)(const JMTX_NAME_TYPED(matrix_cds) * mtx,
                                             JMTX_NAME_TYPED(matrix_cds) * *p_out,
                                             const jmtx_allocator_callbacks *allocator_callbacks)
{
    JMTX_NAME_TYPED(matrix_cds) *cpy = NULL;
    int32_t dummy = 0;
    const JMTX_INDEX_T new_cols = mtx->base.rows;
    const JMTX_INDEX_T new_rows = mtx->base.cols;
    const JMTX_FAST_INT_T dia_count = JMTX_NAME_TYPED(matrix_cds_diagonal_count)(mtx);
    int32_t *p_dia = &dummy;
    if (dia_count != 0)
    {
        p_dia = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*p_dia) * dia_count);
        if (!p_dia)
        {
            return JMTX_RESULT_BAD_ALLOC;
        }
        JMTX_FAST_INT_T k = 0;
        for (JMTX_FAST_INT_T i = 0; i < mtx->sub_diagonals.count; ++i)
        {
            p_dia[k++] = -(int32_t)mtx->sub_diagonals.indices[i];
        }
        if (mtx->main_diagonal)
        {
            p_dia[k++] = 0;
        }
        for (JMTX_FAST_INT_T i = 0; i < mtx->super_diagonals.count; ++i)
        {
            p_dia[k++] = (int32_t)mtx->super_diagonals.indices[i];
        }
        assert(k == dia_count);
    }

    jmtx_result res = JMTX_NAME_TYPED(matrix_cds_new)(&cpy, new_rows, new_cols, dia_count, p_dia, allocator_callbacks);
    allocator_callbacks->free(allocator_callbacks->state, p_dia);
    if (res != JMTX_RESULT_SUCCESS)
    {
        return res;
    }

    for (JMTX_FAST_INT_T i = 0; i < mtx->super_diagonals.count; ++i)
    {
        const int32_t idx = (int32_t)mtx->super_diagonals.indices[i];
        JMTX_SCALAR_T const *dia_ptr = mtx->super_diagonals.diagonals[i];
        res = JMTX_NAME_TYPED(matrix_cds_set_diagonal_full)(cpy, idx, dia_ptr);
        if (res != JMTX_RESULT_SUCCESS)
        {
            JMTX_NAME_TYPED(matrix_cds_destroy)(cpy);
            return res;
        }
    }
    if (mtx->main_diagonal)
    {
        const JMTX_FAST_INT_T len = cds_subdiagonal_length(mtx->base.cols, mtx->base.rows, 0);
        res = JMTX_NAME_TYPED(matrix_cds_set_diagonal_part)(cpy, 0, 0, len, NULL, mtx->main_diagonal);
        if (res != JMTX_RESULT_SUCCESS)
        {
            JMTX_NAME_TYPED(matrix_cds_destroy)(cpy);
            return res;
        }
    }
    for (JMTX_FAST_INT_T i = 0; i < mtx->sub_diagonals.count; ++i)
    {
        const JMTX_FAST_INT_T idx = mtx->sub_diagonals.indices[i];
        const int32_t new_idx = -(int32_t)idx;
        JMTX_SCALAR_T const *dia_ptr = mtx->sub_diagonals.diagonals[i];
        res = JMTX_NAME_TYPED(matrix_cds_set_diagonal_full)(cpy, new_idx, dia_ptr);
        if (res != JMTX_RESULT_SUCCESS)
        {
            JMTX_NAME_TYPED(matrix_cds_destroy)(cpy);
            return res;
        }
    }
    *p_out = cpy;

    return JMTX_RESULT_SUCCESS;
}

/**
 * Sets all the values on the diagonal to same value
 * @param mtx matrix to remove the diagonal from from
 * @param dia_idx diagonal to remove
 * @param v value to set the diagonal to
 */
void JMTX_NAME_TYPED(matrix_cds_set_diagonal)(JMTX_NAME_TYPED(matrix_cds) * mtx, int32_t dia_idx, JMTX_SCALAR_T v)
{
    JMTX_SCALAR_T *const ptr = JMTX_NAME_TYPED(matrix_cds_get_diagonal)(mtx, dia_idx);
    if (ptr)
    {
        JMTX_FAST_INT_T len;
        if (dia_idx < 0)
        {
            len = cds_subdiagonal_length(mtx->base.cols, mtx->base.rows, -dia_idx);
        }
        else
        {
            len = cds_superdiagonal_length(mtx->base.cols, mtx->base.rows, dia_idx);
        }
        for (JMTX_FAST_INT_T i = 0; i < len; ++i)
        {
            ptr[i] = v;
        }
    }
}

/**
 * Sets all the values on the diagonal to zero
 * @param mtx matrix to remove the diagonal from from
 * @param dia_idx diagonal to remove
 */
void JMTX_NAME_TYPED(matrix_cds_zero_diagonal)(JMTX_NAME_TYPED(matrix_cds) * mtx, int32_t dia_idx)
{
    JMTX_SCALAR_T *const ptr = JMTX_NAME_TYPED(matrix_cds_get_diagonal)(mtx, dia_idx);
    if (ptr)
    {
        JMTX_FAST_INT_T len;
        if (dia_idx < 0)
        {
            len = cds_subdiagonal_length(mtx->base.cols, mtx->base.rows, -dia_idx);
        }
        else
        {
            len = cds_superdiagonal_length(mtx->base.cols, mtx->base.rows, dia_idx);
        }
        memset(ptr, 0, sizeof(*ptr) * len);
    }
}

void JMTX_NAME_TYPED(matrix_cds_zero_row)(const JMTX_NAME_TYPED(matrix_cds) * mtx, JMTX_INDEX_T row)
{
    for (JMTX_FAST_INT_T i = 0; i < mtx->sub_diagonals.count; ++i)
    {
        const JMTX_FAST_INT_T dia = mtx->sub_diagonals.indices[i];
        if (dia > row)
        {
            break;
        }
        JMTX_SCALAR_T *const ptr = mtx->sub_diagonals.diagonals[i];
        ptr[row - dia] = 0.0f;
    }
    if (mtx->main_diagonal && row < mtx->base.cols)
    {
        mtx->main_diagonal[row] = 0;
    }
    for (JMTX_FAST_INT_T i = 0; i < mtx->super_diagonals.count; ++i)
    {
        const JMTX_FAST_INT_T dia = mtx->super_diagonals.indices[i];
        if (mtx->base.cols - dia <= row)
        {
            break;
        }
        JMTX_SCALAR_T *const ptr = mtx->super_diagonals.diagonals[i];
        ptr[row] = 0.0f;
        assert(row < cds_superdiagonal_length(mtx->base.cols, mtx->base.rows, dia));
    }
}

void JMTX_NAME_TYPED(matrix_cds_zero_col)(const JMTX_NAME_TYPED(matrix_cds) * mtx, JMTX_INDEX_T col)
{
    for (JMTX_FAST_INT_T i = 0; i < mtx->sub_diagonals.count; ++i)
    {
        const JMTX_FAST_INT_T dia = mtx->sub_diagonals.indices[i];
        if (mtx->base.rows - dia <= col)
        {
            break;
        }
        JMTX_SCALAR_T *const ptr = mtx->sub_diagonals.diagonals[i];
        ptr[col] = 0.0f;
    }
    if (mtx->main_diagonal && col < mtx->base.rows)
    {
        mtx->main_diagonal[col] = 0;
    }
    for (JMTX_FAST_INT_T i = 0; i < mtx->super_diagonals.count; ++i)
    {
        const JMTX_FAST_INT_T dia = mtx->super_diagonals.indices[i];
        if (dia > col)
        {
            break;
        }
        JMTX_SCALAR_T *const ptr = mtx->super_diagonals.diagonals[i];
        ptr[col - dia] = 0.0f;
    }
}
