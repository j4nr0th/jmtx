//
// Created by jan on 27.11.2023.
//

#include <assert.h>
#include <math.h>
#include "sparse_diagonal_compressed_internal.h"


enum {MINIMUM_RESERVED_DIAGONALS = 8};

static inline jmtx_result jmtx_matrix_cds_diagonal_array_init(const jmtx_allocator_callbacks* allocator_callbacks,
                                                              jmtx_matrix_cds_diagonal_array* ptr, uint32_t capacity)
{
    uint32_t* const idx = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*idx) * capacity);
    if (!idx)
    {
        return JMTX_RESULT_BAD_ALLOC;
    }
    float** const diags = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*diags) * capacity);
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

static inline jmtx_result jmtx_matrix_cds_diagonal_array_insert(const jmtx_allocator_callbacks* allocator_callbacks,
                                                                jmtx_matrix_cds_diagonal_array* this, uint32_t idx,
                                                                float* dia)
{
    assert(this->capacity != 0);
    if (this->capacity == this->count)
    {
        const uint32_t new_capacity = this->capacity << 1;
        uint32_t* const new_idx = allocator_callbacks->realloc(allocator_callbacks->state, this->indices, sizeof(*new_idx) * new_capacity);
        if (!idx)
        {
            return JMTX_RESULT_BAD_ALLOC;
        }
        this->indices = new_idx;
        float** const diags = allocator_callbacks->realloc(allocator_callbacks->state, this->diagonals, sizeof(*diags) * new_capacity);
        if (!diags)
        {
            return JMTX_RESULT_BAD_ALLOC;
        }
        this->diagonals = diags;
        this->capacity = new_capacity;
    }

    uint32_t* const indices = this->indices;
    float** const diagonals = this->diagonals;

    uint32_t pos;
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

static inline float* jmtx_matrix_cds_diagonal_array_get_ptr(const jmtx_matrix_cds_diagonal_array* this, uint32_t idx)
{
    for (uint32_t i = 0; i < this->count; ++i)
    {
        if (this->indices[i] == idx)
        {
            return this->diagonals[i];
        }
    }
    return NULL;
}

static inline void jmtx_matrix_cds_diagonal_array_free(const jmtx_allocator_callbacks* allocator_callbacks,
                                                       jmtx_matrix_cds_diagonal_array* ptr)
{
    for (uint32_t i = 0; i < ptr->count; ++i)
    {
        allocator_callbacks->free(allocator_callbacks->state, ptr->diagonals[i]);
    }
    allocator_callbacks->free(allocator_callbacks->state, ptr->diagonals);
    allocator_callbacks->free(allocator_callbacks->state, ptr->indices);
#ifndef NDEBUG
    ptr->diagonals = (void*)0xCCCCCCCCCCCCCCCC;
    ptr->indices = (void*)0xCCCCCCCCCCCCCCCC;
    ptr->count = 0xCCCCCCCC;
    ptr->capacity = 0xCCCCCCCC;
#endif
}


static inline uint_fast32_t cds_subdiagonal_length(const uint_fast32_t cols, const uint_fast32_t rows, const uint_fast32_t dia)
{
    assert(rows > dia);
    if (rows - dia > cols)
    {
        return cols;
    }
    return rows - dia;
}
static inline uint_fast32_t cds_superdiagonal_length(const uint_fast32_t cols, const uint_fast32_t rows, const uint_fast32_t dia)
{
    assert(cols > dia);
    if (cols - dia > rows)
    {
        return rows;
    }
    return cols - dia;
}
//
//static inline uint_fast32_t cds_subdiagonal_first_row(const uint_fast32_t cols, const uint_fast32_t rows, const uint_fast32_t dia)
//{
//    return dia;
//}
//static inline uint_fast32_t cds_subdiagonal_final_row(const uint_fast32_t cols, const uint_fast32_t rows, const uint_fast32_t dia)
//{
//    return dia + cds_subdiagonal_length(cols, rows, dia);
//}
//
//static inline uint_fast32_t cds_subdiagonal_first_col(const uint_fast32_t cols, const uint_fast32_t rows, const uint_fast32_t dia)
//{
//    return 0;
//}
//static inline uint_fast32_t cds_subdiagonal_final_col(const uint_fast32_t cols, const uint_fast32_t rows, const uint_fast32_t dia)
//{
//    return cds_subdiagonal_length(cols, rows, dia);
//}
//
//static inline uint_fast32_t cds_superdiagonal_first_row(const uint_fast32_t cols, const uint_fast32_t rows, const uint_fast32_t dia)
//{
//    return 0;
//}
//static inline uint_fast32_t cds_superdiagonal_final_row(const uint_fast32_t cols, const uint_fast32_t rows, const uint_fast32_t dia)
//{
//    return cds_superdiagonal_length(cols, rows, dia);
//}
//
//static inline uint_fast32_t cds_superdiagonal_first_col(const uint_fast32_t cols, const uint_fast32_t rows, const uint_fast32_t dia)
//{
//    return dia;
//}
//static inline uint_fast32_t cds_superdiagonal_final_col(const uint_fast32_t cols, const uint_fast32_t rows, const uint_fast32_t dia)
//{
//    return dia + cds_superdiagonal_length(cols, rows, dia);
//}

/**
 * Initializes a new Compressed Diagonal Sparse matrix
 * @param p_mtx address that receives the pointer to the matrix
 * @param cols number of columns of the sparse matrix
 * @param rows number of rows of the sparse matrix
 * @param reserved_entries how many entries should the space be reserved for in the matrix initially
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtx_matrix_cds_new(
        jmtx_matrix_cds** p_mtx, uint32_t cols, uint32_t rows, uint32_t n_diagonals,
        const int32_t p_dia_idx[static n_diagonals], const jmtx_allocator_callbacks* allocator_callbacks)
{
    if (!allocator_callbacks)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }
    if (n_diagonals == 0)
    {
        //  Always prepare the main diagonal
    }

    jmtx_result mtx_res;

    jmtx_matrix_cds* const mtx = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*mtx));
    if (!mtx)
    {
        return JMTX_RESULT_BAD_ALLOC;
    }

    uint32_t n_sub = 0, n_sup = 0;
    for (uint32_t i = 0; i < n_diagonals; ++i)
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

    mtx_res = jmtx_matrix_cds_diagonal_array_init(allocator_callbacks, &mtx->sub_diagonals,
                                                  n_sub < MINIMUM_RESERVED_DIAGONALS ? MINIMUM_RESERVED_DIAGONALS : n_sub);
    if (mtx_res != JMTX_RESULT_SUCCESS)
    {
        goto failed_sub;
    }

    mtx_res = jmtx_matrix_cds_diagonal_array_init(allocator_callbacks, &mtx->super_diagonals,
                                                  n_sup < MINIMUM_RESERVED_DIAGONALS ? MINIMUM_RESERVED_DIAGONALS : n_sup);
    if (mtx_res != JMTX_RESULT_SUCCESS)
    {
        goto failed_sup;
    }

    for (uint32_t i = 0; i < n_diagonals; ++i)
    {
        if (p_dia_idx[i] < 0)
        {
            float* const ptr = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*ptr) *
                    cds_subdiagonal_length(cols, rows, -p_dia_idx[i]));
            if (ptr == NULL)
            {
                mtx_res = JMTX_RESULT_BAD_ALLOC;
                goto failed_adding;
            }
            jmtx_matrix_cds_diagonal_array_insert(allocator_callbacks, &mtx->sub_diagonals, -p_dia_idx[i], ptr);
        }
        else if (p_dia_idx[i] > 0)
        {
            float* const ptr = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*ptr) *
                                                                                      cds_superdiagonal_length(cols, rows, p_dia_idx[i]));
            if (ptr == NULL)
            {
                mtx_res = JMTX_RESULT_BAD_ALLOC;
                goto failed_adding;
            }
            jmtx_matrix_cds_diagonal_array_insert(allocator_callbacks, &mtx->super_diagonals, p_dia_idx[i], ptr);
        }
    }
    if (n_diagonals > n_sub + n_sup)
    {
        assert(n_sup + n_sub + 1 == n_diagonals);
        const uint_fast32_t len_main = cds_subdiagonal_length(cols, rows, 0);
        assert(len_main == cds_superdiagonal_length(cols, rows, 0));
        float* const main_ptr = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*main_ptr) * len_main);
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

    mtx->base.type = JMTX_TYPE_CDS;
    mtx->base.cols = cols;
    mtx->base.rows = rows;
    mtx->base.allocator_callbacks = *allocator_callbacks;

    *p_mtx = mtx;

    return JMTX_RESULT_SUCCESS;
failed_adding:
    jmtx_matrix_cds_diagonal_array_free(allocator_callbacks, &mtx->sub_diagonals);
failed_sup:
    jmtx_matrix_cds_diagonal_array_free(allocator_callbacks, &mtx->super_diagonals);
failed_sub:
    allocator_callbacks->free(allocator_callbacks->state, mtx);

    return mtx_res;
}

/**
 * Cleans up the cds matrix and frees all of its memory
 * @param mtx pointer to memory where the matrix is stored
 */
void jmtx_matrix_cds_destroy(jmtx_matrix_cds* mtx)
{
    const jmtx_allocator_callbacks allocator_callbacks = mtx->base.allocator_callbacks;
    allocator_callbacks.free(allocator_callbacks.state, mtx->main_diagonal);
    jmtx_matrix_cds_diagonal_array_free(&allocator_callbacks, &mtx->sub_diagonals);
    jmtx_matrix_cds_diagonal_array_free(&allocator_callbacks, &mtx->super_diagonals);
    allocator_callbacks.free(allocator_callbacks.state, mtx);
}

/**
 * Version of jmtx_matrix_cds_set_row which does not touch the count of entries after the current row. This is useful when
 * building a new matrix, as it avoids unnecessary setting and resetting of these entries. Must be called for each row
 * in order to ensure that the matrix is properly built. Makes no checks on the input parameters
 * @param mtx pointer to the memory where the matrix is stored
 * @param dia_idx index of the diagonal to set
 * @param values values of non-zero values
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtx_matrix_cds_set_diagonal_full(jmtx_matrix_cds* mtx, int32_t dia_idx, const float values[const])
{
    uint32_t len;
    float* const ptr = jmtx_matrix_cds_allocate_diagonal(mtx, dia_idx, &len);
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
jmtx_result jmtx_matrix_cds_set_diagonal_part(jmtx_matrix_cds* mtx, int32_t dia_idx, uint32_t offset, uint32_t n,
                                              uint32_t* p_count, const float values[static n])
{
    uint32_t len;
    float* const ptr = jmtx_matrix_cds_allocate_diagonal(mtx, dia_idx, &len);
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
    float* const p = ptr + offset;
    for (uint_fast32_t i = 0; i < len - offset; ++i)
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
float* jmtx_matrix_cds_allocate_diagonal(jmtx_matrix_cds* mtx, int32_t dia, uint32_t* p_size)
{
    float* ptr;
    uint_fast32_t len;
    uint_fast32_t idx;
    if (dia == 0)
    {
        len = cds_superdiagonal_length(mtx->base.cols, mtx->base.rows, 0);
        //  Main diagonal
        if (mtx->main_diagonal == NULL)
        {
            mtx->main_diagonal = mtx->base.allocator_callbacks.alloc(mtx->base.allocator_callbacks.state, sizeof(*mtx->main_diagonal) * len);
            if (!mtx->main_diagonal)
            {
                NULL;
            }
        }
        ptr = mtx->main_diagonal;
    }
    else if (dia < 0)
    {
        idx = -dia;
        len = cds_subdiagonal_length(mtx->base.cols, mtx->base.rows, idx);
        ptr = jmtx_matrix_cds_diagonal_array_get_ptr(&mtx->sub_diagonals, idx);
        if (!ptr)
        {
            ptr = mtx->base.allocator_callbacks.alloc(mtx->base.allocator_callbacks.state, sizeof(*ptr) * len);
            if (!ptr)
            {
                return NULL;
            }
            const jmtx_result res = jmtx_matrix_cds_diagonal_array_insert(
                    &mtx->base.allocator_callbacks, &mtx->sub_diagonals, idx, ptr);
            if (res != JMTX_RESULT_SUCCESS)
            {
                mtx->base.allocator_callbacks.free(mtx->base.allocator_callbacks.state, ptr);
                return NULL;
            }
        }
    }
    else //if (dia_idx > 0)
    {
        idx = dia;
        len = cds_superdiagonal_length(mtx->base.cols, mtx->base.rows, idx);
        ptr = jmtx_matrix_cds_diagonal_array_get_ptr(&mtx->super_diagonals, idx);
        if (!ptr)
        {
            ptr = mtx->base.allocator_callbacks.alloc(mtx->base.allocator_callbacks.state, sizeof(*ptr) * len);
            if (!ptr)
            {
                return NULL;
            }
            const jmtx_result res = jmtx_matrix_cds_diagonal_array_insert(
                    &mtx->base.allocator_callbacks, &mtx->super_diagonals, idx, ptr);
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
 * Returns the pointer to the diagonal if one already exits, otherwise it returns NULL
 * @param mtx matrix which to allocate the diagonal for
 * @param dia offset of the diagonal from the main diagonal
 * @return pointer to the diagonal, or NULL in case it does not exist
 */
float* jmtx_matrix_cds_get_diagonal(const jmtx_matrix_cds* mtx, int32_t dia)
{
    if (dia == 0)
    {
        return mtx->main_diagonal;
    }
    else if (dia < 0)
    {
        return jmtx_matrix_cds_diagonal_array_get_ptr(&mtx->sub_diagonals, -dia);
    }
    //  else if (dia > 0)
    return jmtx_matrix_cds_diagonal_array_get_ptr(&mtx->super_diagonals, dia);
}

/**
 * Returns the number of entries in the row of the matrix
 * @param mtx pointer to the memory where the matrix is stored
 * @param row row index of the matrix to look at
 * @return number of entries in the row
 */
uint32_t jmtx_matrix_cds_entries_in_row(const jmtx_matrix_cds* mtx, uint32_t row)
{
    uint_fast32_t k = 0;
    for (uint_fast32_t i = 0; i < mtx->sub_diagonals.count; ++i)
    {
        //  Check the subdiagonals
        const uint_fast32_t d = mtx->sub_diagonals.indices[mtx->sub_diagonals.count - 1 - i];
        if (row < d || row >= d + cds_subdiagonal_length(mtx->base.cols, mtx->base.rows, d))
        {
            continue;
        }
        k += 1;
    }
    if (mtx->main_diagonal && row > mtx->base.cols)
    {
        k += 1;
    }
    for (uint_fast32_t i = 0; i < mtx->super_diagonals.count; ++i)
    {
        //  Check the superdiagonals
        const uint_fast32_t idx = mtx->super_diagonals.indices[i];
        if (row < mtx->base.cols - idx)
        {
            continue;
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
 * @param p_rows a buffer of at least n values which receives the row indices of the row
 * @return number of entries that were extracted from the column (may be less than are really in the column if n was too
 * small)
 */
uint32_t jmtx_matrix_cds_get_row(const jmtx_matrix_cds* mtx, uint32_t row, uint32_t n, float p_values[restrict n],
                                 uint32_t p_cols[restrict n])
{
    uint_fast32_t k = 0;
    for (uint_fast32_t i = 0; i < mtx->sub_diagonals.count && k < n; ++i)
    {
        //  Check the subdiagonals
        const uint_fast32_t d = mtx->sub_diagonals.indices[mtx->sub_diagonals.count - 1 - i];
        if (row < d || row >= d + cds_subdiagonal_length(mtx->base.cols, mtx->base.rows, d))
        {
            continue;
        }
        const uint_fast32_t pos = row - d;
        p_cols[k] = pos;
        const float* const diagonal = mtx->sub_diagonals.diagonals[mtx->sub_diagonals.count - 1 - i];
        p_values[k] = diagonal[pos];
        k += 1;
    }
    if (mtx->main_diagonal && row > mtx->base.cols && k < n)
    {
        p_cols[k] = row;
        p_values[k] = mtx->main_diagonal[row];
        k += 1;
    }
    for (uint_fast32_t i = 0; i < mtx->super_diagonals.count && k < n; ++i)
    {
        //  Check the superdiagonals
        const uint_fast32_t idx = mtx->super_diagonals.indices[i];
        if (row < mtx->base.cols - idx)
        {
            continue;
        }

        p_cols[k] = idx;
        const float* const diagonal = mtx->super_diagonals.diagonals[i];
        p_values[k] = diagonal[row];
        k += 1;
    }

    return k;
}

/**
 * Returns the number of entries in the column of the matrix
 * @param mtx pointer to the memory where the matrix is stored
 * @param col column index of the matrix to look at
 * @return number of entries in the column
 */
uint32_t jmtx_matrix_cds_entries_in_col(const jmtx_matrix_cds* mtx, uint32_t col){
    uint_fast32_t k = 0;
    for (uint_fast32_t i = 0; i < mtx->sub_diagonals.count; ++i)
    {
        //  Check the subdiagonals
        const uint_fast32_t idx = mtx->sub_diagonals.indices[i];
        if (col < mtx->base.rows - idx)
        {
            continue;
        }
        k += 1;
    }
    if (col > mtx->base.rows && mtx->main_diagonal)
    {
        k += 1;
    }
    for (uint_fast32_t i = 0; i < mtx->super_diagonals.count; ++i)
    {
        //  Check the superdiagonals
        const uint_fast32_t d = mtx->super_diagonals.indices[mtx->super_diagonals.count - 1 - i];
        if (col < d || col >= d + cds_superdiagonal_length(mtx->base.cols, mtx->base.rows, d))
        {
            continue;
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
uint32_t
jmtx_matrix_cds_get_col(const jmtx_matrix_cds* mtx, uint32_t col, uint32_t n, float p_values[restrict n], uint32_t p_rows[restrict n])
{
    uint_fast32_t k = 0;
    for (uint_fast32_t i = 0; i < mtx->sub_diagonals.count && k < n; ++i)
    {
        //  Check the subdiagonals
        const uint_fast32_t idx = mtx->sub_diagonals.indices[i];
        if (col < mtx->base.rows - idx)
        {
            continue;
        }

        p_rows[k] = idx;
        const float* const diagonal = mtx->sub_diagonals.diagonals[i];
        p_values[k] = diagonal[col];
        k += 1;
    }
    if (mtx->main_diagonal && col > mtx->base.rows && k < n)
    {
        p_rows[k] = col;
        p_values[k] = mtx->main_diagonal[col];
        k += 1;
    }
    for (uint_fast32_t i = 0; i < mtx->super_diagonals.count && k < n; ++i)
    {
        //  Check the superdiagonals
        const uint_fast32_t d = mtx->super_diagonals.indices[mtx->super_diagonals.count - 1 - i];
        if (col < d || col >= d + cds_superdiagonal_length(mtx->base.cols, mtx->base.rows, d))
        {
            continue;
        }
        const uint_fast32_t pos = col - d;
        p_rows[k] = pos;
        const float* const diagonal = mtx->super_diagonals.diagonals[mtx->super_diagonals.count - 1 - i];
        p_values[k] = diagonal[pos];
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
uint_fast32_t jmtx_matrix_cds_entries_in_dia(jmtx_matrix_cds* mtx, int32_t dia)
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
uint_fast32_t jmtx_matrix_cds_diagonal_count(const jmtx_matrix_cds* mtx)
{
    return (mtx->main_diagonal != NULL) + mtx->super_diagonals.count + mtx->sub_diagonals.count;
}

/**
 * Multiplies a dense column vector x by the sparse matrix and stores the result at y
 * @param mtx pointer to the memory where the matrix is stored
 * @param x pointer to vector to be multiplied
 * @param y pointer to vector where the result of multiplication is to be stored
 */
void jmtx_matrix_cds_vector_multiply(const jmtx_matrix_cds* mtx, const float* restrict x, float* restrict y)
{
    if (mtx->base.cols > mtx->base.rows)
    {
        memset(y + mtx->base.rows, 0, sizeof(*y) * (mtx->base.cols - mtx->base.rows));
    }
    for (uint_fast32_t i = 0; i < cds_superdiagonal_length(mtx->base.cols, mtx->base.rows, 0); ++i)
    {
        y[i] = mtx->main_diagonal[i] * x[i];
    }
    for (uint_fast32_t i = 0; i < mtx->super_diagonals.count; ++i)
    {
        const uint_fast32_t offset = mtx->super_diagonals.indices[i];
        const float* const vals = mtx->super_diagonals.diagonals[i];
        for (uint_fast32_t j = 0; j < cds_superdiagonal_length(mtx->base.cols, mtx->base.rows, offset); ++j)
        {
            y[j] += vals[j] * x[offset + j];
        }
    }
    for (uint_fast32_t i = 0; i < mtx->sub_diagonals.count; ++i)
    {
        const uint_fast32_t offset = mtx->sub_diagonals.indices[i];
        const float* const vals = mtx->sub_diagonals.diagonals[i];
        for (uint_fast32_t j = 0; j < cds_subdiagonal_length(mtx->base.cols, mtx->base.rows, offset); ++j)
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
void jmtx_matrix_cds_set_entry(const jmtx_matrix_cds* mtx, uint32_t i, uint32_t j, float value)
{
    if (i == j)
    {
        mtx->main_diagonal[i] = value;
    }
    else if (i > j)
    {
        //  Row idx larger than col idx -> subdiagonal
        float* const ptr = jmtx_matrix_cds_diagonal_array_get_ptr(&mtx->sub_diagonals, i - j);
        ptr[j] = value;
    }
    else// if (i < j)
    {
        //  Row idx smaller than col idx -> superdiagonal
        float* const ptr = jmtx_matrix_cds_diagonal_array_get_ptr(&mtx->super_diagonals, j - i);
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
float jmtx_matrix_cds_get_entry(const jmtx_matrix_cds* mtx, uint32_t i, uint32_t j)
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
        float* const ptr = jmtx_matrix_cds_diagonal_array_get_ptr(&mtx->sub_diagonals, i - j);
        if (ptr)
        {
            return ptr[j];
        }
    }
    else// if (i < j)
    {
        //  Row idx smaller than col idx -> superdiagonal
        float* const ptr = jmtx_matrix_cds_diagonal_array_get_ptr(&mtx->super_diagonals, j - i);
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
jmtx_result jmtx_matrix_cds_insert_entry(jmtx_matrix_cds* mtx, uint32_t i, uint32_t j, float value)
{
    int32_t dia = (int32_t)j - (int32_t)i;
    float* const ptr = jmtx_matrix_cds_allocate_diagonal(mtx, dia, NULL);
    if (!ptr)
    {
        return JMTX_RESULT_BAD_ALLOC;
    }

    if (dia < 0)
    {
        //  Row idx larger than col idx -> subdiagonal
        ptr[j] = value;
    }
    else// if (dia >= 0)
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
jmtx_result jmtx_matrix_cds_add_to_entry(jmtx_matrix_cds* mtx, uint32_t i, uint32_t j, float value)
{
    int32_t dia = (int32_t)j - (int32_t)i;
    float* const ptr = jmtx_matrix_cds_allocate_diagonal(mtx, dia, NULL);
    if (!ptr)
    {
        return JMTX_RESULT_BAD_ALLOC;
    }

    if (dia < 0)
    {
        //  Row idx larger than col idx -> subdiagonal
        ptr[j] += value;
    }
    else// if (dia >= 0)
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
uint32_t jmtx_matrix_cds_count_values(const jmtx_matrix_cds* mtx, float v)
{
    uint_fast32_t k = 0;
    for (uint_fast32_t i = 0; i < mtx->super_diagonals.count; ++i)
    {
        const float* const ptr = mtx->super_diagonals.diagonals[i];
        uint_fast32_t const idx = mtx->super_diagonals.indices[i];
        for (uint_fast32_t j = 0; j < cds_superdiagonal_length(mtx->base.cols, mtx->base.rows, idx); ++j)
        {
            if (v == ptr[j])
            {
                k += 1;
            }
        }
    }
    for (uint_fast32_t i = 0; i < mtx->sub_diagonals.count; ++i)
    {
        const float* const ptr = mtx->sub_diagonals.diagonals[i];
        uint_fast32_t const idx = mtx->sub_diagonals.indices[i];
        for (uint_fast32_t j = 0; j < cds_subdiagonal_length(mtx->base.cols, mtx->base.rows, idx); ++j)
        {
            if (v == ptr[j])
            {
                k += 1;
            }
        }
    }
    if (mtx->main_diagonal)
    {
        for (uint_fast32_t j = 0; j < cds_subdiagonal_length(mtx->base.cols, mtx->base.rows, 0); ++j)
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
void jmtx_matrix_cds_zero_all_entries(const jmtx_matrix_cds* mtx)
{
    for (uint_fast32_t i = 0; i < mtx->super_diagonals.count; ++i)
    {
        float* const ptr = mtx->super_diagonals.diagonals[i];
        uint_fast32_t const idx = mtx->super_diagonals.indices[i];
        memset(ptr, 0, sizeof(*ptr) * cds_superdiagonal_length(mtx->base.cols, mtx->base.rows, idx));
    }
    for (uint_fast32_t i = 0; i < mtx->sub_diagonals.count; ++i)
    {
        float* const ptr = mtx->sub_diagonals.diagonals[i];
        uint_fast32_t const idx = mtx->sub_diagonals.indices[i];
        memset(ptr, 0, sizeof(*ptr) * cds_subdiagonal_length(mtx->base.cols, mtx->base.rows, idx));

    }
    if (mtx->main_diagonal)
    {
        memset(mtx->main_diagonal, 0, sizeof(*mtx->main_diagonal) * cds_subdiagonal_length(mtx->base.cols, mtx->base.rows, 0));
    }
}

/**
 * Similar to jmtx_matrix_cds_zero_all_entries, but slower, since it can not use memset. On the other hand, it allows for
 * the value to be other than 0
 * @param mtx matrix to set
 * @param x value to which to set all entries to
 */
void jmtx_matrix_cds_set_all_entries(const jmtx_matrix_cds* mtx, float x)
{
    for (uint_fast32_t i = 0; i < mtx->super_diagonals.count; ++i)
    {
        float* const ptr = mtx->super_diagonals.diagonals[i];
        uint_fast32_t const idx = mtx->super_diagonals.indices[i];
        for (uint_fast32_t j = 0; j < cds_superdiagonal_length(mtx->base.cols, mtx->base.rows, idx); ++j)
        {
            ptr[j] = x;
        }
    }
    for (uint_fast32_t i = 0; i < mtx->sub_diagonals.count; ++i)
    {
        float* const ptr = mtx->sub_diagonals.diagonals[i];
        uint_fast32_t const idx = mtx->sub_diagonals.indices[i];
        for (uint_fast32_t j = 0; j < cds_subdiagonal_length(mtx->base.cols, mtx->base.rows, idx); ++j)
        {
            ptr[j] = x;
        }
    }
    if (mtx->main_diagonal)
    {
        for (uint_fast32_t j = 0; j < cds_subdiagonal_length(mtx->base.cols, mtx->base.rows, 0); ++j)
        {
            mtx->main_diagonal[j] = x;

        }
    }
}

/**
 * Keeps the memory for the matrix, but sets the entry count to 0, so that matrix can be rebuilt.
 * @param mtx matrix to clear
 */
void jmtx_matrix_cds_clear(jmtx_matrix_cds* mtx)
{
    const jmtx_allocator_callbacks allocator_callbacks = mtx->base.allocator_callbacks;
    uint_fast32_t count = mtx->super_diagonals.count;
    for (uint_fast32_t i = 0; i < count; ++i)
    {
        allocator_callbacks.free(allocator_callbacks.state, mtx->super_diagonals.diagonals[i]);
    }
    mtx->super_diagonals.count = 0;
    count = mtx->sub_diagonals.count;
    for (uint_fast32_t i = 0; i < count; ++i)
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
jmtx_result jmtx_matrix_cds_transpose(
        const jmtx_matrix_cds* mtx, jmtx_matrix_cds** p_out, const jmtx_allocator_callbacks* allocator_callbacks)
{
    jmtx_matrix_cds* tps = NULL;
    int32_t dummy = 0;
    const uint32_t new_cols = mtx->base.rows;
    const uint32_t new_rows = mtx->base.cols;
    const uint_fast32_t dia_count = jmtx_matrix_cds_diagonal_count(mtx);
    int32_t* p_dia = &dummy;
    if (dia_count != 0)
    {
        p_dia = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*p_dia) * dia_count);
        if (!p_dia)
        {
            return JMTX_RESULT_BAD_ALLOC;
        }
        uint_fast32_t k = 0;
        for (uint_fast32_t i = 0; i < mtx->super_diagonals.count; ++i)
        {
            const uint32_t idx = mtx->super_diagonals.indices[i];
            const int32_t new_idx = -(int32_t) idx;
            p_dia[k++] = new_idx;
        }
        if (mtx->main_diagonal)
        {
            p_dia[k++] = 0;
        }
        for (uint_fast32_t i = 0; i < mtx->sub_diagonals.count; ++i)
        {
            p_dia[k++] = (int32_t)mtx->sub_diagonals.indices[i];
        }
        assert(k == dia_count);
    }


    jmtx_result res = jmtx_matrix_cds_new(&tps, new_cols, new_rows, dia_count, p_dia, allocator_callbacks);
    allocator_callbacks->free(allocator_callbacks->state, p_dia);
    if (res != JMTX_RESULT_SUCCESS)
    {
        return res;
    }
    jmtx_matrix_cds_zero_all_entries(tps);
    for (uint_fast32_t i = 0; i < mtx->super_diagonals.count; ++i)
    {
        const uint_fast32_t idx = mtx->super_diagonals.indices[i];
        const int32_t new_idx = -(int32_t) idx;
        const uint_fast32_t len = cds_superdiagonal_length(mtx->base.cols, mtx->base.rows, idx);
        float const* dia_ptr = mtx->super_diagonals.diagonals[i];
        res = jmtx_matrix_cds_set_diagonal_part(tps, new_idx, 0, len, NULL, dia_ptr);
        if (res != JMTX_RESULT_SUCCESS)
        {
            jmtx_matrix_cds_destroy(tps);
            return res;
        }
    }
    if (mtx->main_diagonal)
    {
        const uint_fast32_t len = cds_subdiagonal_length(mtx->base.cols, mtx->base.rows, 0);
        res = jmtx_matrix_cds_set_diagonal_part(tps, 0, 0, len, NULL, mtx->main_diagonal);
        if (res != JMTX_RESULT_SUCCESS)
        {
            jmtx_matrix_cds_destroy(tps);
            return res;
        }
    }
    for (uint_fast32_t i = 0; i < mtx->sub_diagonals.count; ++i)
    {
        const uint_fast32_t idx = mtx->sub_diagonals.indices[i];
        const int32_t new_idx = (int32_t)idx;
        const uint_fast32_t len = cds_subdiagonal_length(mtx->base.cols, mtx->base.rows, idx);
        float const* dia_ptr = mtx->sub_diagonals.diagonals[i];
        res = jmtx_matrix_cds_set_diagonal_part(tps, new_idx, 0, len, NULL, dia_ptr);
        if (res != JMTX_RESULT_SUCCESS)
        {
            jmtx_matrix_cds_destroy(tps);
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
jmtx_result jmtx_matrix_cds_copy(const jmtx_matrix_cds* mtx, jmtx_matrix_cds** p_out,
                                 const jmtx_allocator_callbacks* allocator_callbacks)
{
    jmtx_matrix_cds* cpy = NULL;
    int32_t dummy = 0;
    const uint32_t new_cols = mtx->base.rows;
    const uint32_t new_rows = mtx->base.cols;
    const uint_fast32_t dia_count = jmtx_matrix_cds_diagonal_count(mtx);
    int32_t* p_dia = &dummy;
    if (dia_count != 0)
    {
        p_dia = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*p_dia) * dia_count);
        if (!p_dia)
        {
            return JMTX_RESULT_BAD_ALLOC;
        }
        uint_fast32_t k = 0;
        for (uint_fast32_t i = 0; i < mtx->sub_diagonals.count; ++i)
        {
            p_dia[k++] = -(int32_t)mtx->sub_diagonals.indices[i];
        }
        if (mtx->main_diagonal)
        {
            p_dia[k++] = 0;
        }
        for (uint_fast32_t i = 0; i < mtx->super_diagonals.count; ++i)
        {
            p_dia[k++] = (int32_t)mtx->super_diagonals.indices[i];
        }
        assert(k == dia_count);
    }


    jmtx_result res = jmtx_matrix_cds_new(&cpy, new_cols, new_rows, dia_count, p_dia, allocator_callbacks);
    allocator_callbacks->free(allocator_callbacks->state, p_dia);
    if (res != JMTX_RESULT_SUCCESS)
    {
        return res;
    }

    for (uint_fast32_t i = 0; i < mtx->super_diagonals.count; ++i)
    {
        const int32_t idx = (int32_t)mtx->super_diagonals.indices[i];
        float const* dia_ptr = mtx->super_diagonals.diagonals[i];
        res = jmtx_matrix_cds_set_diagonal_full(cpy, idx, dia_ptr);
        if (res != JMTX_RESULT_SUCCESS)
        {
            jmtx_matrix_cds_destroy(cpy);
            return res;
        }
    }
    if (mtx->main_diagonal)
    {
        const uint_fast32_t len = cds_subdiagonal_length(mtx->base.cols, mtx->base.rows, 0);
        res = jmtx_matrix_cds_set_diagonal_part(cpy, 0, 0, len, NULL, mtx->main_diagonal);
        if (res != JMTX_RESULT_SUCCESS)
        {
            jmtx_matrix_cds_destroy(cpy);
            return res;
        }
    }
    for (uint_fast32_t i = 0; i < mtx->sub_diagonals.count; ++i)
    {
        const uint_fast32_t idx = mtx->sub_diagonals.indices[i];
        const int32_t new_idx = -(int32_t)idx;
        float const* dia_ptr = mtx->sub_diagonals.diagonals[i];
        res = jmtx_matrix_cds_set_diagonal_full(cpy, new_idx, dia_ptr);
        if (res != JMTX_RESULT_SUCCESS)
        {
            jmtx_matrix_cds_destroy(cpy);
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
void jmtx_matrix_cds_set_diagonal(jmtx_matrix_cds* mtx, int32_t dia_idx, float v)
{
    float* const ptr = jmtx_matrix_cds_get_diagonal(mtx, dia_idx);
    if (ptr)
    {
        uint_fast32_t len;
        if (dia_idx < 0)
        {
            len = cds_subdiagonal_length(mtx->base.cols, mtx->base.rows, -dia_idx);
        }
        else
        {
            len = cds_superdiagonal_length(mtx->base.cols, mtx->base.rows, dia_idx);
        }
        for (uint_fast32_t i = 0; i < len; ++i)
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
void jmtx_matrix_cds_zero_diagonal(jmtx_matrix_cds* mtx, int32_t dia_idx)
{
    float* const ptr = jmtx_matrix_cds_get_diagonal(mtx, dia_idx);
    if (ptr)
    {
        uint_fast32_t len;
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


/**
 * Initializes a new Compressed Diagonal Sparse matrix
 * @param p_mtx address that receives the pointer to the matrix
 * @param cols number of columns of the sparse matrix
 * @param rows number of rows of the sparse matrix
 * @param n_diagonals how already filled diagonals are given to the matrix to initialize with
 * @param p_dia_idx SORTED indices of locations of the diagonals to reserve, indices being relative to the main
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result jmtxs_matrix_cds_new(
        jmtx_matrix_cds** p_mtx, uint32_t cols, uint32_t rows, uint32_t n_diagonals,
        const int32_t p_dia_idx[static n_diagonals], const jmtx_allocator_callbacks* allocator_callbacks)
{
    if (!p_mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (n_diagonals != 0)
    {
        if (p_dia_idx == NULL)
        {
            return JMTX_RESULT_NULL_PARAM;
        }
        //  check values are sorted and not duplicate
        for (uint_fast32_t i = 0; i < n_diagonals - 1; ++i)
        {
            if (p_dia_idx[i] >= p_dia_idx[i + 1])
            {
                return JMTX_RESULT_BAD_PARAM;
            }
        }
        //  check the values are within the allowed range
        for (uint_fast32_t i = 0; i < n_diagonals; ++i)
        {
            const int32_t idx = p_dia_idx[i];
            if (idx >= (int32_t)cols || idx <= -(int32_t)rows)
            {
                return JMTX_RESULT_INDEX_OUT_OF_BOUNDS;
            }
        }
    }
    if (allocator_callbacks && (!allocator_callbacks->alloc || !allocator_callbacks->free || !allocator_callbacks->realloc))
    {
        return JMTX_RESULT_BAD_PARAM;
    }
    return jmtx_matrix_cds_new(p_mtx, cols, rows, n_diagonals, p_dia_idx, allocator_callbacks);
}

/**
 * Cleans up the cds matrix and frees all of its memory
 * @param mtx pointer to memory where the matrix is stored
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result jmtxs_matrix_cds_destroy(jmtx_matrix_cds* mtx)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTX_TYPE_CDS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    jmtx_matrix_cds_destroy(mtx);
    return JMTX_RESULT_SUCCESS;
}

/**
 * Returns the total number of diagonals in the matrix
 * @param mtx matrix for which to return this value for
 * @param p_count pointer which receives the number of non-zero diagonals
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result jmtxs_matrix_cds_diagonal_count(const jmtx_matrix_cds* mtx, uint32_t* p_count)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!p_count)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTX_TYPE_CDS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    *p_count = jmtx_matrix_cds_diagonal_count(mtx);
    return JMTX_RESULT_SUCCESS;
}

/**
 * Returns the number of entries in the diagonal of the matrix
 * @param mtx pointer to the memory where the matrix is stored
 * @param dia index of the diagonal of the matrix to look at
 * @param p_count pointer which receives the number of entries in the diagonal
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result jmtxs_matrix_cds_entries_in_dia(jmtx_matrix_cds* mtx, int32_t dia, uint32_t* p_count)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!p_count)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTX_TYPE_CDS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (dia >= (int32_t)mtx->base.cols || dia <= -(int32_t)mtx->base.rows)
    {
        return JMTX_RESULT_INDEX_OUT_OF_BOUNDS;
    }
    *p_count = jmtx_matrix_cds_entries_in_dia(mtx, dia);
    return JMTX_RESULT_SUCCESS;
}

/**
 * Sets the entire diagonal to the values provided in the array. Will allocate a new diagonal if it was previously not
 * @param mtx pointer to the memory where the matrix is stored
 * @param dia_idx index of the diagonal to set
 * @param values values of non-zero values
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result jmtxs_matrix_cds_set_diagonal_full(jmtx_matrix_cds* mtx, int32_t dia_idx, const float values[])
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!values)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTX_TYPE_CDS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (dia_idx >= (int32_t)mtx->base.cols || dia_idx <= -(int32_t)mtx->base.rows)
    {
        return JMTX_RESULT_INDEX_OUT_OF_BOUNDS;
    }
    return jmtx_matrix_cds_set_diagonal_full(mtx, dia_idx, values);
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
jmtx_result jmtxs_matrix_cds_set_diagonal_part(jmtx_matrix_cds* mtx, int32_t dia_idx, uint32_t offset, uint32_t n,
                                               uint32_t* p_count, const float values[static n])
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!values)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTX_TYPE_CDS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (dia_idx >= (int32_t)mtx->base.cols || dia_idx <= -(int32_t)mtx->base.rows)
    {
        return JMTX_RESULT_INDEX_OUT_OF_BOUNDS;
    }
    return jmtx_matrix_cds_set_diagonal_part(mtx, dia_idx, offset, n, p_count, values);
}

/**
 * Allocates a new diagonal if one does not already exits, otherwise it returns the pointer to the existing one.
 * @param mtx matrix which to allocate the diagonal for
 * @param dia offset of the diagonal from the main diagonal
 * @param p_size pointer which receives the number of elements in the allocated diagonal. May be NULL
 * @param p_dia pointer which receives pointer to the newly allocated diagonal
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxs_matrix_cds_allocate_diagonal(jmtx_matrix_cds* mtx, int32_t dia, uint32_t* p_size, float** p_dia)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!p_dia)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTX_TYPE_CDS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (dia >= (int32_t)mtx->base.cols || dia <= -(int32_t)mtx->base.rows)
    {
        return JMTX_RESULT_INDEX_OUT_OF_BOUNDS;
    }
    float* const ptr = jmtx_matrix_cds_allocate_diagonal(mtx, dia, p_size);
    if (!ptr)
    {
        return JMTX_RESULT_BAD_ALLOC;
    }
    *p_dia = ptr;
    return JMTX_RESULT_SUCCESS;
}

/**
 * Returns the pointer to the diagonal if one already exits, otherwise it returns NULL
 * @param mtx matrix which to allocate the diagonal for
 * @param dia offset of the diagonal from the main diagonal
 * @param p_dia pointer which receives the pointer to the diagonal
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result jmtxs_matrix_cds_get_diagonal(const jmtx_matrix_cds* mtx, int32_t dia, float** p_dia)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTX_TYPE_CDS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (!p_dia)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (dia >= (int32_t)mtx->base.cols || dia <= -(int32_t)mtx->base.rows)
    {
        return JMTX_RESULT_INDEX_OUT_OF_BOUNDS;
    }
    *p_dia = jmtx_matrix_cds_get_diagonal(mtx, dia);
    return JMTX_RESULT_SUCCESS;
}

/**
 * Returns the number of entries in the row of the matrix
 * @param mtx pointer to the memory where the matrix is stored
 * @param row row index of the matrix to look at
 * @param p_count pointer which receives number of entries in the row
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result jmtxs_matrix_cds_entries_in_row(const jmtx_matrix_cds* mtx, uint32_t row, uint32_t* p_count)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTX_TYPE_CDS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (row >= mtx->base.rows)
    {
        return JMTX_RESULT_INDEX_OUT_OF_BOUNDS;
    }
    if (!p_count)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    *p_count = jmtx_matrix_cds_entries_in_row(mtx, row);
    return JMTX_RESULT_SUCCESS;
}

/**
 * Returns the values of entries in the matrix, along with what column of the matrix they were located in
 * @param mtx pointer to the memory where the matrix is stored
 * @param row row index of the matrix to look at
 * @param n number of values in the row to be extracted at most
 * @param p_values a buffer of at least n values which receives the values of the row
 * @param p_rows a buffer of at least n values which receives the row indices of the row
 * @param p_count pointer which receives the  number of entries that were extracted from the column (may be less than
 * are really in the column if n was too small)
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result jmtxs_matrix_cds_get_row(const jmtx_matrix_cds* mtx, uint32_t row, uint32_t n, float p_values[restrict n],
                                     uint32_t p_cols[restrict n], uint32_t* p_count)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTX_TYPE_CDS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (row >= mtx->base.rows)
    {
        return JMTX_RESULT_INDEX_OUT_OF_BOUNDS;
    }
    if (!p_count)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!p_values)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!p_cols)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    *p_count = jmtx_matrix_cds_get_row(mtx, row, n, p_values, p_cols);
    return JMTX_RESULT_SUCCESS;
}

/**
 * Returns the number of entries in the column of the matrix
 * @param mtx pointer to the memory where the matrix is stored
 * @param col column index of the matrix to look at
 * @param p_count pointer which receives number of entries in the column
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result jmtxs_matrix_cds_entries_in_col(const jmtx_matrix_cds* mtx, uint32_t col, uint32_t* p_count)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTX_TYPE_CDS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (col >= mtx->base.cols)
    {
        return JMTX_RESULT_INDEX_OUT_OF_BOUNDS;
    }
    if (!p_count)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    *p_count = jmtx_matrix_cds_entries_in_col(mtx, col);
    return JMTX_RESULT_SUCCESS;
}

/**
 * Returns the values of entries in the matrix, along with what row of the matrix they were located in
 * @param mtx pointer to the memory where the matrix is stored
 * @param col column index of the matrix to look at
 * @param n number of values in the column to be extracted at most
 * @param p_values a buffer of at least n values which receives the values of the column
 * @param p_rows a buffer of at least n values which receives the row indices of the column
 * @param p_count pointer which receives the  number of entries that were extracted from the column (may be less than
 * are really in the column if n was too small)
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result
jmtxs_matrix_cds_get_col(const jmtx_matrix_cds* mtx, uint32_t col, uint32_t n, float p_values[restrict n],
                         uint32_t p_rows[restrict n], uint32_t* p_count)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTX_TYPE_CDS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (col >= mtx->base.cols)
    {
        return JMTX_RESULT_INDEX_OUT_OF_BOUNDS;
    }
    if (!p_count)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!p_values)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!p_rows)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    *p_count = jmtx_matrix_cds_get_col(mtx, col, n, p_values, p_rows);
    return JMTX_RESULT_SUCCESS;
}

/**
 * Multiplies a dense column vector x by the sparse matrix and stores the result at y
 * @param mtx pointer to the memory where the matrix is stored
 * @param x pointer to vector to be multiplied
 * @param y pointer to vector where the result of multiplication is to be stored
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result jmtxs_matrix_cds_vector_multiply(const jmtx_matrix_cds* mtx, const float* restrict x, float* restrict y)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTX_TYPE_CDS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (!x)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!y)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    //  Check if pointers alias each other
    if (y < x + mtx->base.cols || x < y + mtx->base.cols)
    {
        return JMTX_RESULT_BAD_PARAM;
    }
    jmtx_matrix_cds_vector_multiply(mtx, x, y);
    return JMTX_RESULT_SUCCESS;
}

/**
 * Sets a single entry in the matrix. The diagonal that it is to be inserted on MUST exist before
 * @param mtx pointer to the memory where the matrix is stored
 * @param i row index
 * @param j column index
 * @param value value to which the value is set
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result jmtxs_matrix_cds_set_entry(const jmtx_matrix_cds* mtx, uint32_t i, uint32_t j, float value)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTX_TYPE_CDS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (i >= mtx->base.rows || j >= mtx->base.cols)
    {
        return JMTX_RESULT_INDEX_OUT_OF_BOUNDS;
    }
    if (!isfinite(value))
    {
        return JMTX_RESULT_BAD_PARAM;
    }
    jmtx_matrix_cds_set_entry(mtx, i, j, value);
    return JMTX_RESULT_SUCCESS;
}

/**
 * Returns a single entry from the matrix.
 * @param mtx pointer to the memory where the matrix is stored
 * @param i row index
 * @param j column index
 * @param p_v pointer which receives the value of the entry (0 if the entry was not manually set to anything else)
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result jmtxs_matrix_cds_get_entry(const jmtx_matrix_cds* mtx, uint32_t i, uint32_t j, float* p_v)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTX_TYPE_CDS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (i >= mtx->base.rows || j >= mtx->base.cols)
    {
        return JMTX_RESULT_INDEX_OUT_OF_BOUNDS;
    }
    if (!p_v)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    *p_v = jmtx_matrix_cds_get_entry(mtx, i, j);
    return JMTX_RESULT_SUCCESS;
}

/**
 * Inserts and entry into the matrix, even if diagonal was not present before in the matrix
 * @param mtx pointer to the memory where the matrix is stored
 * @param i row index
 * @param j column index
 * @param value value to which the value is to be added
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result jmtxs_matrix_cds_insert_entry(jmtx_matrix_cds* mtx, uint32_t i, uint32_t j, float value)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTX_TYPE_CDS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (i >= mtx->base.rows || j >= mtx->base.cols)
    {
        return JMTX_RESULT_INDEX_OUT_OF_BOUNDS;
    }
    if (!isfinite(value))
    {
        return JMTX_RESULT_BAD_PARAM;
    }
    return jmtx_matrix_cds_insert_entry(mtx, i, j, value);
}

/**
 * Adds a value to an entry in the matrix when it exists or sets it to that value if it does not.
 * @param mtx pointer to the memory where the matrix is stored
 * @param i row index
 * @param j column index
 * @param value value to which the value is to be added
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result jmtxs_matrix_cds_add_to_entry(jmtx_matrix_cds* mtx, uint32_t i, uint32_t j, float value)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTX_TYPE_CDS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (i >= mtx->base.rows || j >= mtx->base.cols)
    {
        return JMTX_RESULT_INDEX_OUT_OF_BOUNDS;
    }
    if (!isfinite(value))
    {
        return JMTX_RESULT_BAD_PARAM;
    }
    jmtx_matrix_cds_add_to_entry(mtx, i, j, value);
    return JMTX_RESULT_SUCCESS;
}

/**
 * Counts the number of times a specific value occurs in the matrix
 * @param mtx pointer to the memory where the matrix is stored
 * @param v value which to search for
 * @param p_count pointer which receives the number of times the value appeared in the matrix
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result jmtxs_matrix_cds_count_values(const jmtx_matrix_cds* mtx, float v, uint32_t* p_count)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTX_TYPE_CDS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (!p_count)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    *p_count = jmtx_matrix_cds_count_values(mtx, v);
    return JMTX_RESULT_SUCCESS;
}

/**
 * Zeros all entries within a matrix, but does not remove them in case they need to be reused
 * @param mtx matrix to zero
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result jmtxs_matrix_cds_zero_all_entries(const jmtx_matrix_cds* mtx)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTX_TYPE_CDS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    jmtx_matrix_cds_zero_all_entries(mtx);
    return JMTX_RESULT_SUCCESS;
}

/**
 * Similar to jmtx_matrix_cds_zero_all_entries, but slower, since it can not use memset. On the other hand, it allows for
 * the value to be other than 0
 * @param mtx matrix to set
 * @param x value to which to set all entries to
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result jmtxs_matrix_cds_set_all_entries(const jmtx_matrix_cds* mtx, float x)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTX_TYPE_CDS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (!isfinite(x))
    {
        return JMTX_RESULT_BAD_PARAM;
    }
    jmtx_matrix_cds_set_all_entries(mtx, x);
    return JMTX_RESULT_SUCCESS;
}

/**
 * Keeps the memory for the matrix, but sets the entry count to 0, so that matrix can be rebuilt.
 * @param mtx matrix to clear
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result jmtxs_matrix_cds_clear(jmtx_matrix_cds* mtx)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTX_TYPE_CDS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    jmtx_matrix_cds_clear(mtx);
    return JMTX_RESULT_SUCCESS;
}

/**
 * Creates a transpose of a matrix
 * @param mtx pointer to the memory where the input matrix is stored
 * @param p_out address where the pointer to the output matrix will be returned
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result jmtxs_matrix_cds_transpose(
        const jmtx_matrix_cds* mtx, jmtx_matrix_cds** p_out, const jmtx_allocator_callbacks* allocator_callbacks)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!p_out)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (allocator_callbacks && (!allocator_callbacks->alloc || !allocator_callbacks->free || !allocator_callbacks->realloc))
    {
        return JMTX_RESULT_BAD_PARAM;
    }
    return jmtx_matrix_cds_transpose(mtx, p_out, allocator_callbacks);
}

/**
 * Creates a copy of the matrix
 * @param mtx pointer to the memory where the input matrix is stored
 * @param p_out address where the pointer to the output matrix will be returned
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxs_matrix_cds_copy(const jmtx_matrix_cds* mtx, jmtx_matrix_cds** p_out, const jmtx_allocator_callbacks* allocator_callbacks)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!p_out)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (allocator_callbacks && (!allocator_callbacks->alloc || !allocator_callbacks->free || !allocator_callbacks->realloc))
    {
        return JMTX_RESULT_BAD_PARAM;
    }
    return jmtx_matrix_cds_copy(mtx, p_out, allocator_callbacks);
}


/**
 * Sets all the values on the diagonal to same value
 * @param mtx matrix to remove the diagonal from from
 * @param dia_idx diagonal to remove
 * @param v value to set the diagonal to
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result jmtxs_matrix_cds_set_diagonal(jmtx_matrix_cds* mtx, int32_t dia_idx, float v)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTX_TYPE_CDS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (dia_idx >= (int32_t)mtx->base.cols || dia_idx <= -(int32_t)mtx->base.rows)
    {
        return JMTX_RESULT_INDEX_OUT_OF_BOUNDS;
    }
    if (!isfinite(v))
    {
        return JMTX_RESULT_BAD_PARAM;
    }
    jmtx_matrix_cds_set_diagonal(mtx, dia_idx, v);
    return JMTX_RESULT_SUCCESS;
}


/**
 * Sets all the values on the diagonal to zero
 * @param mtx matrix to remove the diagonal from from
 * @param dia_idx diagonal to remove
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result jmtxs_matrix_cds_zero_diagonal(jmtx_matrix_cds* mtx, int32_t dia_idx)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTX_TYPE_CDS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (dia_idx >= (int32_t)mtx->base.cols || dia_idx <= -(int32_t)mtx->base.rows)
    {
        return JMTX_RESULT_INDEX_OUT_OF_BOUNDS;
    }
    jmtx_matrix_cds_zero_diagonal(mtx, dia_idx);
    return JMTX_RESULT_SUCCESS;
}





