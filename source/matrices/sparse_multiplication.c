// Automatically generated from source/float/matrices/sparse_multiplication.c on Fri Dec  1 06:43:01 2023
//
// Created by jan on 2.11.2023.
//

#include "sparse_multiplication.h"
#include "band_row_major.h"
#include "sparse_column_compressed.h"
#include "sparse_diagonal_compressed.h"
#include "sparse_row_compressed.h"

/**
 * Multiplies CRS and CCS matrix together and saves the result into a CRS matrix
 * @param a CRS matrix
 * @param b CCS matrix
 * @param p_out pointer which receives the resulting CRS matrix
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on allocation failure
 */
jmtx_result JMTX_NAME_TYPED(multiply_matrix_crs)(const JMTX_NAME_TYPED(matrix_crs) * a,
                                                 const JMTX_NAME_TYPED(matrix_ccs) * b,
                                                 JMTX_NAME_TYPED(matrix_crs) * *p_out,
                                                 const jmtx_allocator_callbacks *allocator_callbacks)
{

    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }

    const JMTX_INDEX_T r_out = a->base.rows;
    const JMTX_INDEX_T c_out = b->base.cols;

    JMTX_NAME_TYPED(matrix_crs) * out;
    jmtx_result res =
        JMTX_NAME_TYPED(matrix_crs_new)(&out, r_out, c_out, a->n_entries + b->n_entries, allocator_callbacks);
    if (res != JMTX_RESULT_SUCCESS)
    {
        return res;
    }

    JMTX_INDEX_T capacity = 1 + a->n_entries / a->base.rows > 1 + b->n_entries / b->base.cols
                                ? 1 + a->n_entries / a->base.rows
                                : 1 + b->n_entries / b->base.cols;
    JMTX_INDEX_T count = 0;
    JMTX_INDEX_T *indices = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*indices) * capacity);
    if (!indices)
    {
        JMTX_NAME_TYPED(matrix_crs_destroy)(out);
        return JMTX_RESULT_BAD_ALLOC;
    }
    JMTX_SCALAR_T *values = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*values) * capacity);
    if (!values)
    {
        allocator_callbacks->free(allocator_callbacks->state, indices);
        JMTX_NAME_TYPED(matrix_crs_destroy)(out);
        return JMTX_RESULT_BAD_ALLOC;
    }

    //  Go row by row, so that the building of output matrix is done more efficiently
    for (JMTX_INDEX_T i = 0; i < r_out; ++i)
    {
        count = 0;
        JMTX_INDEX_T n_a, n_b;
        JMTX_INDEX_T *i_a, *i_b;
        JMTX_SCALAR_T *v_a, *v_b;
        n_a = JMTX_NAME_TYPED(matrix_crs_get_row)(a, i, &i_a, &v_a);
        for (JMTX_INDEX_T j = 0; j < c_out; ++j)
        {
            n_b = JMTX_NAME_TYPED(matrix_ccs_get_col)(b, j, &i_b, &v_b);
            JMTX_SCALAR_T v = 0;
            for (JMTX_INDEX_T k_a = 0, k_b = 0; k_a < n_a && k_b < n_b;)
            {
                if (i_a[k_a] == i_b[k_b])
                {
                    v += v_a[k_a] * v_b[k_b];
                    k_a += 1;
                    k_b += 1;
                }
                else if (i_a[k_a] < i_b[k_b])
                {
                    k_a += 1;
                }
                else
                {
                    // i_a[k_a] > i_b[k_b]
                    k_b += 1;
                }
            }
            if (v != 0)
            {
                if (count == capacity)
                {
                    const JMTX_INDEX_T new_capacity = capacity << 1;
                    JMTX_SCALAR_T *const new_v = allocator_callbacks->realloc(allocator_callbacks->state, values,
                                                                              sizeof(*values) * new_capacity);
                    if (!new_v)
                    {
                        allocator_callbacks->free(allocator_callbacks->state, values);
                        allocator_callbacks->free(allocator_callbacks->state, indices);
                        JMTX_NAME_TYPED(matrix_crs_destroy)(out);
                        return JMTX_RESULT_BAD_ALLOC;
                    }
                    values = new_v;
                    JMTX_INDEX_T *const new_i = allocator_callbacks->realloc(allocator_callbacks->state, indices,
                                                                             sizeof(*indices) * new_capacity);
                    if (!new_i)
                    {
                        allocator_callbacks->free(allocator_callbacks->state, values);
                        allocator_callbacks->free(allocator_callbacks->state, indices);
                        JMTX_NAME_TYPED(matrix_crs_destroy)(out);
                        return JMTX_RESULT_BAD_ALLOC;
                    }
                    indices = new_i;

                    capacity = new_capacity;
                }
                indices[count] = j;
                values[count] = v;
                count += 1;
            }
        }
        JMTX_NAME_TYPED(matrix_crs_build_row)(out, i, count, indices, values);
    }

    allocator_callbacks->free(allocator_callbacks->state, values);
    allocator_callbacks->free(allocator_callbacks->state, indices);
    *p_out = out;
    return JMTX_RESULT_SUCCESS;
}

/**
 * Multiplies CRS and CCS matrix together and saves the result into a CCS matrix
 * @param a CRS matrix
 * @param b CCS matrix
 * @param p_out pointer which receives the resulting CCS matrix
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on allocation failure
 */
jmtx_result JMTX_NAME_TYPED(multiply_matrix_ccs)(const JMTX_NAME_TYPED(matrix_crs) * a,
                                                 const JMTX_NAME_TYPED(matrix_ccs) * b,
                                                 JMTX_NAME_TYPED(matrix_ccs) * *p_out,
                                                 const jmtx_allocator_callbacks *allocator_callbacks)
{
    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }

    const JMTX_INDEX_T r_out = a->base.rows;
    const JMTX_INDEX_T c_out = b->base.cols;

    JMTX_NAME_TYPED(matrix_ccs) * out;
    jmtx_result res =
        JMTX_NAME_TYPED(matrix_ccs_new)(&out, r_out, c_out, a->n_entries + b->n_entries, allocator_callbacks);
    if (res != JMTX_RESULT_SUCCESS)
    {
        return res;
    }

    JMTX_INDEX_T capacity = 1 + a->n_entries / a->base.rows > 1 + b->n_entries / b->base.cols
                                ? a->n_entries / a->base.rows
                                : b->n_entries / b->base.cols;
    JMTX_INDEX_T count;
    JMTX_INDEX_T *indices = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*indices) * capacity);
    if (!indices)
    {
        JMTX_NAME_TYPED(matrix_ccs_destroy)(out);
        return JMTX_RESULT_BAD_ALLOC;
    }
    JMTX_SCALAR_T *values = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*values) * capacity);
    if (!values)
    {
        allocator_callbacks->free(allocator_callbacks->state, indices);
        JMTX_NAME_TYPED(matrix_ccs_destroy)(out);
        return JMTX_RESULT_BAD_ALLOC;
    }

    //  Go column by column, so that the building of output matrix is done more efficiently
    for (JMTX_INDEX_T j = 0; j < c_out; ++j)
    {
        count = 0;
        JMTX_INDEX_T n_a, n_b;
        JMTX_INDEX_T *i_a, *i_b;
        JMTX_SCALAR_T *v_a, *v_b;
        n_b = JMTX_NAME_TYPED(matrix_ccs_get_col)(b, j, &i_b, &v_b);
        for (JMTX_INDEX_T i = 0; i < r_out; ++i)
        {
            n_a = JMTX_NAME_TYPED(matrix_crs_get_row)(a, i, &i_a, &v_a);
            JMTX_SCALAR_T v = 0;
            for (JMTX_INDEX_T k_a = 0, k_b = 0; k_a < n_a && k_b < n_b;)
            {
                if (i_a[k_a] == i_b[k_b])
                {
                    v += v_a[k_a] * v_b[k_b];
                    k_a += 1;
                    k_b += 1;
                }
                else if (i_a[k_a] < i_b[k_b])
                {
                    k_b += 1;
                }
                else
                {
                    // i_a[k_a] > i_b[k_b]
                    k_a += 1;
                }
            }
            if (v != 0)
            {
                if (count == capacity)
                {
                    const JMTX_INDEX_T new_capacity = capacity;
                    JMTX_SCALAR_T *const new_v = allocator_callbacks->realloc(allocator_callbacks->state, values,
                                                                              sizeof(*values) * new_capacity);
                    if (!new_v)
                    {
                        allocator_callbacks->free(allocator_callbacks->state, values);
                        allocator_callbacks->free(allocator_callbacks->state, indices);
                        JMTX_NAME_TYPED(matrix_ccs_destroy)(out);
                        return JMTX_RESULT_BAD_ALLOC;
                    }
                    values = new_v;
                    JMTX_INDEX_T *const new_i = allocator_callbacks->realloc(allocator_callbacks->state, indices,
                                                                             sizeof(*indices) * new_capacity);
                    if (!new_i)
                    {
                        allocator_callbacks->free(allocator_callbacks->state, values);
                        allocator_callbacks->free(allocator_callbacks->state, indices);
                        JMTX_NAME_TYPED(matrix_ccs_destroy)(out);
                        return JMTX_RESULT_BAD_ALLOC;
                    }
                    indices = new_i;

                    capacity = new_capacity;
                }
                indices[count] = j;
                values[count] = v;
                count += 1;
            }
        }
        JMTX_NAME_TYPED(matrix_ccs_build_col)(out, j, count, indices, values);
    }

    allocator_callbacks->free(allocator_callbacks->state, values);
    allocator_callbacks->free(allocator_callbacks->state, indices);
    *p_out = out;
    return JMTX_RESULT_SUCCESS;
}

JMTX_SCALAR_T JMTX_NAME_TYPED(multiply_matrix_sparse_vectors)(
    JMTX_INDEX_T n_a, const JMTX_INDEX_T i_a[JMTX_ARRAY_ATTRIB(const static n_a)],
    const JMTX_SCALAR_T v_a[JMTX_ARRAY_ATTRIB(const static n_a)], JMTX_INDEX_T n_b,
    const JMTX_INDEX_T i_b[JMTX_ARRAY_ATTRIB(const static n_b)],
    const JMTX_SCALAR_T v_b[JMTX_ARRAY_ATTRIB(const static n_b)])
{
    JMTX_SCALAR_T v = 0;
    JMTX_INDEX_T ia = 0, ib = 0;
    while (ia < n_a && ib < n_b)
    {
        if (i_a[ia] == i_b[ib])
        {
            v += v_a[ia] * v_b[ib];
            ia += 1;
            ib += 1;
        }
        else if (i_a[ia] < i_b[ib])
        {
            ia += 1;
        }
        else // if (i_a[ia] > i_b[ib])
        {
            ib += 1;
        }
    }
    return v;
}

JMTX_SCALAR_T JMTX_NAME_TYPED(multiply_matrix_sparse_vectors_limit)(
    JMTX_INDEX_T max_a, JMTX_INDEX_T max_b, JMTX_INDEX_T n_a, const JMTX_INDEX_T i_a[JMTX_ARRAY_ATTRIB(static n_a)],
    const JMTX_SCALAR_T v_a[JMTX_ARRAY_ATTRIB(static max_a)], JMTX_INDEX_T n_b,
    const JMTX_INDEX_T i_b[JMTX_ARRAY_ATTRIB(static n_b)], const JMTX_SCALAR_T v_b[JMTX_ARRAY_ATTRIB(static max_b)])
{
    JMTX_SCALAR_T v = 0;
    JMTX_INDEX_T ia = 0, ib = 0;
    while (ia < n_a && ib < n_b && i_a[ia] < max_a && i_b[ib] < max_b)
    {
        if (i_a[ia] == i_b[ib])
        {
            v += v_a[ia] * v_b[ib];
            i_a += 1;
            i_b += 1;
        }
        else if (i_a[ia] < i_b[ib])
        {
            i_a += 1;
        }
        else // if (i_a[ia] > i_b[ib])
        {
            i_b += 1;
        }
    }
    return v;
}

/**
 * Multiplies two BRM matrices together and produces a BRM matrix with the result of the matrix multiplication.
 * @param a BRM matrix
 * @param b BRM matrix
 * @param p_out pointer which receives the resulting BRM matrix
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on allocation failure
 */
jmtx_result JMTX_NAME_TYPED(multiply_matrix_brm)(const JMTX_NAME_TYPED(matrix_brm) * a,
                                                 const JMTX_NAME_TYPED(matrix_brm) * b,
                                                 JMTX_NAME_TYPED(matrix_brm) * *p_out,
                                                 const jmtx_allocator_callbacks *allocator_callbacks)
{

    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }

    const JMTX_INDEX_T r_out = a->base.rows;
    const JMTX_INDEX_T c_out = b->base.cols;

    JMTX_FAST_INT_T max_ubw;
    //    JMTX_FAST_INT_T min_ubw;
    if (a->upper_bandwidth > b->upper_bandwidth)
    {
        max_ubw = a->upper_bandwidth;
        //        min_ubw = b->upper_bandwidth;
    }
    else
    {
        max_ubw = b->upper_bandwidth;
        //        min_ubw = a->upper_bandwidth;
    }

    JMTX_FAST_INT_T max_lbw;
    //    JMTX_FAST_INT_T min_lbw;
    if (a->lower_bandwidth > b->lower_bandwidth)
    {
        max_lbw = a->lower_bandwidth;
        //        min_lbw = b->lower_bandwidth;
    }
    else
    {
        max_lbw = b->lower_bandwidth;
        //        min_lbw = a->lower_bandwidth;
    }

    JMTX_NAME_TYPED(matrix_brm) * out;
    jmtx_result res = JMTX_NAME_TYPED(matrix_brm_new)(&out, r_out, c_out, max_ubw, max_lbw, NULL, allocator_callbacks);
    if (res != JMTX_RESULT_SUCCESS)
    {
        return res;
    }

    const JMTX_INDEX_T capacity = 1 + max_ubw + max_lbw;
    JMTX_SCALAR_T *col_vals = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*col_vals) * capacity);
    if (!col_vals)
    {
        JMTX_NAME_TYPED(matrix_brm_destroy)(out);
        return JMTX_RESULT_BAD_ALLOC;
    }
    JMTX_SCALAR_T *values = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*col_vals) * capacity);
    if (!values)
    {
        allocator_callbacks->free(allocator_callbacks->state, col_vals);
        JMTX_NAME_TYPED(matrix_brm_destroy)(out);
        return JMTX_RESULT_BAD_ALLOC;
    }

    //  Go row by row, so that the building of output matrix is done more efficiently
    for (JMTX_INDEX_T i = 0; i < r_out; ++i)
    {
        JMTX_FAST_INT_T k = 0;
        JMTX_INDEX_T n_a, n_b;
        //        JMTX_INDEX_T* i_a, *i_b;
        JMTX_SCALAR_T *v_a; //, *v_b = col_vals;
        n_a = JMTX_NAME_TYPED(matrix_brm_get_row)(a, i, &v_a);
        (void)n_a, (void)n_b;
        const JMTX_FAST_INT_T first_a = JMTX_NAME_TYPED(matrix_brm_first_pos_in_row)(a, i);
        const JMTX_FAST_INT_T last_a = JMTX_NAME_TYPED(matrix_brm_last_pos_in_row)(a, i);
        for (JMTX_FAST_INT_T j = JMTX_NAME_TYPED(matrix_brm_first_pos_in_row)(out, i);
             j < JMTX_NAME_TYPED(matrix_brm_last_pos_in_row)(out, i) + 1; ++j, ++k)
        {
            n_b = JMTX_NAME_TYPED(matrix_brm_get_col)(b, j, col_vals);
            const JMTX_FAST_INT_T first_b = JMTX_NAME_TYPED(matrix_brm_first_pos_in_col)(b, j);
            const JMTX_FAST_INT_T last_b = JMTX_NAME_TYPED(matrix_brm_last_pos_in_col)(b, j);

            JMTX_FAST_INT_T first;

            JMTX_FAST_INT_T off_a = 0;
            JMTX_FAST_INT_T off_b = 0;

            int_fast32_t len;

            if (first_a > first_b)
            {
                first = first_a;
                off_b = first_b - first_a;
            }
            else
            {
                first = first_b;
                off_a = first_a - first_b;
            }

            if (last_a > last_b)
            {
                len = 1 + (int_fast32_t)(last_b) - (int_fast32_t)first;
            }
            else
            {
                len = 1 + (int_fast32_t)last_a - (int_fast32_t)first;
            }

            JMTX_SCALAR_T v = 0;
            for (int_fast32_t p = 0; p < len; ++p)
            {
                v += v_a[off_a + p] * col_vals[off_b + p];
            }
            values[k] = v;
        }

        JMTX_NAME_TYPED(matrix_brm_set_row)(out, i, values);
    }

    allocator_callbacks->free(allocator_callbacks->state, values);
    allocator_callbacks->free(allocator_callbacks->state, col_vals);
    *p_out = out;
    return JMTX_RESULT_SUCCESS;
}

/**
 * Multiplies two CDS matrices together and produces a CDS matrix with the result of the matrix multiplication.
 * @param a CDS matrix
 * @param b CDS matrix
 * @param p_out pointer which receives the resulting CDS matrix
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on allocation failure
 */
jmtx_result JMTX_NAME_TYPED(multiply_matrix_cds)(const JMTX_NAME_TYPED(matrix_cds) * a,
                                                 const JMTX_NAME_TYPED(matrix_cds) * b,
                                                 JMTX_NAME_TYPED(matrix_cds) * *p_out,
                                                 const jmtx_allocator_callbacks *allocator_callbacks)
{
    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }

    const JMTX_INDEX_T r_out = a->base.rows;
    const JMTX_INDEX_T c_out = b->base.cols;

    const JMTX_FAST_INT_T cnt_a = JMTX_NAME_TYPED(matrix_cds_diagonal_count)(a);
    const JMTX_FAST_INT_T cnt_b = JMTX_NAME_TYPED(matrix_cds_diagonal_count)(b);

    JMTX_INDEX_T *const idx_a = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*idx_a) * cnt_a);
    if (!idx_a)
    {
        return JMTX_RESULT_BAD_ALLOC;
    }
    JMTX_INDEX_T *const idx_b = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*idx_b) * cnt_b);
    if (!idx_b)
    {
        allocator_callbacks->free(allocator_callbacks->state, idx_a);
        return JMTX_RESULT_BAD_ALLOC;
    }
    JMTX_SCALAR_T *const val_a = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*val_a) * cnt_a);
    if (!val_a)
    {
        allocator_callbacks->free(allocator_callbacks->state, idx_b);
        allocator_callbacks->free(allocator_callbacks->state, idx_a);
        return JMTX_RESULT_BAD_ALLOC;
    }
    JMTX_SCALAR_T *const val_b = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*val_b) * cnt_b);
    if (!val_b)
    {
        allocator_callbacks->free(allocator_callbacks->state, val_a);
        allocator_callbacks->free(allocator_callbacks->state, idx_b);
        allocator_callbacks->free(allocator_callbacks->state, idx_a);
        return JMTX_RESULT_BAD_ALLOC;
    }

    JMTX_NAME_TYPED(matrix_cds) * out;
    jmtx_result res = JMTX_NAME_TYPED(matrix_cds_new)(&out, r_out, c_out, 0, (int32_t[]){0}, allocator_callbacks);
    if (res != JMTX_RESULT_SUCCESS)
    {
        allocator_callbacks->free(allocator_callbacks->state, val_b);
        allocator_callbacks->free(allocator_callbacks->state, val_a);
        allocator_callbacks->free(allocator_callbacks->state, idx_b);
        allocator_callbacks->free(allocator_callbacks->state, idx_a);
        return res;
    }

    //  CDS matrices are continuous in diagonals, meaning that there's no good way to build them, so just loop over each
    //  element. This might take long, but there's no better option for CDS matrices
    for (JMTX_INDEX_T i = 0; i < r_out; ++i)
    {
        const JMTX_INDEX_T c_a = JMTX_NAME_TYPED(matrix_cds_get_row)(a, i, cnt_a, val_a, idx_a);
        for (JMTX_FAST_INT_T j = 0; j < c_out; ++j)
        {
            const JMTX_INDEX_T c_b = JMTX_NAME_TYPED(matrix_cds_get_col)(b, j, cnt_b, val_b, idx_b);
            const JMTX_SCALAR_T val =
                JMTX_NAME_TYPED(multiply_matrix_sparse_vectors)(c_a, idx_a, val_a, c_b, idx_b, val_b);
            if (val != 0)
            {
                res = JMTX_NAME_TYPED(matrix_cds_insert_entry)(out, i, j, val);
                if (res != JMTX_RESULT_SUCCESS)
                {
                    allocator_callbacks->free(allocator_callbacks->state, val_b);
                    allocator_callbacks->free(allocator_callbacks->state, val_a);
                    allocator_callbacks->free(allocator_callbacks->state, idx_b);
                    allocator_callbacks->free(allocator_callbacks->state, idx_a);
                    return res;
                }
            }
        }
    }

    allocator_callbacks->free(allocator_callbacks->state, val_b);
    allocator_callbacks->free(allocator_callbacks->state, val_a);
    allocator_callbacks->free(allocator_callbacks->state, idx_b);
    allocator_callbacks->free(allocator_callbacks->state, idx_a);
    *p_out = out;
    return JMTX_RESULT_SUCCESS;
}
