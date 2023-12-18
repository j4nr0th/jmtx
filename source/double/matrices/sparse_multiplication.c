// Automatically generated from source/float/matrices/sparse_multiplication.c on Fri Dec  1 06:43:01 2023
//
// Created by jan on 2.11.2023.
//

#include "sparse_column_compressed_internal.h"
#include "sparse_row_compressed_internal.h"
#include "band_row_major_internal.h"
#include "sparse_diagonal_compressed_internal.h"
#include "../../../include/jmtx/double/matrices/sparse_multiplication.h"
#include "../../../include/jmtx/double/matrices/sparse_multiplication_safe.h"

/**
 * Multiplies CRS and CCS matrix together and saves the result into a CRS matrix
 * @param a CRS matrix
 * @param b CCS matrix
 * @param p_out pointer which receives the resulting CRS matrix
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on allocation failure
 */
jmtx_result jmtxd_multiply_matrix_crs(
        const jmtxd_matrix_crs* a, const jmtxd_matrix_ccs* b, jmtxd_matrix_crs** p_out,
        const jmtx_allocator_callbacks* allocator_callbacks)
{

    if (!allocator_callbacks)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }

    const uint32_t r_out = a->base.rows;
    const uint32_t c_out = b->base.cols;

    jmtxd_matrix_crs* out;
    jmtx_result res = jmtxd_matrix_crs_new(&out, c_out, r_out, a->n_entries + b->n_entries, allocator_callbacks);
    if (res != JMTX_RESULT_SUCCESS)
    {
        return res;
    }

    uint32_t capacity = 1 + a->n_entries / a->base.rows > 1 + b->n_entries / b->base.cols ?
                              1 + a->n_entries / a->base.rows :
                              1 + b->n_entries / b->base.cols;
    uint32_t count = 0;
    uint32_t* indices = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*indices) * capacity);
    if (!indices)
    {
        jmtxd_matrix_crs_destroy(out);
        return JMTX_RESULT_BAD_ALLOC;
    }
    double* values = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*values) * capacity);
    if (!values)
    {
        allocator_callbacks->free(allocator_callbacks->state, indices);
        jmtxd_matrix_crs_destroy(out);
        return JMTX_RESULT_BAD_ALLOC;
    }


    //  Go row by row, so that the building of output matrix is done more efficiently
    for (uint32_t i = 0; i < r_out; ++i)
    {
        count = 0;
        uint32_t n_a, n_b;
        uint32_t* i_a, *i_b;
        double* v_a, *v_b;
        n_a = jmtxd_matrix_crs_get_row(a, i, &i_a, &v_a);
        for (uint32_t j = 0; j < c_out; ++j)
        {
            n_b = jmtxd_matrix_ccs_get_col(b, j, &i_b, &v_b);
            double v = 0;
            for (uint32_t k_a = 0, k_b = 0; k_a < n_a && k_b < n_b;)
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
                    //i_a[k_a] > i_b[k_b]
                    k_b += 1;
                }
            }
            if (v != 0)
            {
                if (count == capacity)
                {
                    const uint32_t new_capacity = capacity << 1;
                    double* const new_v = allocator_callbacks->realloc(allocator_callbacks->state, values, sizeof(*values) * new_capacity);
                    if (!new_v)
                    {
                        allocator_callbacks->free(allocator_callbacks->state, values);
                        allocator_callbacks->free(allocator_callbacks->state, indices);
                        jmtxd_matrix_crs_destroy(out);
                        return JMTX_RESULT_BAD_ALLOC;
                    }
                    values = new_v;
                    uint32_t* const new_i = allocator_callbacks->realloc(allocator_callbacks->state, indices, sizeof(*indices) * new_capacity);
                    if (!new_i)
                    {
                        allocator_callbacks->free(allocator_callbacks->state, values);
                        allocator_callbacks->free(allocator_callbacks->state, indices);
                        jmtxd_matrix_crs_destroy(out);
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
        jmtxd_matrix_crs_build_row(out, i, count, indices, values);
    }


    allocator_callbacks->free(allocator_callbacks->state, values);
    allocator_callbacks->free(allocator_callbacks->state, indices);
    *p_out = out;
    return JMTX_RESULT_SUCCESS;
}

/**
 * Multiplies CRS and CCS matrix together and saves the result into a CRS matrix
 * @param a CRS matrix
 * @param b CCS matrix
 * @param p_out pointer which receives the resulting CRS matrix
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result jmtxds_multiply_matrix_crs(const jmtxd_matrix_crs* a, const jmtxd_matrix_ccs* b, jmtxd_matrix_crs** p_out,
                                       const jmtx_allocator_callbacks* allocator_callbacks)
{
    if (!a)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (a->base.type != JMTXD_TYPE_CRS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (!b)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (b->base.type != JMTXD_TYPE_CCS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (a->base.cols != b->base.rows)
    {
        //  can't do multiplication
        return JMTX_RESULT_BAD_MATRIX;
    }
    if (!p_out)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (allocator_callbacks && (!allocator_callbacks->alloc || !allocator_callbacks->free))
    {
        return JMTX_RESULT_NULL_PARAM;
    }

    return jmtxd_multiply_matrix_crs(a, b, p_out, allocator_callbacks);
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
jmtx_result jmtxd_multiply_matrix_ccs(
        const jmtxd_matrix_crs* a, const jmtxd_matrix_ccs* b, jmtxd_matrix_ccs** p_out,
        const jmtx_allocator_callbacks* allocator_callbacks)
{
    if (!allocator_callbacks)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }

    const uint32_t r_out = a->base.rows;
    const uint32_t c_out = b->base.cols;

    jmtxd_matrix_ccs* out;
    jmtx_result res = jmtxd_matrix_ccs_new(&out, c_out, r_out, a->n_entries + b->n_entries, allocator_callbacks);
    if (res != JMTX_RESULT_SUCCESS)
    {
        return res;
    }

    uint32_t capacity = 1 + a->n_entries / a->base.rows > 1 + b->n_entries / b->base.cols ?
                        a->n_entries / a->base.rows :
                        b->n_entries / b->base.cols;
    uint32_t count;
    uint32_t* indices = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*indices) * capacity);
    if (!indices)
    {
        jmtxd_matrix_ccs_destroy(out);
        return JMTX_RESULT_BAD_ALLOC;
    }
    double* values = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*values) * capacity);
    if (!values)
    {
        allocator_callbacks->free(allocator_callbacks->state, indices);
        jmtxd_matrix_ccs_destroy(out);
        return JMTX_RESULT_BAD_ALLOC;
    }


    //  Go column by column, so that the building of output matrix is done more efficiently
    for (uint32_t j = 0; j < c_out; ++j)
    {
        count = 0;
        uint32_t n_a, n_b;
        uint32_t* i_a, *i_b;
        double* v_a, *v_b;
        n_b = jmtxd_matrix_ccs_get_col(b, j, &i_b, &v_b);
        for (uint32_t i = 0; i < r_out; ++i)
        {
            n_a = jmtxd_matrix_crs_get_row(a, i, &i_a, &v_a);
            double v = 0;
            for (uint32_t k_a = 0, k_b = 0; k_a < n_a && k_b < n_b;)
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
                    //i_a[k_a] > i_b[k_b]
                    k_a += 1;
                }
            }
            if (v != 0)
            {
                if (count == capacity)
                {
                    const uint32_t new_capacity = capacity;
                    double* const new_v = allocator_callbacks->realloc(allocator_callbacks->state, values, sizeof(*values) * new_capacity);
                    if (!new_v)
                    {
                        allocator_callbacks->free(allocator_callbacks->state, values);
                        allocator_callbacks->free(allocator_callbacks->state, indices);
                        jmtxd_matrix_ccs_destroy(out);
                        return JMTX_RESULT_BAD_ALLOC;
                    }
                    values = new_v;
                    uint32_t* const new_i = allocator_callbacks->realloc(allocator_callbacks->state, indices, sizeof(*indices) * new_capacity);
                    if (!new_i)
                    {
                        allocator_callbacks->free(allocator_callbacks->state, values);
                        allocator_callbacks->free(allocator_callbacks->state, indices);
                        jmtxd_matrix_ccs_destroy(out);
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
        jmtxd_matrix_ccs_build_col(out, j, count, indices, values);
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
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result jmtxds_multiply_matrix_ccs(const jmtxd_matrix_crs* a, const jmtxd_matrix_ccs* b, jmtxd_matrix_ccs** p_out,
                                       const jmtx_allocator_callbacks* allocator_callbacks)
{
    if (!a)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (a->base.type != JMTXD_TYPE_CRS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (!b)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (b->base.type != JMTXD_TYPE_CCS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (a->base.cols != b->base.rows)
    {
        //  can't do multiplication
        return JMTX_RESULT_BAD_MATRIX;
    }
    if (!p_out)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (allocator_callbacks && (!allocator_callbacks->alloc || !allocator_callbacks->free))
    {
        return JMTX_RESULT_NULL_PARAM;
    }

    return jmtxd_multiply_matrix_ccs(a, b, p_out, allocator_callbacks);
}


double jmtxd_multiply_matrix_sparse_vectors(uint32_t n_a, const uint32_t i_a[const static n_a], const double v_a[const static n_a],
                                          uint32_t n_b, const uint32_t i_b[const static n_b], const double v_b[const static n_b])
{
    double v = 0;
    uint32_t ia = 0, ib = 0;
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
        else //if (i_a[ia] > i_b[ib])
        {
            ib += 1;
        }
    }
    return v;
}

double jmtxd_multiply_matrix_sparse_vectors_limit(uint32_t max_a, uint32_t max_b, uint32_t n_a,
                                                const uint32_t i_a[static n_a], const double v_a[static max_a],
                                                uint32_t n_b, const uint32_t i_b[static n_b],
                                                const double v_b[static max_b])
{
    double v = 0;
    uint32_t ia = 0, ib = 0;
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
        else //if (i_a[ia] > i_b[ib])
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
jmtx_result jmtxd_multiply_matrix_brm(
        const jmtxd_matrix_brm* a, const jmtxd_matrix_brm* b, jmtxd_matrix_brm** p_out,
        const jmtx_allocator_callbacks* allocator_callbacks)
{

    if (!allocator_callbacks)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }

    const uint32_t r_out = a->base.rows;
    const uint32_t c_out = b->base.cols;

    uint_fast32_t max_ubw;
//    uint_fast32_t min_ubw;
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

    uint_fast32_t max_lbw;
//    uint_fast32_t min_lbw;
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

    jmtxd_matrix_brm* out;
    jmtx_result res = jmtxd_matrix_brm_new(&out, c_out, r_out, max_ubw, max_lbw, NULL, allocator_callbacks);
    if (res != JMTX_RESULT_SUCCESS)
    {
        return res;
    }

    const uint32_t capacity = 1 + max_ubw + max_lbw;
    double* col_vals = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*col_vals) * capacity);
    if (!col_vals)
    {
        jmtxd_matrix_brm_destroy(out);
        return JMTX_RESULT_BAD_ALLOC;
    }
    double* values = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*col_vals) * capacity);
    if (!values)
    {
        allocator_callbacks->free(allocator_callbacks->state, col_vals);
        jmtxd_matrix_brm_destroy(out);
        return JMTX_RESULT_BAD_ALLOC;
    }


    //  Go row by row, so that the building of output matrix is done more efficiently
    for (uint32_t i = 0; i < r_out; ++i)
    {
        uint_fast32_t k = 0;
        uint32_t n_a, n_b;
//        uint32_t* i_a, *i_b;
        double* v_a;//, *v_b = col_vals;
        n_a = jmtxd_matrix_brm_get_row(a, i, &v_a);
        (void) n_a, (void)n_b;
        const uint_fast32_t first_a = jmtxd_matrix_brm_first_pos_in_row(a, i);
        const uint_fast32_t last_a = jmtxd_matrix_brm_last_pos_in_row(a, i);
        for (uint_fast32_t j = jmtxd_matrix_brm_first_pos_in_row(out, i);
             j < jmtxd_matrix_brm_last_pos_in_row(out, i) + 1;
             ++j, ++k)
        {
            n_b = jmtxd_matrix_brm_get_col(b, j, col_vals);
            const uint_fast32_t first_b = jmtxd_matrix_brm_first_pos_in_col(b, j);
            const uint_fast32_t last_b = jmtxd_matrix_brm_last_pos_in_col(b, j);

            uint_fast32_t first;

            uint_fast32_t off_a = 0;
            uint_fast32_t off_b = 0;

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

            double v = 0;
            for (int_fast32_t p = 0; p < len; ++p)
            {
                v += v_a[off_a + p] * col_vals[off_b + p];
            }
            values[k] = v;
        }
        
        jmtxd_matrix_brm_set_row(out, i, values);
    }


    allocator_callbacks->free(allocator_callbacks->state, values);
    allocator_callbacks->free(allocator_callbacks->state, col_vals);
    *p_out = out;
    return JMTX_RESULT_SUCCESS;
}

/**
 * Multiplies two BRM matrices together and produces a BRM matrix with the result of the matrix multiplication.
 * @param a BRM matrix
 * @param b BRM matrix
 * @param p_out pointer which receives the resulting BRM matrix
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result jmtxds_multiply_matrix_brm(const jmtxd_matrix_brm* a, const jmtxd_matrix_brm* b, jmtxd_matrix_brm** p_out,
                                       const jmtx_allocator_callbacks* allocator_callbacks)
{
    if (!a)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (a->base.type != JMTXD_TYPE_BRM)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (!b)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (b->base.type != JMTXD_TYPE_BRM)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (a->base.cols != b->base.rows)
    {
        //  can't do multiplication
        return JMTX_RESULT_BAD_MATRIX;
    }
    if (!p_out)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (allocator_callbacks && (!allocator_callbacks->alloc || !allocator_callbacks->free))

    {
        return JMTX_RESULT_NULL_PARAM;
    }

    return jmtxd_multiply_matrix_brm(a, b, p_out, allocator_callbacks);
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
jmtx_result jmtxd_multiply_matrix_cds(const jmtxd_matrix_cds* a, const jmtxd_matrix_cds* b, jmtxd_matrix_cds** p_out,
                                     const jmtx_allocator_callbacks* allocator_callbacks)
{

    if (!allocator_callbacks)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }

    const uint32_t r_out = a->base.rows;
    const uint32_t c_out = b->base.cols;

    const uint_fast32_t cnt_a = jmtxd_matrix_cds_diagonal_count(a);
    const uint_fast32_t cnt_b = jmtxd_matrix_cds_diagonal_count(b);

    uint32_t* const idx_a = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*idx_a) * cnt_a);
    if (!idx_a)
    {
        return JMTX_RESULT_BAD_ALLOC;
    }
    uint32_t* const idx_b = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*idx_b) * cnt_b);
    if (!idx_b)
    {
        allocator_callbacks->free(allocator_callbacks->state, idx_a);
        return JMTX_RESULT_BAD_ALLOC;
    }
    double* const val_a = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*val_a) * cnt_a);
    if (!val_a)
    {
        allocator_callbacks->free(allocator_callbacks->state, idx_b);
        allocator_callbacks->free(allocator_callbacks->state, idx_a);
        return JMTX_RESULT_BAD_ALLOC;
    }
    double* const val_b = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*val_b) * cnt_b);
    if (!val_b)
    {
        allocator_callbacks->free(allocator_callbacks->state, val_a);
        allocator_callbacks->free(allocator_callbacks->state, idx_b);
        allocator_callbacks->free(allocator_callbacks->state, idx_a);
        return JMTX_RESULT_BAD_ALLOC;
    }


    jmtxd_matrix_cds* out;
    jmtx_result res = jmtxd_matrix_cds_new(&out, c_out, r_out, 0, (int32_t[]){0}, allocator_callbacks);
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
    for (uint32_t i = 0; i < r_out; ++i)
    {
        const uint32_t c_a = jmtxd_matrix_cds_get_row(a, i, cnt_a, val_a, idx_a);
        for (uint_fast32_t j = 0; j < c_out; ++j)
        {
            const uint32_t c_b = jmtxd_matrix_cds_get_col(b, j, cnt_b, val_b, idx_b);
            const double val = jmtxd_multiply_matrix_sparse_vectors(c_a, idx_a, val_a, c_b, idx_b, val_b);
            if (val != 0)
            {
                res = jmtxd_matrix_cds_insert_entry(out, i, j, val);
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

/**
 * Multiplies two CDS matrices together and produces a CDS matrix with the result of the matrix multiplication.
 * @param a CDS matrix
 * @param b CDS matrix
 * @param p_out pointer which receives the resulting CDS matrix
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result jmtxds_multiply_matrix_cds(const jmtxd_matrix_cds* a, const jmtxd_matrix_cds* b, jmtxd_matrix_cds** p_out,
                                       const jmtx_allocator_callbacks* allocator_callbacks)
{
    if (!a)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (a->base.type != JMTXD_TYPE_CDS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (!b)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (b->base.type != JMTXD_TYPE_CDS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (a->base.cols != b->base.rows)
    {
        //  can't do multiplication
        return JMTX_RESULT_BAD_MATRIX;
    }
    if (!p_out)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (allocator_callbacks && (!allocator_callbacks->alloc || !allocator_callbacks->free))
    {
        return JMTX_RESULT_NULL_PARAM;
    }

    return jmtxd_multiply_matrix_cds(a, b, p_out, allocator_callbacks);
}
