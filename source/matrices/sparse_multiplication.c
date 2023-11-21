//
// Created by jan on 2.11.2023.
//

#include "../../include/jmtx/matrices/sparse_multiplication.h"
#include "sparse_column_compressed_internal.h"
#include "sparse_row_compressed_internal.h"

jmtx_result jmtx_matrix_multiply_crs(
        const jmtx_matrix_crs* a, const jmtx_matrix_ccs* b, jmtx_matrix_crs** p_out,
        const jmtx_allocator_callbacks* allocator_callbacks)
{
    if (!a)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (a->base.type != JMTX_TYPE_CRS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (!b)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (b->base.type != JMTX_TYPE_CCS)
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

    if (!allocator_callbacks)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }
    else if (!allocator_callbacks->alloc || !allocator_callbacks->free)
    {
        return JMTX_RESULT_NULL_PARAM;
    }

    const uint32_t r_out = a->base.rows;
    const uint32_t c_out = b->base.cols;

    jmtx_matrix_crs* out;
    jmtx_result res = jmtx_matrix_crs_new(&out, c_out, r_out, a->n_entries + b->n_entries, allocator_callbacks);
    if (res != JMTX_RESULT_SUCCESS)
    {
        return res;
    }

    uint32_t capacity = 1 + a->n_entries / a->base.rows > 1 + b->n_entries / b->base.cols ?
                               a->n_entries / a->base.rows :
                               b->n_entries / b->base.cols;
    uint32_t count = 0;
    uint32_t* indices = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*indices) * capacity);
    if (!indices)
    {
        jmtx_matrix_crs_destroy(out);
        return JMTX_RESULT_BAD_ALLOC;
    }
    float* values = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*indices) * capacity);
    if (!values)
    {
        allocator_callbacks->free(allocator_callbacks->state, indices);
        jmtx_matrix_crs_destroy(out);
        return JMTX_RESULT_BAD_ALLOC;
    }


    //  Go row by row, so that the building of output matrix is done more efficiently
    for (uint32_t i = 0; i < r_out; ++i)
    {
        count = 0;
        uint32_t n_a, n_b;
        uint32_t* i_a, *i_b;
        float* v_a, *v_b;
        n_a = jmtx_matrix_crs_get_row(a, i, &i_a, &v_a);
        for (uint32_t j = 0; j < c_out; ++j)
        {
            n_b = jmtx_matrix_ccs_get_col(b, j, &i_b, &v_b);
            float v = 0;
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
                    const uint32_t new_capacity = capacity << 1;
                    float* const new_v = allocator_callbacks->realloc(allocator_callbacks->state, values, sizeof(*values) * new_capacity);
                    if (!new_v)
                    {
                        allocator_callbacks->free(allocator_callbacks->state, values);
                        allocator_callbacks->free(allocator_callbacks->state, indices);
                        jmtx_matrix_crs_destroy(out);
                        return JMTX_RESULT_BAD_ALLOC;
                    }
                    values = new_v;
                    uint32_t* const new_i = allocator_callbacks->realloc(allocator_callbacks->state, indices, sizeof(*indices) * new_capacity);
                    if (!new_i)
                    {
                        allocator_callbacks->free(allocator_callbacks->state, values);
                        allocator_callbacks->free(allocator_callbacks->state, indices);
                        jmtx_matrix_crs_destroy(out);
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
        jmtx_matrix_crs_build_row(out, i, count, indices, values);
    }


    allocator_callbacks->free(allocator_callbacks->state, values);
    allocator_callbacks->free(allocator_callbacks->state, indices);
    *p_out = out;
    return JMTX_RESULT_SUCCESS;
}

jmtx_result jmtx_matrix_multiply_ccs(
        const jmtx_matrix_crs* a, const jmtx_matrix_ccs* b, jmtx_matrix_ccs** p_out,
        const jmtx_allocator_callbacks* allocator_callbacks)
{
    if (!a)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (a->base.type != JMTX_TYPE_CRS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (!b)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (b->base.type != JMTX_TYPE_CCS)
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

    if (!allocator_callbacks)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }
    else if (!allocator_callbacks->alloc || !allocator_callbacks->free)
    {
        return JMTX_RESULT_NULL_PARAM;
    }

    const uint32_t r_out = a->base.rows;
    const uint32_t c_out = b->base.cols;

    jmtx_matrix_ccs* out;
    jmtx_result res = jmtx_matrix_ccs_new(&out, c_out, r_out, a->n_entries + b->n_entries, allocator_callbacks);
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
        jmtx_matrix_ccs_destroy(out);
        return JMTX_RESULT_BAD_ALLOC;
    }
    float* values = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*indices) * capacity);
    if (!values)
    {
        allocator_callbacks->free(allocator_callbacks->state, indices);
        jmtx_matrix_ccs_destroy(out);
        return JMTX_RESULT_BAD_ALLOC;
    }


    //  Go column by column, so that the building of output matrix is done more efficiently
    for (uint32_t j = 0; j < c_out; ++j)
    {
        count = 0;
        uint32_t n_a, n_b;
        uint32_t* i_a, *i_b;
        float* v_a, *v_b;
        n_b = jmtx_matrix_ccs_get_col(b, j, &i_b, &v_b);
        for (uint32_t i = 0; i < r_out; ++i)
        {
            n_a = jmtx_matrix_crs_get_row(a, i, &i_a, &v_a);
            float v = 0;
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
                    float* const new_v = allocator_callbacks->realloc(allocator_callbacks->state, values, sizeof(*values) * new_capacity);
                    if (!new_v)
                    {
                        allocator_callbacks->free(allocator_callbacks->state, values);
                        allocator_callbacks->free(allocator_callbacks->state, indices);
                        jmtx_matrix_ccs_destroy(out);
                        return JMTX_RESULT_BAD_ALLOC;
                    }
                    values = new_v;
                    uint32_t* const new_i = allocator_callbacks->realloc(allocator_callbacks->state, indices, sizeof(*indices) * new_capacity);
                    if (!new_i)
                    {
                        allocator_callbacks->free(allocator_callbacks->state, values);
                        allocator_callbacks->free(allocator_callbacks->state, indices);
                        jmtx_matrix_ccs_destroy(out);
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
        jmtx_matrix_ccs_build_col(out, j, count, indices, values);
    }


    allocator_callbacks->free(allocator_callbacks->state, values);
    allocator_callbacks->free(allocator_callbacks->state, indices);
    *p_out = out;
    return JMTX_RESULT_SUCCESS;
}
