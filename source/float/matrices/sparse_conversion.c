//
// Created by jan on 2.11.2023.
//

#include <assert.h>
#include "../../../include/jmtx/float/matrices/sparse_conversion.h"
#include "sparse_row_compressed_internal.h"
#include "sparse_column_compressed_internal.h"


jmtx_result jmtx_convert_crs_to_ccs_inplace_transpose(jmtx_matrix_crs* in, jmtx_matrix_ccs** p_out)
{
    if (!in)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (in->base.type != JMTX_TYPE_CRS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (!p_out)
    {
        return JMTX_RESULT_NULL_PARAM;
    }

    //  Lol, lmao even
    in->base.type = JMTX_TYPE_CCS;

    *p_out = (jmtx_matrix_ccs*)in;

    return JMTX_RESULT_SUCCESS;
}

jmtx_result jmtx_convert_ccs_to_crs_inplace_transpose(jmtx_matrix_ccs* in, jmtx_matrix_crs** p_out)
{
    if (!in)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (in->base.type != JMTX_TYPE_CCS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (!p_out)
    {
        return JMTX_RESULT_NULL_PARAM;
    }


    //  Lol, lmao even
    in->base.type = JMTX_TYPE_CRS;
    *p_out = (jmtx_matrix_crs*)in;

    return JMTX_RESULT_SUCCESS;
}

jmtx_result jmtx_convert_crs_to_ccs(
        const jmtx_matrix_crs* in, jmtx_matrix_ccs** p_out, const jmtx_allocator_callbacks* allocator_callbacks)
{
    jmtx_matrix_crs* cpy;
    jmtx_result res = jmtx_matrix_crs_transpose(in, &cpy, allocator_callbacks);
    if (res != JMTX_RESULT_SUCCESS)
    {
        return res;
    }
    res = jmtx_convert_crs_to_ccs_inplace_transpose(cpy, p_out);
    assert(res == JMTX_RESULT_SUCCESS);
    return JMTX_RESULT_SUCCESS;
}

jmtx_result jmtx_convert_ccs_to_crs(
        const jmtx_matrix_ccs* in, jmtx_matrix_crs** p_out, const jmtx_allocator_callbacks* allocator_callbacks)
{
    jmtx_matrix_ccs* cpy;
    jmtx_result res = jmtx_matrix_ccs_transpose(in, &cpy, allocator_callbacks);
    if (res != JMTX_RESULT_SUCCESS)
    {
        return res;
    }
    res = jmtx_convert_ccs_to_crs_inplace_transpose(cpy, p_out);
    assert(res == JMTX_RESULT_SUCCESS);
    return JMTX_RESULT_SUCCESS;
//    jmtx_matrix_crs* out;
//    if (!allocator_callbacks)
//    {
//        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
//    }
//    jmtx_result res = jmtx_matrix_crs_new(&out, in->base.cols, in->base.rows, in->n_entries, allocator_callbacks);
//    if (res != JMTX_RESULT_SUCCESS)
//    {
//        return res;
//    }
//    uint32_t capacity = in->n_entries / in->base.rows + 1;
//    float* values = allocator_callbacks->alloc(allocator_callbacks->state, capacity * sizeof(*values));
//    uint32_t* indices = allocator_callbacks->alloc(allocator_callbacks->state, capacity * sizeof(*indices));
//    if (!values || !indices)
//    {
//        goto failed_alloc;
//    }
//    for (uint32_t i = 0; i < in->base.rows; ++i)
//    {
//        uint32_t count;
//        jmtx_matrix_ccs_elements_in_row(in, i, &count);
//        if (count > capacity)
//        {
//            const uint32_t new_capacity = count;
//            float* const new_f = allocator_callbacks->alloc(allocator_callbacks->state, new_capacity * sizeof(*new_f));
//            if (!new_f)
//            {
//                goto failed_alloc;
//            }
//            values = new_f;
//            uint32_t* const new_u = allocator_callbacks->alloc(allocator_callbacks->state, new_capacity * sizeof(*new_u));
//            if (!new_u)
//            {
//                goto failed_alloc;
//            }
//            indices = new_u;
//            capacity = new_capacity;
//        }
//        uint32_t real_count;
//        jmtx_matrix_ccs_get_row(in, i, capacity, values, &real_count, indices);
//        jmtx_matrix_crs_build_row(out, i, real_count, indices, values);
//    }
//
//    allocator_callbacks->free(allocator_callbacks->state, indices);
//    allocator_callbacks->free(allocator_callbacks->state, values);
//
//    *p_out = out;
//    return JMTX_RESULT_SUCCESS;
//
//failed_alloc:
//    allocator_callbacks->free(allocator_callbacks->state, indices);
//    allocator_callbacks->free(allocator_callbacks->state, values);
//    jmtx_matrix_crs_destroy(out);
//    return JMTX_RESULT_BAD_ALLOC;
}
