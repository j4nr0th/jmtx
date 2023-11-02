//
// Created by jan on 2.11.2023.
//

#include <assert.h>
#include "sparse_conversion.h"
#include "sparse_row_compressed_internal.h"
#include "sparse_column_compressed_internal.h"


jmtx_result jmtx_convert_crs_to_ccs_inplace_transpose(jmtx_matrix_crs* in, jmtx_matrix_ccs** p_out)
{
#ifndef JMTX_NO_VERIFY_PARAMS
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
#endif

    //  Lol, lmao even
    in->base.type = JMTX_TYPE_CCS;
    in->base.get_element = (jmtx_result (*)(
            jmtx_matrix*, jmtx_index_t, jmtx_index_t, float*)) jmtx_matrix_ccs_get_entry;
    *p_out = (jmtx_matrix_ccs*)in;

    return JMTX_RESULT_SUCCESS;
}

jmtx_result jmtx_convert_ccs_to_crs_inplace_transpose(jmtx_matrix_ccs* in, jmtx_matrix_crs** p_out)
{
#ifndef JMTX_NO_VERIFY_PARAMS
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
#endif

    //  Lol, lmao even
    in->base.type = JMTX_TYPE_CRS;
    in->base.get_element = (jmtx_result (*)(
            jmtx_matrix*, jmtx_index_t, jmtx_index_t, float*)) jmtx_matrix_crs_get_entry;
    *p_out = (jmtx_matrix_crs*)in;

    return JMTX_RESULT_SUCCESS;
}

jmtx_result jmtx_convert_crs_to_ccs(
        const jmtx_matrix_crs* in, jmtx_matrix_ccs** p_out, const jmtx_allocator_callbacks* allocator_callbacks)
{
    jmtx_matrix_crs* cpy;
    jmtx_result res = jmtx_matrix_crs_copy(in, &cpy, allocator_callbacks);
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
    jmtx_result res = jmtx_matrix_ccs_copy(in, &cpy, allocator_callbacks);
    if (res != JMTX_RESULT_SUCCESS)
    {
        return res;
    }
    res = jmtx_convert_ccs_to_crs_inplace_transpose(cpy, p_out);
    assert(res == JMTX_RESULT_SUCCESS);
    return JMTX_RESULT_SUCCESS;
}
