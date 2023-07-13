//
// Created by jan on 13.7.2023.
//

#include "dense_col_major.h"

jmtx_result jmtx_matrix_crm_new(
        jmtx_matrix_crm** p_out, uint32_t rows, uint32_t cols, const jmtx_allocator_callbacks* allocator_callbacks)
{
#ifndef JMTX_NO_VERIFY_PARAMS
    if (!p_out)
    {
//        REPORT_ERROR_MESSAGE("Matrix pointer was null");
//        LEAVE_FUNCTION();
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!rows)
    {
//        REPORT_ERROR_MESSAGE("Number of rows was 0");
//        LEAVE_FUNCTION();
        return JMTX_RESULT_BAD_PARAM;
    }
    if (!cols)
    {
//        REPORT_ERROR_MESSAGE("Number of columns was 0");
//        LEAVE_FUNCTION();
        return JMTX_RESULT_BAD_PARAM;
    }
#endif
    if (!allocator_callbacks)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }
    else if (!allocator_callbacks->alloc || !allocator_callbacks->realloc || !allocator_callbacks->free)
    {
        return JMTX_RESULT_BAD_PARAM;
    }

    float* p_elements = allocator_callbacks->alloc(allocator_callbacks->state, (cols * rows) *  sizeof(*p_elements));
    if (!p_elements)
    {
        return JMTX_RESULT_BAD_ALLOC;
    }
    jmtx_matrix_crm* mtx = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*mtx));
    if (!mtx)
    {
        allocator_callbacks->free(allocator_callbacks->state, p_elements);
        return JMTX_RESULT_BAD_ALLOC;
    }
    memset(p_elements, 0, (rows * cols) *  sizeof(*p_elements));
    mtx->base.allocator_callbacks = *allocator_callbacks;
    mtx->base.rows = rows;
    mtx->base.cols = cols;
    mtx->base.type = JMTX_TYPE_DRM;
    mtx->elements = p_elements;
    *p_out = mtx;
    return JMTX_RESULT_SUCCESS;
}

jmtx_result jmtx_matrix_crm_get_element(const jmtx_matrix_crm* mtx, uint32_t row, uint32_t col, float* x)
{
#ifndef JMTX_NO_VERIFY_PARAMS
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTX_TYPE_DRM)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (row >= mtx->base.rows)
    {
        return JMTX_RESULT_INDEX_OUT_OF_BOUNDS;
    }
    if (col >= mtx->base.cols)
    {
        return JMTX_RESULT_INDEX_OUT_OF_BOUNDS;
    }
    if (!x)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
#endif
    *x = mtx->elements[row + col * mtx->base.rows];
    return JMTX_RESULT_SUCCESS;
}

jmtx_result jmtx_matrix_crm_set_element(jmtx_matrix_crm* mtx, uint32_t row, uint32_t col, float x)
{
#ifndef JMTX_NO_VERIFY_PARAMS
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTX_TYPE_DRM)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (row >= mtx->base.rows)
    {
        return JMTX_RESULT_INDEX_OUT_OF_BOUNDS;
    }
    if (col >= mtx->base.cols)
    {
        return JMTX_RESULT_INDEX_OUT_OF_BOUNDS;
    }
#endif
    mtx->elements[row + col * mtx->base.rows] = x;
    return JMTX_RESULT_SUCCESS;
}

jmtx_result jmtx_matrix_crm_add_element(jmtx_matrix_crm* mtx, uint32_t row, uint32_t col, float x)
{
#ifndef JMTX_NO_VERIFY_PARAMS
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTX_TYPE_DRM)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (row >= mtx->base.rows)
    {
        return JMTX_RESULT_INDEX_OUT_OF_BOUNDS;
    }
    if (col >= mtx->base.cols)
    {
        return JMTX_RESULT_INDEX_OUT_OF_BOUNDS;
    }
#endif
    mtx->elements[row + col * mtx->base.rows] += x;
    return JMTX_RESULT_SUCCESS;
}

jmtx_result jmtx_matrix_crm_destroy(jmtx_matrix_crm* mtx)
{
    mtx->base.allocator_callbacks.free(mtx->base.allocator_callbacks.state, mtx->elements);
    mtx->base.allocator_callbacks.free(mtx->base.allocator_callbacks.state, mtx);
    return JMTX_RESULT_SUCCESS;
}

jmtx_result jmtx_matrix_crm_set_col(jmtx_matrix_crm* mtx, uint32_t col, const float* elements)
{
#ifndef JMTX_NO_VERIFY_PARAMS
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTX_TYPE_DRM)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (col >= mtx->base.cols)
    {
        return JMTX_RESULT_INDEX_OUT_OF_BOUNDS;
    }
    if (!elements)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
#endif
    memcpy(mtx->elements + mtx->base.rows * col, elements, sizeof(*elements) * mtx->base.rows);

    return JMTX_RESULT_SUCCESS;
}

jmtx_result jmtx_matrix_crm_set_row(jmtx_matrix_crm* mtx, uint32_t row, const float* elements)
{
#ifndef JMTX_NO_VERIFY_PARAMS
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTX_TYPE_DRM)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (row >= mtx->base.rows)
    {
        return JMTX_RESULT_INDEX_OUT_OF_BOUNDS;
    }
    if (!elements)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
#endif
    for (uint32_t i = 0; i < mtx->base.cols; ++i)
    {
        mtx->elements[row + i * mtx->base.rows] = elements[i];
    }
    return JMTX_RESULT_SUCCESS;
}

jmtx_result jmtx_matrix_crm_set_all_rm(jmtx_matrix_crm* mtx, const float* elements)
{
#ifndef JMTX_NO_VERIFY_PARAMS
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTX_TYPE_DRM)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (!elements)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
#endif
    for (uint32_t col = 0; col < mtx->base.cols; ++col)
    {
        for (uint32_t row = 0; row < mtx->base.rows; ++row)
        {
            mtx->elements[col * mtx->base.rows + row] = elements[row * mtx->base.cols + col];
        }
    }
    return JMTX_RESULT_SUCCESS;
}

jmtx_result jmtx_matrix_crm_set_all_cm(jmtx_matrix_crm* mtx, const float* elements)
{
#ifndef JMTX_NO_VERIFY_PARAMS
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTX_TYPE_DRM)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (!elements)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
#endif
    memcpy(mtx->elements , elements, sizeof(*elements) * mtx->base.rows * mtx->base.cols);
    return JMTX_RESULT_SUCCESS;
}

jmtx_result jmtx_matrix_crm_transpose(jmtx_matrix_crm* mtx, jmtx_matrix_crm** p_out)
{
#ifndef JMTX_NO_VERIFY_PARAMS
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTX_TYPE_DRM)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
#endif
    if (!p_out)
    {
        //  Transpose in place
        float* const elements = mtx->elements;
        if (mtx->base.rows == mtx->base.cols)
        {
            //  Transpose for a square matrix
            for (uint32_t row = 0; row < mtx->base.rows - 1; ++row)
            {
                for (uint32_t col = row + 1; col < mtx->base.cols; ++col)
                {
                    const float tmp = elements[row + col * mtx->base.cols];
                    elements[row + col * mtx->base.cols] = elements[row + col * mtx->base.rows];
                    elements[row + col * mtx->base.rows] = tmp;
                }
            }
        }
        else
        {
            //  Transpose for a non-square matrix
            const jmtx_allocator_callbacks* allocator_callbacks = &mtx->base.allocator_callbacks;
            float* const new_elements = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*new_elements) * mtx->base.rows * mtx->base.cols);
            if (!new_elements)
            {
                return JMTX_RESULT_BAD_ALLOC;
            }
            for (uint32_t row = 0; row < mtx->base.rows; ++row)
            {
                for (uint32_t col = 0; col < mtx->base.cols; ++col)
                {
                    new_elements[row + mtx->base.cols * col] = elements[col + mtx->base.cols * row];
                }
            }
            allocator_callbacks->free(allocator_callbacks->state, mtx->elements);
            mtx->elements = new_elements;
            const uint32_t tmp = mtx->base.rows;
            mtx->base.rows = mtx->base.cols;
            mtx->base.cols = tmp;
        }
    }
    else
    {
        const jmtx_allocator_callbacks* allocator_callbacks = &mtx->base.allocator_callbacks;
        jmtx_matrix_crm* out = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*out));
        if (!out)
        {
            return JMTX_RESULT_BAD_ALLOC;
        }
        float* const elements = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*elements) * mtx->base.rows * mtx->base.cols);
        if (!elements)
        {
            allocator_callbacks->free(allocator_callbacks->state, out);
            return JMTX_RESULT_BAD_ALLOC;
        }
        out->base = mtx->base;
        out->base.cols = mtx->base.rows;
        out->base.rows = mtx->base.cols;
        out->elements = elements;
        for (uint32_t row = 0; row < mtx->base.rows; ++row)
        {
            for (uint32_t col = 0; col < mtx->base.cols; ++col)
            {
                out->elements[row * mtx->base.cols + col] = mtx->elements[col * mtx->base.rows + row];
            }
        }
        *p_out = out;
    }

    return JMTX_RESULT_SUCCESS;
}

jmtx_result jmtx_matrix_crm_elements_in_row(jmtx_matrix_crm* mtx, uint32_t* p_count)
{
#ifndef JMTX_NO_VERIFY_PARAMS
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTX_TYPE_DRM)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (!p_count)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
#endif
    *p_count = mtx->base.cols;
    return JMTX_RESULT_SUCCESS;
}

jmtx_result jmtx_matrix_crm_elements_in_col(jmtx_matrix_crm* mtx, uint32_t* p_count)
{
#ifndef JMTX_NO_VERIFY_PARAMS
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTX_TYPE_DRM)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (!p_count)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
#endif
    *p_count = mtx->base.rows;
    return JMTX_RESULT_SUCCESS;
}

jmtx_result jmtx_matrix_crm_get_row(jmtx_matrix_crm* mtx, uint32_t row, float* p_out)
{
#ifndef JMTX_NO_VERIFY_PARAMS
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTX_TYPE_DRM)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (row >= mtx->base.rows)
    {
        return JMTX_RESULT_INDEX_OUT_OF_BOUNDS;
    }
    if (!p_out)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
#endif
    for (uint32_t i = 0; i < mtx->base.cols; ++i)
    {
        p_out[i] = mtx->elements[row + mtx->base.rows * i];
    }
    return JMTX_RESULT_SUCCESS;
}

jmtx_result jmtx_matrix_crm_get_col(jmtx_matrix_crm* mtx, uint32_t col, float* p_out)
{
#ifndef JMTX_NO_VERIFY_PARAMS
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTX_TYPE_DRM)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (col >= mtx->base.cols)
    {
        return JMTX_RESULT_INDEX_OUT_OF_BOUNDS;
    }
    if (!p_out)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
#endif
    memcpy(p_out, mtx->elements + mtx->base.rows * col, sizeof(*mtx->elements) * mtx->base.rows);
    return JMTX_RESULT_SUCCESS;
}
