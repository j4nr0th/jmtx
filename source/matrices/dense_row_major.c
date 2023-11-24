//
// Created by jan on 13.7.2023.
//
#include "matrix_base_internal.h"

#include "dense_row_major.h"

jmtx_result jmtx_matrix_drm_new(
        jmtx_matrix_drm** p_out, uint32_t rows, uint32_t cols, int zero,
        const jmtx_allocator_callbacks* allocator_callbacks)
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
    jmtx_matrix_drm* mtx = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*mtx));
    if (!mtx)
    {
        allocator_callbacks->free(allocator_callbacks->state, p_elements);
        return JMTX_RESULT_BAD_ALLOC;
    }
    if (zero)
    {
        memset(p_elements, 0, (rows * cols) *  sizeof(*p_elements));
    }
    mtx->base.allocator_callbacks = *allocator_callbacks;
    mtx->base.rows = rows;
    mtx->base.cols = cols;
    mtx->base.type = JMTX_TYPE_DRM;
    mtx->elements = p_elements;
    *p_out = mtx;
    return JMTX_RESULT_SUCCESS;
}

jmtx_result jmtx_matrix_drm_get_element(const jmtx_matrix_drm* mtx, uint32_t row, uint32_t col, float* x)
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
    *x = mtx->elements[col + row * mtx->base.cols];
    return JMTX_RESULT_SUCCESS;
}

jmtx_result jmtx_matrix_drm_set_element(jmtx_matrix_drm* mtx, uint32_t row, uint32_t col, float x)
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
    mtx->elements[col + row * mtx->base.cols] = x;
    return JMTX_RESULT_SUCCESS;
}

jmtx_result jmtx_matrix_drm_add_element(jmtx_matrix_drm* mtx, uint32_t row, uint32_t col, float x)
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
    mtx->elements[col + row * mtx->base.cols] += x;
    return JMTX_RESULT_SUCCESS;
}

jmtx_result jmtx_matrix_drm_destroy(jmtx_matrix_drm* mtx)
{
    mtx->base.allocator_callbacks.free(mtx->base.allocator_callbacks.state, mtx->elements);
    mtx->base.allocator_callbacks.free(mtx->base.allocator_callbacks.state, mtx);
    return JMTX_RESULT_SUCCESS;
}

jmtx_result jmtx_matrix_drm_set_col(jmtx_matrix_drm* mtx, uint32_t col, const float* elements)
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
    for (uint32_t i = 0; i < mtx->base.rows; ++i)
    {
        mtx->elements[col + i * mtx->base.cols] = elements[i];
    }

    return JMTX_RESULT_SUCCESS;
}

jmtx_result jmtx_matrix_drm_set_row(jmtx_matrix_drm* mtx, uint32_t row, const float* elements)
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
    memcpy(mtx->elements + mtx->base.cols * row, elements, sizeof(*elements) * mtx->base.cols);
    return JMTX_RESULT_SUCCESS;
}

jmtx_result jmtx_matrix_drm_set_all_rm(jmtx_matrix_drm* mtx, const float* elements)
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
    memcpy(mtx->elements , elements, sizeof(*elements) * mtx->base.cols * mtx->base.rows);
    return JMTX_RESULT_SUCCESS;
}

jmtx_result jmtx_matrix_drm_set_all_cm(jmtx_matrix_drm* mtx, const float* elements)
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
    for (uint32_t row = 0; row < mtx->base.rows; ++row)
    {
        for (uint32_t col = 0; col < mtx->base.cols; ++col)
        {
            mtx->elements[row * mtx->base.cols + col] = elements[col * mtx->base.rows + row];
        }
    }
    return JMTX_RESULT_SUCCESS;
}

jmtx_result jmtx_matrix_drm_transpose(jmtx_matrix_drm* mtx, jmtx_matrix_drm* out)
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
    if (out)
    {
        if (out->base.type != JMTX_TYPE_DRM)
        {
            return JMTX_RESULT_WRONG_TYPE;
        }
        if (out->base.cols != mtx->base.rows || out->base.rows != mtx->base.cols)
        {
            return JMTX_RESULT_DIMS_MISMATCH;
        }
    }
#endif
    if (!out)
    {
        out = mtx;
    }

    const float* const elements_in = mtx->elements;
    float* const elements_out = out->elements;

    for (uint32_t row = 0; row < mtx->base.rows - 1; ++row)
    {
        for (uint32_t col = row + 1; col < mtx->base.cols; ++col)
        {
            const float tmp = elements_in[row + col * mtx->base.cols];
            elements_out[row + col * mtx->base.cols] = elements_in[col + row * mtx->base.cols];
            elements_out[col + row * mtx->base.cols] = tmp;
        }
    }
    const uint32_t tmp = mtx->base.rows;
    out->base.rows = mtx->base.cols;
    out->base.cols = tmp;

    return JMTX_RESULT_SUCCESS;
}

jmtx_result jmtx_matrix_drm_elements_in_row(jmtx_matrix_drm* mtx, uint32_t* p_count)
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

jmtx_result jmtx_matrix_drm_elements_in_col(jmtx_matrix_drm* mtx, uint32_t* p_count)
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

jmtx_result jmtx_matrix_drm_get_row(jmtx_matrix_drm* mtx, uint32_t row, float* p_out)
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
    memcpy(p_out, mtx->elements + mtx->base.cols * row, sizeof(*mtx->elements) * mtx->base.cols);
    return JMTX_RESULT_SUCCESS;
}

jmtx_result jmtx_matrix_drm_get_col(jmtx_matrix_drm* mtx, uint32_t col, float* p_out)
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
    for (uint32_t i = 0; i < mtx->base.rows; ++i)
    {
        p_out[i] = mtx->elements[col + mtx->base.cols * i];
    }
    return JMTX_RESULT_SUCCESS;
}

jmtx_result jmtx_matrix_drm_multiply(jmtx_matrix_drm* out, const jmtx_matrix_drm* first, const jmtx_matrix_drm* second)
{
#ifndef JMTX_NO_VERIFY_PARAMS
    if (!first)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (first->base.type != JMTX_TYPE_DRM)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (!second)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (second->base.type != JMTX_TYPE_DRM)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (first->base.cols != second->base.rows)
    {
        return JMTX_RESULT_DIMS_MISMATCH;
    }
    if (!out)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (out->base.type != JMTX_TYPE_DRM)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (out->base.cols != second->base.cols)
    {
        return JMTX_RESULT_DIMS_MISMATCH;
    }
    if (out->base.rows != first->base.rows)
    {
        return JMTX_RESULT_DIMS_MISMATCH;
    }
    //  No duplicates allowed
    if (out == second || out == first)
    {
        return JMTX_RESULT_BAD_PARAM;
    }
#endif

    float* const c = out->elements;
    const float* const a = first->elements;
    const float* const b = second->elements;

    for (uint32_t row = 0; row < out->base.rows; ++row)
    {
        for (uint32_t col = 0; col < out->base.cols; ++col)
        {
            float res = 0.0f;
            for (uint32_t k = 0; k < first->base.cols; ++row) //  first->base.cols == second->base.rows
            {
                res += a[k + first->base.cols * row] * b[col + second->base.cols * k];
            }
            c[col + out->base.cols * row] = res;
        }
    }

    return JMTX_RESULT_SUCCESS;
}
