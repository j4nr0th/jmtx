//
// Created by jan on 13.7.2023.
//

#ifndef JMTX_DENSE_ROW_MAJOR_H
#define JMTX_DENSE_ROW_MAJOR_H
#include "matrix_base.h"

typedef struct jmtx_matrix_drm_struct jmtx_matrix_drm;
struct jmtx_matrix_drm_struct
{
    jmtx_matrix base;
    float* elements;
};

jmtx_result jmtx_matrix_drm_new(jmtx_matrix_drm** p_out, uint32_t rows, uint32_t cols, const jmtx_allocator_callbacks* allocator_callbacks);

jmtx_result jmtx_matrix_drm_get_element(const jmtx_matrix_drm* mtx, uint32_t row, uint32_t col, float* x);

jmtx_result jmtx_matrix_drm_set_element(jmtx_matrix_drm* mtx, uint32_t row, uint32_t col, float x);

jmtx_result jmtx_matrix_drm_add_element(jmtx_matrix_drm* mtx, uint32_t row, uint32_t col, float x);

jmtx_result jmtx_matrix_drm_destroy(jmtx_matrix_drm* mtx);

jmtx_result jmtx_matrix_drm_set_col(jmtx_matrix_drm* mtx, uint32_t col, const float* elements);

jmtx_result jmtx_matrix_drm_set_row(jmtx_matrix_drm* mtx, uint32_t row, const float* elements);

jmtx_result jmtx_matrix_drm_set_all_rm(jmtx_matrix_drm* mtx, const float* elements);

jmtx_result jmtx_matrix_drm_set_all_cm(jmtx_matrix_drm* mtx, const float* elements);

jmtx_result jmtx_matrix_drm_transpose(jmtx_matrix_drm* mtx, jmtx_matrix_drm** p_out);

jmtx_result jmtx_matrix_drm_elements_in_row(jmtx_matrix_drm* mtx, uint32_t* p_count);

jmtx_result jmtx_matrix_drm_elements_in_col(jmtx_matrix_drm* mtx, uint32_t* p_count);

jmtx_result jmtx_matrix_drm_get_row(jmtx_matrix_drm* mtx, uint32_t row, float* p_out);

jmtx_result jmtx_matrix_drm_get_col(jmtx_matrix_drm* mtx, uint32_t col, float* p_out);

#endif //JMTX_DENSE_ROW_MAJOR_H
