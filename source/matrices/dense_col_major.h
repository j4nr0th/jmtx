//
// Created by jan on 13.7.2023.
//

#ifndef JMTX_DENSE_COL_MAJOR_H
#define JMTX_DENSE_COL_MAJOR_H
#ifndef JMTX_MATRIX_BASE_H
    #include "../../include/jmtx/matrices/matrix_base.h"
#endif

#ifndef JMTX_MATRIX_BASE_INTERNAL_H
    #include "matrix_base_internal.h"
#endif

typedef struct jmtx_matrix_dcm_struct jmtx_matrix_dcm;
struct jmtx_matrix_dcm_struct
{
    jmtx_matrix_base base;
    float* elements;
};

jmtx_result jmtx_matrix_dcm_new(
        jmtx_matrix_dcm** p_out, uint32_t rows, uint32_t cols, int zero,
        const jmtx_allocator_callbacks* memory_allocators);

jmtx_result jmtx_matrix_dcm_get_element(const jmtx_matrix_dcm* mtx, uint32_t row, uint32_t col, float* x);

jmtx_result jmtx_matrix_dcm_set_element(jmtx_matrix_dcm* mtx, uint32_t row, uint32_t col, float x);

jmtx_result jmtx_matrix_dcm_destroy(jmtx_matrix_dcm* mtx);

jmtx_result jmtx_matrix_dcm_set_col(jmtx_matrix_dcm* mtx, uint32_t col, const float* elements);

jmtx_result jmtx_matrix_dcm_set_row(jmtx_matrix_dcm* mtx, uint32_t row, const float* elements);

jmtx_result jmtx_matrix_dcm_set_rm(jmtx_matrix_dcm* mtx, float* elements);

jmtx_result jmtx_matrix_dcm_set_cm(jmtx_matrix_dcm* mtx, float* elements);

jmtx_result jmtx_matrix_dcm_transpose(jmtx_matrix_dcm* mtx, jmtx_matrix_dcm* out);

jmtx_result jmtx_matrix_dcm_elements_in_row(jmtx_matrix_dcm* mtx, uint32_t* p_count);

jmtx_result jmtx_matrix_dcm_elements_in_col(jmtx_matrix_dcm* mtx, uint32_t* p_count);

jmtx_result jmtx_matrix_dcm_get_row(jmtx_matrix_dcm* mtx, uint32_t row, float* p_out);

jmtx_result jmtx_matrix_dcm_get_col(jmtx_matrix_dcm* mtx, uint32_t col, float* p_out);

jmtx_result
jmtx_matrix_dcm_multiply(jmtx_matrix_dcm* out, const jmtx_matrix_dcm* first, const jmtx_matrix_dcm* second);

#endif //JMTX_DENSE_COL_MAJOR_H
