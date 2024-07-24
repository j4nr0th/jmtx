//
// Created by jan on 9.1.2024.
//

#ifndef JMTX_INDUCED_DIMENSION_REDUCTION_ITERATION_H
#define JMTX_INDUCED_DIMENSION_REDUCTION_ITERATION_H
#ifndef JMTX_SOLVER_BASE_H
    #include "../../solver_base.h"
#endif

#include "../matrices/sparse_row_compressed.h"
#include "../matrices/dense_row_major.h"

jmtx_result jmtx_solve_iterative_idrs_crs(const jmtx_matrix_crs* mtx, const float* restrict y, float* restrict x,
                                          float* restrict aux_vec1, float* restrict aux_vec2, uint32_t s,
                                          float* restrict aux_vec3, float* restrict aux_vec4,
                                          float* restrict aux_vec5, float* restrict aux_vec6,
                                          float* restrict aux_vec7, const jmtx_matrix_drm* p_mtx,
                                          jmtx_matrix_drm* aux_mtx1, jmtx_matrix_drm* aux_mtx2,
                                          jmtx_matrix_drm* aux_mtx3, jmtx_matrix_drm* aux_mtx4,
                                          jmtx_matrix_drm* aux_mtx5, jmtx_solver_arguments* args);

#endif //JMTX_INDUCED_DIMENSION_REDUCTION_ITERATION_H
