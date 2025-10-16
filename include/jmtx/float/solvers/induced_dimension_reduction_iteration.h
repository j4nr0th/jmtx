//
// Created by jan on 9.1.2024.
//

#ifndef JMTX_INDUCED_DIMENSION_REDUCTION_ITERATION_H
#define JMTX_INDUCED_DIMENSION_REDUCTION_ITERATION_H
#ifndef JMTX_SOLVER_BASE_H
#    include "../../solver_base.h"
#endif

#include "../matrices/dense_row_major.h"
#include "../matrices/sparse_row_compressed.h"

jmtx_result jmtx_solve_iterative_idrs_crs(const jmtxf_matrix_crs *mtx, const float *restrict y, float *restrict x,
                                          float *restrict aux_vec1, float *restrict aux_vec2, uint32_t s,
                                          float *restrict aux_vec3, float *restrict aux_vec4, float *restrict aux_vec5,
                                          float *restrict aux_vec6, float *restrict aux_vec7,
                                          const jmtxf_matrix_drm *p_mtx, jmtxf_matrix_drm *aux_mtx1,
                                          jmtxf_matrix_drm *aux_mtx2, jmtxf_matrix_drm *aux_mtx3,
                                          jmtxf_matrix_drm *aux_mtx4, jmtxf_matrix_drm *aux_mtx5,
                                          jmtxf_solver_arguments *args);

#endif // JMTX_INDUCED_DIMENSION_REDUCTION_ITERATION_H
