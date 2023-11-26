//
// Created by jan on 6.11.2023.
//

#ifndef JMTX_LU_SOLVING_H
#define JMTX_LU_SOLVING_H
#ifndef JMTX_BAND_ROW_MAJOR_H
    #include "../matrices/band_row_major.h"
#endif
#ifndef JMTX_SPARSE_ROW_COMPRESSED_H
    #include "../matrices/sparse_row_compressed.h"
#endif
#ifndef JMTX_SOLVER_BASE_H
    #include "solver_base.h"
#endif

void jmtx_lu_solve_crs(const jmtx_matrix_crs* l, const jmtx_matrix_crs* u, const float* restrict y, float* restrict x);

void jmtx_lu_solve_inplace_crs(const jmtx_matrix_crs* l, const jmtx_matrix_crs* u, float* restrict x);

void jmtx_lu_solve_brm(const jmtx_matrix_brm* l, const jmtx_matrix_brm* u, const float* restrict y, float* restrict x);

void jmtx_lu_solve_inplace_brm(const jmtx_matrix_brm* l, const jmtx_matrix_brm* u, float* restrict x);


jmtx_result jmtx_incomplete_lu_decomposition_solve_crs(
        const jmtx_matrix_crs* mtx, const float* y, float* x, float* aux_vec, jmtx_solver_arguments* args,
        const jmtx_allocator_callbacks* allocator_callbacks);

jmtx_result jmtx_incomplete_lu_decomposition_solve_precomputed_crs(
        const jmtx_matrix_crs* mtx, const jmtx_matrix_crs* l, const jmtx_matrix_crs* u, const float* y, float* x,
        float* aux_vec, jmtx_solver_arguments* args);

jmtx_result jmtx_incomplete_lu_decomposition_solve_crs_parallel(
        const jmtx_matrix_crs* mtx, const float* y, float* x, float* aux_vec, jmtx_solver_arguments* args,
        const jmtx_allocator_callbacks* allocator_callbacks);

jmtx_result jmtx_incomplete_lu_decomposition_solve_precomputed_crs_parallel(
        const jmtx_matrix_crs* mtx, const jmtx_matrix_crs* l, const jmtx_matrix_crs* u, const float* y, float* x,
        float* aux_vec, jmtx_solver_arguments* args);


jmtx_result jmtx_lu_solve_iterative_bmr(const jmtx_matrix_brm* a, const jmtx_matrix_brm* l, const jmtx_matrix_brm* u,
                                        const float y[const restrict], float x[const restrict],
                                        float aux_vec[const restrict], jmtx_solver_arguments* args);

jmtx_result jmtx_lu_solve_iterative_bmr_parallel(const jmtx_matrix_brm* a, const jmtx_matrix_brm* l,
                                                 const jmtx_matrix_brm* u,  const float y[const restrict],
                                                 float x[const restrict], float aux_vec[const restrict],
                                                 jmtx_solver_arguments* args);

#endif //JMTX_LU_SOLVING_H
