// Automatically generated from include/jmtx/float/solvers/lu_solving.h on Fri Dec  1 17:35:57 2023
//
// Created by jan on 6.11.2023.
//

#ifndef JMTXC_LU_SOLVING_H
#define JMTXC_LU_SOLVING_H
#ifndef JMTXC_BAND_ROW_MAJOR_H
    #include "../matrices/band_row_major.h"
#endif
#ifndef JMTXC_SPARSE_ROW_COMPRESSED_H
    #include "../matrices/sparse_row_compressed.h"
#endif
#ifndef JMTX_SOLVER_BASE_H
    #include "../../solver_base.h"
#endif

void jmtxc_lu_solve_crs(const jmtxc_matrix_crs* l, const jmtxc_matrix_crs* u, const _Complex float* restrict y, _Complex float* restrict x);

void jmtxc_lu_solve_inplace_crs(const jmtxc_matrix_crs* l, const jmtxc_matrix_crs* u, _Complex float* restrict x);

void jmtxc_lu_solve_brm(const jmtxc_matrix_brm* l, const jmtxc_matrix_brm* u, const _Complex float* restrict y, _Complex float* restrict x);

void jmtxc_lu_solve_inplace_brm(const jmtxc_matrix_brm* l, const jmtxc_matrix_brm* u, _Complex float* restrict x);


jmtx_result jmtxc_incomplete_lu_decomposition_solve_crs(
        const jmtxc_matrix_crs* mtx, const _Complex float* y, _Complex float* x, _Complex float* aux_vec, jmtx_solver_arguments* args,
        const jmtx_allocator_callbacks* allocator_callbacks);

jmtx_result jmtxc_incomplete_lu_decomposition_solve_precomputed_crs(
        const jmtxc_matrix_crs* mtx, const jmtxc_matrix_crs* l, const jmtxc_matrix_crs* u, const _Complex float* y, _Complex float* x,
        _Complex float* aux_vec, jmtx_solver_arguments* args);

jmtx_result jmtxc_incomplete_lu_decomposition_solve_crs_parallel(
        const jmtxc_matrix_crs* mtx, const _Complex float* y, _Complex float* x, _Complex float* aux_vec, jmtx_solver_arguments* args,
        const jmtx_allocator_callbacks* allocator_callbacks);

jmtx_result jmtxc_incomplete_lu_decomposition_solve_precomputed_crs_parallel(
        const jmtxc_matrix_crs* mtx, const jmtxc_matrix_crs* l, const jmtxc_matrix_crs* u, const _Complex float* y, _Complex float* x,
        _Complex float* aux_vec, jmtx_solver_arguments* args);


jmtx_result jmtxc_lu_solve_iterative_bmr(const jmtxc_matrix_brm* a, const jmtxc_matrix_brm* l, const jmtxc_matrix_brm* u,
                                        const _Complex float y[const restrict], _Complex float x[const restrict],
                                        _Complex float aux_vec[const restrict], jmtx_solver_arguments* args);

jmtx_result jmtxc_lu_solve_iterative_bmr_parallel(const jmtxc_matrix_brm* a, const jmtxc_matrix_brm* l,
                                                 const jmtxc_matrix_brm* u,  const _Complex float y[const restrict],
                                                 _Complex float x[const restrict], _Complex float aux_vec[const restrict],
                                                 jmtx_solver_arguments* args);

#endif //JMTXC_LU_SOLVING_H
