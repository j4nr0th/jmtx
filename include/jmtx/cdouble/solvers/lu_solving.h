// Automatically generated from include/jmtx/cfloat/solvers/lu_solving.h on Fri Dec  1 18:48:13 2023
// Automatically generated from include/jmtx/cdouble/solvers/lu_solving.h on Fri Dec  1 17:35:57 2023
//
// Created by jan on 6.11.2023.
//

#ifndef JMTXZ_LU_SOLVING_H
#define JMTXZ_LU_SOLVING_H
#ifndef JMTXZ_BAND_ROW_MAJOR_H
    #include "../matrices/band_row_major.h"
#endif
#ifndef JMTXZ_SPARSE_ROW_COMPRESSED_H
    #include "../matrices/sparse_row_compressed.h"
#endif
#ifndef JMTX_SOLVER_BASE_H
    #include "../../solver_base.h"
#endif

void jmtxz_lu_solve_crs(const jmtxz_matrix_crs* l, const jmtxz_matrix_crs* u, const _Complex double* restrict y, _Complex double* restrict x);

void jmtxz_lu_solve_inplace_crs(const jmtxz_matrix_crs* l, const jmtxz_matrix_crs* u, _Complex double* restrict x);

void jmtxz_lu_solve_brm(const jmtxz_matrix_brm* l, const jmtxz_matrix_brm* u, const _Complex double* restrict y, _Complex double* restrict x);

void jmtxz_lu_solve_inplace_brm(const jmtxz_matrix_brm* l, const jmtxz_matrix_brm* u, _Complex double* restrict x);


jmtx_result jmtxz_incomplete_lu_decomposition_solve_crs(
        const jmtxz_matrix_crs* mtx, const _Complex double* y, _Complex double* x, _Complex double* aux_vec, jmtxd_solver_arguments* args,
        const jmtx_allocator_callbacks* allocator_callbacks);

jmtx_result jmtxz_incomplete_lu_decomposition_solve_precomputed_crs(
        const jmtxz_matrix_crs* mtx, const jmtxz_matrix_crs* l, const jmtxz_matrix_crs* u, const _Complex double* y, _Complex double* x,
        _Complex double* aux_vec, jmtxd_solver_arguments* args);

jmtx_result jmtxz_incomplete_lu_decomposition_solve_crs_parallel(
        const jmtxz_matrix_crs* mtx, const _Complex double* y, _Complex double* x, _Complex double* aux_vec, jmtxd_solver_arguments* args,
        const jmtx_allocator_callbacks* allocator_callbacks);

jmtx_result jmtxz_incomplete_lu_decomposition_solve_precomputed_crs_parallel(
        const jmtxz_matrix_crs* mtx, const jmtxz_matrix_crs* l, const jmtxz_matrix_crs* u, const _Complex double* y, _Complex double* x,
        _Complex double* aux_vec, jmtxd_solver_arguments* args);


jmtx_result jmtxz_lu_solve_iterative_bmr(const jmtxz_matrix_brm* a, const jmtxz_matrix_brm* l, const jmtxz_matrix_brm* u,
                                        const _Complex double y[const restrict], _Complex double x[const restrict],
                                        _Complex double aux_vec[const restrict], jmtxd_solver_arguments* args);

jmtx_result jmtxz_lu_solve_iterative_bmr_parallel(const jmtxz_matrix_brm* a, const jmtxz_matrix_brm* l,
                                                 const jmtxz_matrix_brm* u,  const _Complex double y[const restrict],
                                                 _Complex double x[const restrict], _Complex double aux_vec[const restrict],
                                                 jmtxd_solver_arguments* args);

#endif //JMTXZ_LU_SOLVING_H
