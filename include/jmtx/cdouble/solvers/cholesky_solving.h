// Automatically generated from include/jmtx/cfloat/solvers/cholesky_solving.h on Fri Dec  1 18:48:13 2023
// Automatically generated from include/jmtx/cdouble/solvers/cholesky_solving.h on Fri Dec  1 17:35:57 2023
//
// Created by jan on 6.11.2023.
//

#ifndef JMTXZ_CHOLESKY_SOLVING_H
#define JMTXZ_CHOLESKY_SOLVING_H

#ifndef JMTXZ_SPARSE_ROW_COMPRESSED_H
    #include "../matrices/sparse_row_compressed.h"
#endif
#ifndef JMTX_SOLVER_BASE_H
    #include "../../solver_base.h"
#endif

void jmtxz_cholesky_solve(const jmtxz_matrix_crs* c, const jmtxz_matrix_crs* ct, const _Complex double* restrict y, _Complex double* restrict x);

void jmtxz_cholesky_solve_inplace(const jmtxz_matrix_crs* c, const jmtxz_matrix_crs* ct, _Complex double* restrict x);

//jmtx_result jmtxz_incomplete_cholesky_decomposition_solve(
//        const jmtxz_matrix_crs* mtx, const _Complex double* y, _Complex double* x, _Complex double* aux_vec, jmtxd_solver_arguments* args,
//        const jmtx_allocator_callbacks* allocator_callbacks);
//
//jmtx_result jmtxz_incomplete_cholesky_decomposition_solve_precomputed(
//        const jmtxz_matrix_crs* mtx, const jmtxz_matrix_crs* c, const jmtxz_matrix_crs* ct, const _Complex double* y, _Complex double* x,
//        _Complex double* aux_vec, jmtxd_solver_arguments* args);
//
//jmtx_result jmtxz_incomplete_cholesky_decomposition_solve_parallel(
//        const jmtxz_matrix_crs* mtx, const _Complex double* y, _Complex double* x, _Complex double* aux_vec, jmtxd_solver_arguments* args,
//        const jmtx_allocator_callbacks* allocator_callbacks);
//
//jmtx_result jmtxz_incomplete_cholesky_decomposition_solve_precomputed_parallel(
//        const jmtxz_matrix_crs* mtx, const jmtxz_matrix_crs* l, const jmtxz_matrix_crs* u, const _Complex double* y, _Complex double* x,
//        _Complex double* aux_vec, jmtxd_solver_arguments* args);

#endif //JMTXZ_CHOLESKY_SOLVING_H
