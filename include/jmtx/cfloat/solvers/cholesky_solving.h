// Automatically generated from include/jmtx/float/solvers/cholesky_solving.h on Fri Dec  1 17:35:57 2023
//
// Created by jan on 6.11.2023.
//

#ifndef JMTXC_CHOLESKY_SOLVING_H
#define JMTXC_CHOLESKY_SOLVING_H

#ifndef JMTXC_SPARSE_ROW_COMPRESSED_H
    #include "../matrices/sparse_row_compressed.h"
#endif
#ifndef JMTX_SOLVER_BASE_H
    #include "../../solver_base.h"
#endif

void jmtxc_cholesky_solve(const jmtxc_matrix_crs* c, const jmtxc_matrix_crs* ct, const _Complex float* restrict y, _Complex float* restrict x);

void jmtxc_cholesky_solve_inplace(const jmtxc_matrix_crs* c, const jmtxc_matrix_crs* ct, _Complex float* restrict x);

//jmtx_result jmtxc_incomplete_cholesky_decomposition_solve(
//        const jmtxc_matrix_crs* mtx, const _Complex float* y, _Complex float* x, _Complex float* aux_vec, jmtx_solver_arguments* args,
//        const jmtx_allocator_callbacks* allocator_callbacks);
//
//jmtx_result jmtxc_incomplete_cholesky_decomposition_solve_precomputed(
//        const jmtxc_matrix_crs* mtx, const jmtxc_matrix_crs* c, const jmtxc_matrix_crs* ct, const _Complex float* y, _Complex float* x,
//        _Complex float* aux_vec, jmtx_solver_arguments* args);
//
//jmtx_result jmtxc_incomplete_cholesky_decomposition_solve_parallel(
//        const jmtxc_matrix_crs* mtx, const _Complex float* y, _Complex float* x, _Complex float* aux_vec, jmtx_solver_arguments* args,
//        const jmtx_allocator_callbacks* allocator_callbacks);
//
//jmtx_result jmtxc_incomplete_cholesky_decomposition_solve_precomputed_parallel(
//        const jmtxc_matrix_crs* mtx, const jmtxc_matrix_crs* l, const jmtxc_matrix_crs* u, const _Complex float* y, _Complex float* x,
//        _Complex float* aux_vec, jmtx_solver_arguments* args);

#endif //JMTXC_CHOLESKY_SOLVING_H
