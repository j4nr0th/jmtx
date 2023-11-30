// Automatically generated from include/jmtx/float/solvers/cholesky_solving.h on Thu Nov 30 19:26:51 2023
//
// Created by jan on 6.11.2023.
//

#ifndef JMTX_CHOLESKY_SOLVING_H
#define JMTX_CHOLESKY_SOLVING_H

#ifndef JMTX_SPARSE_ROW_COMPRESSED_H
    #include "../matrices/sparse_row_compressed.h"
#endif
#ifndef JMTX_SOLVER_BASE_H
    #include "solver_base.h"
#endif

void jmtxd_cholesky_solve(const jmtxd_matrix_crs* c, const jmtxd_matrix_crs* ct, const double* restrict y, double* restrict x);

void jmtxd_cholesky_solve_inplace(const jmtxd_matrix_crs* c, const jmtxd_matrix_crs* ct, double* restrict x);

//jmtx_result jmtxd_incomplete_cholesky_decomposition_solve(
//        const jmtxd_matrix_crs* mtx, const double* y, double* x, double* aux_vec, jmtxd_solver_arguments* args,
//        const jmtx_allocator_callbacks* allocator_callbacks);
//
//jmtx_result jmtxd_incomplete_cholesky_decomposition_solve_precomputed(
//        const jmtxd_matrix_crs* mtx, const jmtxd_matrix_crs* c, const jmtxd_matrix_crs* ct, const double* y, double* x,
//        double* aux_vec, jmtxd_solver_arguments* args);
//
//jmtx_result jmtxd_incomplete_cholesky_decomposition_solve_parallel(
//        const jmtxd_matrix_crs* mtx, const double* y, double* x, double* aux_vec, jmtxd_solver_arguments* args,
//        const jmtx_allocator_callbacks* allocator_callbacks);
//
//jmtx_result jmtxd_incomplete_cholesky_decomposition_solve_precomputed_parallel(
//        const jmtxd_matrix_crs* mtx, const jmtxd_matrix_crs* l, const jmtxd_matrix_crs* u, const double* y, double* x,
//        double* aux_vec, jmtxd_solver_arguments* args);

#endif //JMTX_CHOLESKY_SOLVING_H
