// Automatically generated from include/jmtx/float/solvers/lu_solving.h on Fri Dec  1 06:43:05 2023
//
// Created by jan on 6.11.2023.
//

#ifndef JMTXD_LU_SOLVING_H
#define JMTXD_LU_SOLVING_H
#ifndef JMTXD_BAND_ROW_MAJOR_H
    #include "../matrices/band_row_major.h"
#endif
#ifndef JMTXD_SPARSE_ROW_COMPRESSED_H
    #include "../matrices/sparse_row_compressed.h"
#endif
#ifndef JMTXD_SOLVER_BASE_H
    #include "../../solver_base.h"
#endif

void jmtxd_lu_solve_crs(const jmtxd_matrix_crs* l, const jmtxd_matrix_crs* u, const double* restrict y, double* restrict x);

void jmtxd_lu_solve_inplace_crs(const jmtxd_matrix_crs* l, const jmtxd_matrix_crs* u, double* restrict x);

void jmtxd_lu_solve_brm(const jmtxd_matrix_brm* l, const jmtxd_matrix_brm* u, const double* restrict y, double* restrict x);

void jmtxd_lu_solve_inplace_brm(const jmtxd_matrix_brm* l, const jmtxd_matrix_brm* u, double* restrict x);


jmtx_result jmtxd_incomplete_lu_decomposition_solve_crs(
        const jmtxd_matrix_crs* mtx, const double* y, double* x, double* aux_vec, jmtxd_solver_arguments* args,
        const jmtx_allocator_callbacks* allocator_callbacks);

jmtx_result jmtxd_incomplete_lu_decomposition_solve_precomputed_crs(
        const jmtxd_matrix_crs* mtx, const jmtxd_matrix_crs* l, const jmtxd_matrix_crs* u, const double* y, double* x,
        double* aux_vec, jmtxd_solver_arguments* args);

jmtx_result jmtxd_incomplete_lu_decomposition_solve_crs_parallel(
        const jmtxd_matrix_crs* mtx, const double* y, double* x, double* aux_vec, jmtxd_solver_arguments* args,
        const jmtx_allocator_callbacks* allocator_callbacks);

jmtx_result jmtxd_incomplete_lu_decomposition_solve_precomputed_crs_parallel(
        const jmtxd_matrix_crs* mtx, const jmtxd_matrix_crs* l, const jmtxd_matrix_crs* u, const double* y, double* x,
        double* aux_vec, jmtxd_solver_arguments* args);


jmtx_result jmtxd_lu_solve_iterative_bmr(const jmtxd_matrix_brm* a, const jmtxd_matrix_brm* l, const jmtxd_matrix_brm* u,
                                        const double y[const restrict], double x[const restrict],
                                        double aux_vec[const restrict], jmtxd_solver_arguments* args);

jmtx_result jmtxd_lu_solve_iterative_bmr_parallel(const jmtxd_matrix_brm* a, const jmtxd_matrix_brm* l,
                                                 const jmtxd_matrix_brm* u,  const double y[const restrict],
                                                 double x[const restrict], double aux_vec[const restrict],
                                                 jmtxd_solver_arguments* args);

#endif //JMTXD_LU_SOLVING_H
