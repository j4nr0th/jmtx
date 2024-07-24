// Automatically generated from include/jmtx/float/solvers/gauss_seidel_iteration.h on Sun Dec 17 16:45:42 2023
//
// Created by jan on 16.6.2022.
//

#ifndef JMTXZ_GAUSS_SEIDEL_ITERATION_H
#define JMTXZ_GAUSS_SEIDEL_ITERATION_H
#ifndef JMTX_SOLVER_BASE_H
    #include "../../solver_base.h"
#endif

#ifdef JMTXZ_SPARSE_ROW_COMPRESSED_H
/*
 * Gauss-Seidel is an iterative method for solving the system Ax = y. It works by splitting the matrix A into
 * matrices D (diagonal), L (lower triangular with zero diagonal), and U (upper triangular with zero diagonal), such
 * that A = D + L + U. The equation is then expressed as x_(n+1) = (D + L)^{-1} (y - U x_(n)). This means that all the
 * functions in this file require that the diagonals of the matrices are non-zero.
 */

/**
 * Uses Gauss-Seidel (https://en.wikipedia.org/wiki/Gauss%E2%80%93Seidel_method)
 * to solve the linear system Ax = y
 *
 * @param mtx pointer to the memory where matrix A is stored as a band row major matrix
 * @param y pointer to the memory where the vector y is stored
 * @param x pointer to the memory where the solution vector x will be stored
 * @param aux_vec auxiliary memory for a vector of the same size as x and y
 * @param args::in_convergence_criterion tolerance to determine if the solution is close enough
 * @param args::in_max_iterations number of iterations to stop at
 * @param args::out_last_error receives the value of the error criterion at the final iteration
 * @param args::out_last_iteration receives the number of the final iteration
 * @param args::opt_error_evolution (optional) pointer to an array of length max_iterations, that receives the error value of each
 * iteration
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_NOT_CONVERGED if it hasn't reached given stopping criterion,
 * in case of failure it returns the associated error code
 */
jmtx_result jmtxz_solve_iterative_gauss_seidel_crs(const jmtxz_matrix_crs* mtx, const _Complex double* restrict y, _Complex double* restrict x,
                                  _Complex double* restrict aux_vec1, jmtxd_solver_arguments* args);

/**
 * Uses Gauss-Seidel (https://en.wikipedia.org/wiki/Gauss%E2%80%93Seidel_method)
 * to solve the linear system Ax = y
 *
 * @param mtx pointer to the memory where matrix A is stored as a band row major matrix
 * @param y pointer to the memory where the vector y is stored
 * @param x pointer to the memory where the solution vector x will be stored
 * @param aux_vec auxiliary memory for a vector of the same size as x and y
 * @param args::in_convergence_criterion tolerance to determine if the solution is close enough
 * @param args::in_max_iterations number of iterations to stop at
 * @param args::out_last_error receives the value of the error criterion at the final iteration
 * @param args::out_last_iteration receives the number of the final iteration
 * @param args::opt_error_evolution (optional) pointer to an array of length max_iterations, that receives the error value of each
 * iteration
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_NOT_CONVERGED if it hasn't reached given stopping criterion,
 * in case of failure it returns the associated error code
 */
jmtx_result jmtxzs_solve_iterative_gauss_seidel_crs(const jmtxz_matrix_crs* mtx, uint32_t n, const _Complex double y[JMTX_ARRAY_ATTRIB(static restrict n)],
                                   _Complex double x[JMTX_ARRAY_ATTRIB(restrict n)], _Complex double aux_vec1[JMTX_ARRAY_ATTRIB(restrict n)], jmtxd_solver_arguments* args);

#endif

#endif //JMTXZ_GAUSS_SEIDEL_ITERATION_H
