// Automatically generated from include/jmtx/cfloat/solvers/gauss_seidel_iteration.h on Fri Dec  1 18:48:13 2023
// Automatically generated from include/jmtx/cdouble/solvers/gauss_seidel_iteration.h on Fri Dec  1 17:35:57 2023
//
// Created by jan on 16.6.2022.
//

#ifndef JMTXZ_GAUSS_SEIDEL_ITERATION_H
#define JMTXZ_GAUSS_SEIDEL_ITERATION_H
#ifndef JMTXZ_SPARSE_ROW_COMPRESSED_H
    #include "../matrices/sparse_row_compressed.h"
#endif
#ifndef JMTXZ_SOLVER_BASE_H
    #include "../../solver_base.h"
#endif
/*
 * Gauss-Seidel is an iterative method for solving the system Ax = y. It works by splitting the matrix A into
 * matrices D (diagonal), L (lower triangular with zero diagonal), and U (upper triangular with zero diagonal), such
 * that A = D + L + U. The equation is then expressed as x_(n+1) = (D + L)^{-1} (y - U x_(n)). This means that all the
 * functions in this file require that the diagonals of the matrices are non-zero.
 * Note: None of these functions do any safety checks. They are built to be fast and straight forward, not to hold your
 * hand. This will likely be a pain in the ass, but as long as they are tested well enough, it should be fine.
 */

/**
 * Uses Gauss-Seidel (https://en.wikipedia.org/wiki/Gauss%E2%80%93Seidel_method)
 * to solve the linear system Ax = y
 * @param mtx pointer to the memory where matrix A is stored as a compressed row sparse matrix
 * @param y pointer to the memory where the vector y is stored
 * @param x pointer to the memory where the solution vector x will be stored
 * @param convergence_dif when the largest value of change per iteration for an element in x is less than this,
 * the solution is considered found
 * @param n_max_iter maximum number of iterations to perform before giving up
 * @param p_iter pointer which if not null receives the number of iterations that were performed
 * @param p_final_error pointer which receives the final error (may be null)
 * @param p_error pointer to array of error evolution (may be null)
 * @param allocator_callbacks pointer to allocation callbacks used internally by the function (may be null)
 * @return zero if successful
 */
jmtx_result jmtxz_gauss_seidel_crs(const jmtxz_matrix_crs* mtx, const _Complex double* restrict y, _Complex double* restrict x, jmtxd_solver_arguments* args);


/**
 * Uses Gauss-Seidel (https://en.wikipedia.org/wiki/Gauss%E2%80%93Seidel_method)
 * to solve the linear system Ax = y. Uses modified version when running in parallel, randomly reading and writing
 * updated values, which (in best case) give same behaviour as Gauss-Seidel, but in worst case lead to behaviour as
 * expected by Point Jacobi.
 * @param mtx pointer to the memory where matrix A is stored as a compressed row sparse matrix
 * @param y pointer to the memory where the vector y is stored
 * @param x pointer to the memory where the solution vector x will be stored
 * @param convergence_dif when the ratio ||Ax|| / ||y|| is less than this value, iterations stop
 * @param n_max_iter maximum number of iterations to perform before giving up
 * @param p_final_iterations pointer which if not null receives the number of iterations that were performed
 * @param p_final_error pointer which receives the final error (may be null)
 * @param p_error_evolution pointer to array of error evolution (may be null)
 * @param aux_vector memory the size of vectors x and y to use as additional memory
 * @return zero if successful
 */
jmtx_result jmtxz_gauss_seidel_crs_parallel(
        const jmtxz_matrix_crs* mtx, const _Complex double* restrict y, _Complex double* restrict x, _Complex double* restrict aux_vector,
        jmtxd_solver_arguments* args);

#endif //JMTXZ_GAUSS_SEIDEL_ITERATION_H
