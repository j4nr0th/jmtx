// Automatically generated from include/jmtx/float/solvers/jacobi_point_iteration.h on Thu Dec 14 17:47:50 2023
//
// Created by jan on 15.6.2022.
//

#ifndef JMTXZ_JACOBI_POINT_ITERATION_H
#define JMTXZ_JACOBI_POINT_ITERATION_H
#ifndef JMTXZ_SOLVER_BASE_H
    #include "../../solver_base.h"
#endif

/*
 * Jacobi Point Iteration is an iterative method for solving the system Ax = y. It works by splitting the matrix A into
 * matrices D (diagonal), L (lower triangular with zero diagonal), and U (upper triangular with zero diagonal), such
 * that A = D + L + U. The equation is then expressed as x_(n+1) = D^{-1} (y - (L + U)x_(n)). This means that all the
 * functions in this file require that the diagonals of the matrices are non-zero.
 */


#ifdef JMTXZ_SPARSE_ROW_COMPRESSED_H
/**
 * Uses Jacobi point iteration (also known as Jacobi method: https://en.wikipedia.org/wiki/Jacobi_method)
 * to solve the linear system Ax = y
 *
 * @param mtx pointer to the memory where matrix A is stored as a compressed row sparse matrix
 * @param y pointer to the memory where the vector y is stored
 * @param x pointer to the memory where the solution vector x will be written
 * @param aux_vec1 auxiliary memory for a vector of the same size as x and y
 * @param aux_vec2 auxiliary memory for a vector of the same size as x and y
 * @param args::in_convergence_criterion tolerance to determine if the solution is close enough
 * @param args::in_max_iterations number of iterations to stop at
 * @param args::out_last_error receives the value of the error criterion at the final iteration
 * @param args::out_last_iteration receives the number of the final iteration
 * @param args::opt_error_evolution (optional) pointer to an array of length max_iterations, that receives the error value of each
 * iteration
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_NOT_CONVERGED if it hasn't reached given stopping criterion,
 * in case of failure it returns the associated error code
 */
jmtx_result jmtxz_jacobi_crs(
        const jmtxz_matrix_crs* mtx, const _Complex double* restrict y, _Complex double* restrict x, _Complex double* restrict aux_vec1, _Complex double* restrict aux_vec2,
        jmtxd_solver_arguments* args);

/**
 * Uses Jacobi point iteration (also known as Jacobi method: https://en.wikipedia.org/wiki/Jacobi_method)
 * to solve the linear system Ax = y
 *
 * @param mtx pointer to the memory where matrix A is stored as a compressed row sparse matrix
 * @param n size of the matrix and input vectors
 * @param y pointer to the memory where the vector y is stored
 * @param x pointer to the memory where the solution vector x will be stored
 * @param aux_vec1 auxiliary memory for a vector of the same size as x and y
 * @param aux_vec2 auxiliary memory for a vector of the same size as x and y
 * @param args::in_convergence_criterion tolerance to determine if the solution is close enough
 * @param args::in_max_iterations number of iterations to stop at
 * @param args::out_last_error receives the value of the error criterion at the final iteration
 * @param args::out_last_iteration receives the number of the final iteration
 * @param args::opt_error_evolution (optional) pointer to an array of length max_iterations, that receives the error value of each
 * iteration
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_NOT_CONVERGED if it hasn't reached given stopping criterion,
 * in case of failure it returns the associated error code
 */
jmtx_result jmtxzs_jacobi_crs(
        const jmtxz_matrix_crs* mtx, uint32_t n, const _Complex double y[static restrict n], _Complex double x[restrict n], _Complex double aux_vec1[restrict n], _Complex double aux_vec2[restrict n],
        jmtxd_solver_arguments* args);

/**
 * Uses Jacobi point iteration (also known as Jacobi method: https://en.wikipedia.org/wiki/Jacobi_method)
 * to solve the linear system Ax = y. Uses a relaxation factor ω for the following relation:
 *
 *  x(n + 1) = ω D⁻¹ (y - A x(n)) + x(n)
 *
 *  This may allow for better convergence
 *
 * @param mtx pointer to the memory where matrix A is stored as a compressed row sparse matrix
 * @param y pointer to the memory where the vector y is stored
 * @param x pointer to the memory where the solution vector x will be stored
 * @param relaxation_factor relaxation factor used for the iterations, which must be greater than zero
 * @param aux_vec1 auxiliary memory for a vector of the same size as x and y
 * @param aux_vec2 auxiliary memory for a vector of the same size as x and y
 * @param args::in_convergence_criterion tolerance to determine if the solution is close enough
 * @param args::in_max_iterations number of iterations to stop at
 * @param args::out_last_error receives the value of the error criterion at the final iteration
 * @param args::out_last_iteration receives the number of the final iteration
 * @param args::opt_error_evolution (optional) pointer to an array of length max_iterations, that receives the error value of each
 * iteration
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_NOT_CONVERGED if it hasn't reached given stopping criterion,
 * in case of failure it returns the associated error code
 */
jmtx_result jmtxz_jacobi_relaxed_crs(
        const jmtxz_matrix_crs* mtx, const _Complex double* restrict y, _Complex double* restrict x, _Complex double relaxation_factor, _Complex double* restrict aux_vec1,
        _Complex double* restrict aux_vec2, jmtxd_solver_arguments* args);


/**
 * Uses Jacobi point iteration (also known as Jacobi method: https://en.wikipedia.org/wiki/Jacobi_method)
 * to solve the linear system Ax = y.
 *
 * This version of the function uses OpenMP to solve the system in parallel.
 *
 * @param mtx pointer to the memory where matrix A is stored as a compressed row sparse matrix
 * @param y pointer to the memory where the vector y is stored
 * @param x pointer to the memory where the solution vector x will be stored
 * @param aux_vec1 auxiliary memory for a vector of the same size as x and y
 * @param aux_vec2 auxiliary memory for a vector of the same size as x and y
 * @param args::in_convergence_criterion tolerance to determine if the solution is close enough
 * @param args::in_max_iterations number of iterations to stop at
 * @param args::out_last_error receives the value of the error criterion at the final iteration
 * @param args::out_last_iteration receives the number of the final iteration
 * @param args::opt_error_evolution (optional) pointer to an array of length max_iterations, that receives the error value of each
 * iteration
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_NOT_CONVERGED if it hasn't reached given stopping criterion,
 * in case of failure it returns the associated error code
 */
jmtx_result jmtxz_jacobi_crs_parallel(
        const jmtxz_matrix_crs* mtx, const _Complex double* restrict y, _Complex double* restrict x, _Complex double* restrict aux_vector1, _Complex double* restrict aux_vector2,
        jmtxd_solver_arguments* args);
#endif



#ifdef JMTXZ_SPARSE_DIAGONAL_COMPRESSED_H
/**
 * Uses Jacobi point iteration (also known as Jacobi method: https://en.wikipedia.org/wiki/Jacobi_method)
 * to solve the linear system Ax = y
 *
 * @param mtx pointer to the memory where matrix A is stored as a compressed diagonal sparse matrix
 * @param y pointer to the memory where the vector y is stored
 * @param x pointer to the memory where the solution vector x will be stored
 * @param aux_vec1 auxiliary memory for a vector of the same size as x and y
 * @param aux_vec2 auxiliary memory for a vector of the same size as x and y
 * @param args::in_convergence_criterion tolerance to determine if the solution is close enough
 * @param args::in_max_iterations number of iterations to stop at
 * @param args::out_last_error receives the value of the error criterion at the final iteration
 * @param args::out_last_iteration receives the number of the final iteration
 * @param args::opt_error_evolution (optional) pointer to an array of length max_iterations, that receives the error value of each
 * iteration
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_NOT_CONVERGED if it hasn't reached given stopping criterion,
 * in case of failure it returns the associated error code
 */
jmtx_result jmtxz_jacobi_cds(
        const jmtxz_matrix_cds* mtx, const _Complex double* restrict y, _Complex double* restrict x, _Complex double* restrict aux_vec1, _Complex double* restrict aux_vec2,
        jmtxd_solver_arguments* args);

/**
 * Uses Jacobi point iteration (also known as Jacobi method: https://en.wikipedia.org/wiki/Jacobi_method)
 * to solve the linear system Ax = y. Uses a relaxation factor ω for the following relation:
 *
 *  x(n + 1) = ω D⁻¹ (y - A x(n)) + x(n)
 *
 *  This may allow for better convergence
 *
 * @param mtx pointer to the memory where matrix A is stored as a compressed diagonal sparse matrix
 * @param y pointer to the memory where the vector y is stored
 * @param x pointer to the memory where the solution vector x will be stored
 * @param relaxation_factor relaxation factor used for the iterations, which must be greater than zero
 * @param aux_vec1 auxiliary memory for a vector of the same size as x and y
 * @param aux_vec2 auxiliary memory for a vector of the same size as x and y
 * @param args::in_convergence_criterion tolerance to determine if the solution is close enough
 * @param args::in_max_iterations number of iterations to stop at
 * @param args::out_last_error receives the value of the error criterion at the final iteration
 * @param args::out_last_iteration receives the number of the final iteration
 * @param args::opt_error_evolution (optional) pointer to an array of length max_iterations, that receives the error value of each
 * iteration
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_NOT_CONVERGED if it hasn't reached given stopping criterion,
 * in case of failure it returns the associated error code
 */
jmtx_result jmtxz_jacobi_relaxed_cds(
        const jmtxz_matrix_cds* mtx, const _Complex double* restrict y, _Complex double* restrict x, _Complex double relaxation_factor, _Complex double* restrict aux_vec1,
        _Complex double* restrict aux_vec2, jmtxd_solver_arguments* args);
#endif




#ifdef JMTXZ_BAND_ROW_MAJOR_H
/**
 * Uses Jacobi point iteration (also known as Jacobi method: https://en.wikipedia.org/wiki/Jacobi_method)
 * to solve the linear system Ax = y
 *
 * @param mtx pointer to the memory where matrix A is stored as a band row major matrix
 * @param y pointer to the memory where the vector y is stored
 * @param x pointer to the memory where the solution vector x will be stored
 * @param aux_vec1 auxiliary memory for a vector of the same size as x and y
 * @param aux_vec2 auxiliary memory for a vector of the same size as x and y
 * @param args::in_convergence_criterion tolerance to determine if the solution is close enough
 * @param args::in_max_iterations number of iterations to stop at
 * @param args::out_last_error receives the value of the error criterion at the final iteration
 * @param args::out_last_iteration receives the number of the final iteration
 * @param args::opt_error_evolution (optional) pointer to an array of length max_iterations, that receives the error value of each
 * iteration
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_NOT_CONVERGED if it hasn't reached given stopping criterion,
 * in case of failure it returns the associated error code
 */
jmtx_result jmtxz_jacobi_brm(
        const jmtxz_matrix_brm* mtx, const _Complex double* restrict y, _Complex double* restrict x, _Complex double* restrict aux_vec1, _Complex double* restrict aux_vec2,
        jmtxd_solver_arguments* args);

/**
 * Uses Jacobi point iteration (also known as Jacobi method: https://en.wikipedia.org/wiki/Jacobi_method)
 * to solve the linear system Ax = y. Uses a relaxation factor ω for the following relation:
 *
 *  x(n + 1) = ω D⁻¹ (y - A x(n)) + x(n)
 *
 *  This may allow for better convergence
 *
 * @param mtx pointer to the memory where matrix A is stored as a band row major matrix
 * @param y pointer to the memory where the vector y is stored
 * @param x pointer to the memory where the solution vector x will be stored
 * @param aux_vec1 auxiliary memory for a vector of the same size as x and y
 * @param aux_vec2 auxiliary memory for a vector of the same size as x and y
 * @param relaxation_factor relaxation factor used for the iterations, which must be greater than zero
 * @param args::in_convergence_criterion tolerance to determine if the solution is close enough
 * @param args::in_max_iterations number of iterations to stop at
 * @param args::out_last_error receives the value of the error criterion at the final iteration
 * @param args::out_last_iteration receives the number of the final iteration
 * @param args::opt_error_evolution (optional) pointer to an array of length max_iterations, that receives the error value of each
 * iteration
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_NOT_CONVERGED if it hasn't reached given stopping criterion,
 * in case of failure it returns the associated error code
 */
jmtx_result jmtxz_jacobi_relaxed_brm(
        const jmtxz_matrix_brm* mtx, const _Complex double* restrict y, _Complex double* restrict x, _Complex double relaxation_factor, _Complex double* restrict aux_vec1,
        _Complex double* restrict aux_vec2, jmtxd_solver_arguments* args);


/**
 * Uses Jacobi point iteration (also known as Jacobi method: https://en.wikipedia.org/wiki/Jacobi_method)
 * to solve the linear system Ax = y.
 *
 * This version of the function uses OpenMP to solve the system in parallel.
 *
 * @param mtx pointer to the memory where matrix A is stored as a band row major matrix
 * @param y pointer to the memory where the vector y is stored
 * @param x pointer to the memory where the solution vector x will be stored
 * @param aux_vec1 auxiliary memory for a vector of the same size as x and y
 * @param aux_vec2 auxiliary memory for a vector of the same size as x and y
 * @param args::in_convergence_criterion tolerance to determine if the solution is close enough
 * @param args::in_max_iterations number of iterations to stop at
 * @param args::out_last_error receives the value of the error criterion at the final iteration
 * @param args::out_last_iteration receives the number of the final iteration
 * @param args::opt_error_evolution (optional) pointer to an array of length max_iterations, that receives the error value of each
 * iteration
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_NOT_CONVERGED if it hasn't reached given stopping criterion,
 * in case of failure it returns the associated error code
 */
jmtx_result jmtxz_jacobi_brm_parallel(
        const jmtxz_matrix_brm* mtx, const _Complex double* restrict y, _Complex double* restrict x, _Complex double* restrict aux_vector1, _Complex double* restrict aux_vector2,
        jmtxd_solver_arguments* args);
#endif

#endif //JMTXZ_JACOBI_POINT_ITERATION_H
