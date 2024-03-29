// Automatically generated from include/jmtx/float/solvers/conjugate_gradient_iteration.h on Sun Dec 17 20:59:03 2023
//
// Created by jan on 30.10.2023.
//

#ifndef JMTXZ_CONJUGATE_GRADIENT_ITERATION_H
#define JMTXZ_CONJUGATE_GRADIENT_ITERATION_H

#ifndef JMTX_SOLVER_BASE_H
    #include "../../solver_base.h"
#endif

#ifdef JMTXZ_SPARSE_ROW_COMPRESSED_H
/**
 * Iterative solution method used to solve systems of equations Ax = y where A is symmetric positive definite (SPD).
 * Directly solves and N x N system in N iterations, but converges closely in fewer. Convergence speed is better for
 * systems with a lower condition number, which can be achieved by preconditioning. Stopping criterion is determined as
 * ||y - Au|| / ||y|| < e, where u is the iterative solution and e is the error tolerance.
 * @param mtx system matrix, that must be SPD
 * @param y result vector
 * @param x vector which receives the iterative solution
 * @param aux_vec1 auxiliary memory for a vector of the same size as x and y
 * @param aux_vec2 auxiliary memory for a vector of the same size as x and y
 * @param aux_vec3 auxiliary memory for a vector of the same size as x and y
 * @param args::in_convergence_criterion tolerance to determine if the solution is close enough
 * @param args::in_max_iterations number of iterations to stop at
 * @param args::out_last_error receives the value of the error criterion at the final iteration
 * @param args::out_last_iteration receives the number of the final iteration
 * @param args::opt_error_evolution (optional) pointer to an array of length max_iterations, that receives the error value of each
 * iteration
 * @return JMTX_RESULT_SUCCESS if solution converged, JMTX_RESULT_NOT_CONVERGED if solution did not converge in the
 * given number of iterations, JMTX_RESULT_STAGNATED if stagnation was detected, other error codes for other errors
 */
jmtx_result jmtxz_solve_iterative_conjugate_gradient_crs(
        const jmtxz_matrix_crs* mtx, const _Complex double* restrict y, _Complex double* restrict x,
        _Complex double* restrict aux_vec1, _Complex double* restrict aux_vec2,
        _Complex double* restrict aux_vec3, jmtxd_solver_arguments* args);


/**
 * Parallel version of the function jmtxz_solve_iterative_conjugate_gradient_crs using OpenMP
 *
 * Iterative solution method used to solve systems of equations Ax = y where A is symmetric positive definite (SPD).
 * Directly solves and N x N system in N iterations, but converges closely in fewer. Convergence speed is better for
 * systems with a lower condition number, which can be achieved by preconditioning. Stopping criterion is determined as
 * ||y - Au|| / ||y|| < e, where u is the iterative solution and e is the value of args::in_convergence_criterion.
 * @param mtx system matrix, that must be SPD
 * @param y result vector
 * @param x vector which receives the iterative solution
 * @param stagnation error reduction is less than this, give up
 * @param recalculation_interval after this many calculations, residual will be computed explicitly, instead of implicit
 * updates performed otherwise (value of 0 means to only do it when necessary to confirm convergence)
 * @param aux_vec1 auxiliary memory for a vector of the same size as x and y
 * @param aux_vec2 auxiliary memory for a vector of the same size as x and y
 * @param aux_vec3 auxiliary memory for a vector of the same size as x and y
 * @param args::in_convergence_criterion tolerance to determine if the solution is close enough
 * @param args::in_max_iterations number of iterations to stop at
 * @param args::out_last_error receives the value of the error criterion at the final iteration
 * @param args::out_last_iteration receives the number of the final iteration
 * @param args::opt_error_evolution (optional) pointer to an array of length max_iterations, that receives the error value of each
 * iteration
 * @return JMTX_RESULT_SUCCESS if solution converged, JMTX_RESULT_NOT_CONVERGED if solution did not converge in the
 * given number of iterations, JMTX_RESULT_STAGNATED if stagnation was detected, other error codes for other errors
 */
jmtx_result jmtxz_solve_iterative_conjugate_gradient_crs_parallel(
        const jmtxz_matrix_crs* mtx, const _Complex double* restrict y, _Complex double* restrict x,
        _Complex double* restrict aux_vec1, _Complex double* restrict aux_vec2,
        _Complex double* restrict aux_vec3, jmtxd_solver_arguments* args);



/**
 * Iterative solution method used to solve systems of equations Ax = y where A is symmetric positive definite (SPD).
 * Directly solves and N x N system in N iterations, but converges closely in fewer. Convergence speed is better for
 * systems with a lower condition number, which can be achieved by preconditioning. Stopping criterion is determined as
 * ||y - Au|| / ||y|| < e, where u is the iterative solution and e is the error tolerance.
 *
 * For this version preconditioning by incomplete Cholesky decomposition is used. The decomposition is given by the
 * lower triangular matrix C, such that C C^T = M is close to A, but is still sparse. The system to solve then is
 * M^-1 A x = M^-1 y, where M^-1 is not directly computed, but M^-1 a = b can be easily computed. It improves the
 * convergence by reducing the condition number of the matrix A, since the spectral radius of M^-1 A is less.
 *
 * @param mtx system matrix, that must be SPD
 * @param cho Incomplete Cholesky decomposition lower triangular matrix C
 * @param cho_t Transpose of incomplete Cholesky decomposition matrix C^T
 * @param y result vector
 * @param x vector which receives the iterative solution
 * @param stagnation error reduction is less than this, give up
 * @param recalculation_interval after this many calculations, residual will be computed explicitly, instead of implicit
 * updates performed otherwise (value of 0 means to only do it when necessary to confirm convergence)
 * @param aux_vec1 auxiliary memory for a vector of the same size as x and y
 * @param aux_vec2 auxiliary memory for a vector of the same size as x and y
 * @param aux_vec3 auxiliary memory for a vector of the same size as x and y
 * @param aux_vec4 auxiliary memory for a vector of the same size as x and y
 * @param args::in_convergence_criterion tolerance to determine if the solution is close enough
 * @param args::in_max_iterations number of iterations to stop at
 * @param args::out_last_error receives the value of the error criterion at the final iteration
 * @param args::out_last_iteration receives the number of the final iteration
 * @param args::opt_error_evolution (optional) pointer to an array of length max_iterations, that receives the error value of each
 * iteration
 * @return JMTX_RESULT_SUCCESS if solution converged, JMTX_RESULT_NOT_CONVERGED if solution did not converge in the
 * given number of iterations, JMTX_RESULT_STAGNATED if stagnation was detected, other error codes for other errors
 */
jmtx_result jmtxz_incomplete_cholesky_preconditioned_solve_iterative_conjugate_gradient_crs(
        const jmtxz_matrix_crs* mtx, const jmtxz_matrix_crs* cho, const jmtxz_matrix_crs* cho_t,  const _Complex double* restrict y,
        _Complex double* restrict x, _Complex double* restrict aux_vec1, _Complex double* restrict aux_vec2,
        _Complex double* restrict aux_vec3, _Complex double* restrict aux_vec4, jmtxd_solver_arguments* args);
#endif // JMTXZ_SPARSE_ROW_COMPRESSED_H

#ifdef JMTXZ_SPARSE_DIAGONAL_COMPRESSED_H
/**
 * Iterative solution method used to solve systems of equations Ax = y where A is symmetric positive definite (SPD).
 * Directly solves and N x N system in N iterations, but converges closely in fewer. Convergence speed is better for
 * systems with a lower condition number, which can be achieved by preconditioning. Stopping criterion is determined as
 * ||y - Au|| / ||y|| < e, where u is the iterative solution and e is the error tolerance.
 * @param mtx system matrix, that must be SPD
 * @param y result vector
 * @param x vector which receives the iterative solution
 * @param stagnation error reduction is less than this, give up
 * @param recalculation_interval after this many calculations, residual will be computed explicitly, instead of implicit
 * updates performed otherwise (value of 0 means to only do it when necessary to confirm convergence)
 * @param aux_vec1 auxiliary memory for a vector of the same size as x and y
 * @param aux_vec2 auxiliary memory for a vector of the same size as x and y
 * @param aux_vec3 auxiliary memory for a vector of the same size as x and y
 * @param args::in_convergence_criterion tolerance to determine if the solution is close enough
 * @param args::in_max_iterations number of iterations to stop at
 * @param args::out_last_error receives the value of the error criterion at the final iteration
 * @param args::out_last_iteration receives the number of the final iteration
 * @param args::opt_error_evolution (optional) pointer to an array of length max_iterations, that receives the error value of each
 * iteration
 * @return JMTX_RESULT_SUCCESS if solution converged, JMTX_RESULT_NOT_CONVERGED if solution did not converge in the
 * given number of iterations, JMTX_RESULT_STAGNATED if stagnation was detected, other error codes for other errors
 */
jmtx_result jmtxz_solve_iterative_conjugate_gradient_cds(const jmtxz_matrix_cds* mtx, const _Complex double* restrict y, _Complex double* restrict x,
                                        _Complex double* restrict aux_vec1, _Complex double* restrict aux_vec2, _Complex double* restrict aux_vec3,
                                        jmtxd_solver_arguments* args);
#endif //JMTXZ_SPARSE_DIAGONAL_COMPRESSED_H

#endif //JMTXZ_CONJUGATE_GRADIENT_ITERATION_H
