//
// Created by jan on 30.10.2023.
//

#ifndef JMTX_CONJUGATE_GRADIENT_ITERATION_H
#define JMTX_CONJUGATE_GRADIENT_ITERATION_H
#include "../matrices/sparse_row_compressed.h"
#include "solver_base.h"

/**
 * Iterative solution method used to solve systems of equations Ax = y where A is symmetric positive definite (SPD).
 * Directly solves and N x N system in N iterations, but converges closely in fewer. Convergence speed is better for
 * systems with a lower condition number, which can be achieved by preconditioning. Stopping criterion is determined as
 * ||y - Au|| / ||y|| < e, where u is the iterative solution and e is the error tolerance.
 * @param mtx system matrix, that must be SPD
 * @param y result vector
 * @param x vector which receives the iterative solution
 * @param tolerance tolerance to determine if the solution is close enough
 * @param stagnation error reduction is less than this, give up
 * @param recalculation_interval after this many calculations, residual will be computed explicitly, instead of implicit
 * updates performed otherwise (value of 0 means to only do it when necessary to confirm convergence)
 * @param max_iterations number of iterations to stop at
 * @param p_final_iteration (optional) pointer that receives the number of the final iteration
 * @param err_evolution (optional) pointer to an array of length max_iterations, that receives the error value of each
 * iteration
 * @param final_error (optional) pointer that receives the value of the error criterion at the final iteration
 * @param aux_vec1 auxiliary memory for a vector of the same size as x and y
 * @param aux_vec2 auxiliary memory for a vector of the same size as x and y
 * @param aux_vec3 auxiliary memory for a vector of the same size as x and y
 * @return JMTX_RESULT_SUCCESS if solution converged, JMTX_RESULT_NOT_CONVERGED if solution did not converge in the
 * given number of iterations, JMTX_RESULT_STAGNATED if stagnation was detected, other error codes for other errors
 */
jmtx_result jmtx_conjugate_gradient_crs(
        const jmtx_matrix_crs* mtx, const float* y, float* x, const float stagnation,
        const uint32_t recalculation_interval, float* restrict aux_vec1, float* restrict aux_vec2,
        float* restrict aux_vec3, jmtx_solver_arguments* args);


/**
 * Parallel version of the function jmtx_conjugate_gradient_crs using OpenMP
 *
 * Iterative solution method used to solve systems of equations Ax = y where A is symmetric positive definite (SPD).
 * Directly solves and N x N system in N iterations, but converges closely in fewer. Convergence speed is better for
 * systems with a lower condition number, which can be achieved by preconditioning. Stopping criterion is determined as
 * ||y - Au|| / ||y|| < e, where u is the iterative solution and e is the error tolerance.
 * @param mtx system matrix, that must be SPD
 * @param y result vector
 * @param x vector which receives the iterative solution
 * @param tolerance tolerance to determine if the solution is close enough
 * @param stagnation error reduction is less than this, give up
 * @param recalculation_interval after this many calculations, residual will be computed explicitly, instead of implicit
 * updates performed otherwise (value of 0 means to only do it when necessary to confirm convergence)
 * @param max_iterations number of iterations to stop at
 * @param p_final_iteration (optional) pointer that receives the number of the final iteration
 * @param err_evolution (optional) pointer to an array of length max_iterations, that receives the error value of each
 * iteration
 * @param final_error (optional) pointer that receives the value of the error criterion at the final iteration
 * @param aux_vec1 auxiliary memory for a vector of the same size as x and y
 * @param aux_vec2 auxiliary memory for a vector of the same size as x and y
 * @return JMTX_RESULT_SUCCESS if solution converged, JMTX_RESULT_NOT_CONVERGED if solution did not converge in the
 * given number of iterations, JMTX_RESULT_STAGNATED if stagnation was detected, other error codes for other errors
 */
jmtx_result jmtx_conjugate_gradient_crs_parallel(
        const jmtx_matrix_crs* mtx, const float* y, float* x, const float stagnation,
        const uint32_t recalculation_interval, float* restrict aux_vec1, float* restrict aux_vec2,
        float* restrict aux_vec3, jmtx_solver_arguments* args);

#endif //JMTX_CONJUGATE_GRADIENT_ITERATION_H
