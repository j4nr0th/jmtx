//
// Created by jan on 17.6.2022.
//
// Properly done and understood on the 17.12.2023

#ifndef JMTX_BICGSTAB_ITERATION_H
#define JMTX_BICGSTAB_ITERATION_H
#ifndef JMTX_SOLVER_BASE_H
    #include "../../solver_base.h"
#endif

#ifdef JMTX_SPARSE_ROW_COMPRESSED_H
/**
 *  Solves the linear problem A x = y for a general matrix A by using the relations used for Bi-CG, but does not
 *  explicitly solve the adjoint problem, instead computing values by computing results of polynomial relations for it.
 *  Stabilized method also computes these indirectly by using a polynomial with a lower condition number, giving better
 *  convergence behaviour.
 *
 *  This version of the funciton does not check if its inputs are valid and just assumes they are.
 *
 * @param mtx system matrix A
 * @param y solution to the system A x = y
 * @param x the initial guess of x, which receives the final solution
 * @param aux_vec1 auxiliary memory used by the algorithm which needs to be the same size as x and y
 * @param aux_vec2 auxiliary memory used by the algorithm which needs to be the same size as x and y
 * @param aux_vec3 auxiliary memory used by the algorithm which needs to be the same size as x and y
 * @param aux_vec4 auxiliary memory used by the algorithm which needs to be the same size as x and y
 * @param aux_vec5 auxiliary memory used by the algorithm which needs to be the same size as x and y
 * @param aux_vec6 auxiliary memory used by the algorithm which needs to be the same size as x and y
 * @param args::in_convergence_criterion tolerance to determine if the solution is close enough
 * @param args::in_max_iterations number of iterations to stop at
 * @param args::out_last_error receives the value of the error criterion at the final iteration
 * @param args::out_last_iteration receives the number of the final iteration
 * @param args::opt_error_evolution (optional) pointer to an array of length max_iterations, that receives the error
 * value of each iteration
 * @return JMTX_RESULT_SUCCESS if solution converged, JMTX_RESULT_NOT_CONVERGED if solution did not converge in the
 * given number of iterations
 */
jmtx_result jmtx_solve_iterative_bicgstab_crs(
        const jmtx_matrix_crs* mtx, const float* restrict y, float* restrict x, float* restrict aux_vec1,
        float* restrict aux_vec2, float* restrict aux_vec3, float* restrict aux_vec4, float* restrict aux_vec5,
        float* restrict aux_vec6, jmtx_solver_arguments* args);

/**
 *  Solves the linear problem A x = y for a general matrix A by using the relations used for Bi-CG, but does not
 *  explicitly solve the adjoint problem, instead computing values by computing results of polynomial relations for it.
 *  Stabilized method also computes these indirectly by using a polynomial with a lower condition number, giving better
 *  convergence behaviour.
 *
 *  This version of the funciton checks for appropriate matrix type and dimensions, as well as for memory not
 *  overlapping.
 *
 * @param mtx system matrix A
 * @param y solution to the system A x = y
 * @param x the initial guess of x, which receives the final solution
 * @param aux_vec1 auxiliary memory used by the algorithm which needs to be the same size as x and y
 * @param aux_vec2 auxiliary memory used by the algorithm which needs to be the same size as x and y
 * @param aux_vec3 auxiliary memory used by the algorithm which needs to be the same size as x and y
 * @param aux_vec4 auxiliary memory used by the algorithm which needs to be the same size as x and y
 * @param aux_vec5 auxiliary memory used by the algorithm which needs to be the same size as x and y
 * @param aux_vec6 auxiliary memory used by the algorithm which needs to be the same size as x and y
 * @param args::in_convergence_criterion tolerance to determine if the solution is close enough
 * @param args::in_max_iterations number of iterations to stop at
 * @param args::out_last_error receives the value of the error criterion at the final iteration
 * @param args::out_last_iteration receives the number of the final iteration
 * @param args::opt_error_evolution (optional) pointer to an array of length max_iterations, that receives the error
 * value of each iteration
 * @return JMTX_RESULT_SUCCESS if solution converged, JMTX_RESULT_NOT_CONVERGED if solution did not converge in the
 * given number of iterations, other error codes in case of other errors
 */
jmtx_result jmtxs_solve_iterative_bicgstab_crs(
        const jmtx_matrix_crs* mtx, uint32_t n, const float y[restrict static n], float x[restrict n], float aux_vec1[restrict n],
        float aux_vec2[restrict n], float aux_vec3[restrict n], float aux_vec4[restrict n], float aux_vec5[restrict n],
        float aux_vec6[restrict n], jmtx_solver_arguments* args);
/**
 *  Solves the linear problem A x = y for a general matrix A by using the relations used for Bi-CG, but does not
 *  explicitly solve the adjoint problem, instead computing values by computing results of polynomial relations for it.
 *  Stabilized method also computes these indirectly by using a polynomial with a lower condition number, giving better
 *  convergence behaviour.
 *
 *  This version uses incomplete LU decomposition (ILU) of the matrix, which then allows for better convergence
 *  properties. The decomposition must be given to the function.
 *
 *  This version of the funciton does not check if its inputs are valid and just assumes they are.
 *
 * @param mtx system matrix A
 * @param l lower triangular matrix
 * @param u upper triangular matrix
 * @param y solution to the system A x = y
 * @param x the initial guess of x, which receives the final solution
 * @param aux_vec1 auxiliary memory used by the algorithm which needs to be the same size as x and y
 * @param aux_vec2 auxiliary memory used by the algorithm which needs to be the same size as x and y
 * @param aux_vec3 auxiliary memory used by the algorithm which needs to be the same size as x and y
 * @param aux_vec4 auxiliary memory used by the algorithm which needs to be the same size as x and y
 * @param aux_vec5 auxiliary memory used by the algorithm which needs to be the same size as x and y
 * @param aux_vec6 auxiliary memory used by the algorithm which needs to be the same size as x and y
 * @param aux_vec7 auxiliary memory used by the algorithm which needs to be the same size as x and y
 * @param aux_vec8 auxiliary memory used by the algorithm which needs to be the same size as x and y
 * @param args::in_convergence_criterion tolerance to determine if the solution is close enough
 * @param args::in_max_iterations number of iterations to stop at
 * @param args::out_last_error receives the value of the error criterion at the final iteration
 * @param args::out_last_iteration receives the number of the final iteration
 * @param args::opt_error_evolution (optional) pointer to an array of length max_iterations, that receives the error
 * value of each iteration
 * @return JMTX_RESULT_SUCCESS if solution converged, JMTX_RESULT_NOT_CONVERGED if solution did not converge in the
 * given number of iterations
 */
jmtx_result jmtx_solve_iterative_pilubicgstab_crs(
        const jmtx_matrix_crs* mtx, const jmtx_matrix_crs* l, const jmtx_matrix_crs* u, const float* restrict y,
        float* restrict x, float* restrict aux_vec1, float* restrict aux_vec2, float* restrict aux_vec3,
        float* restrict aux_vec4, float* restrict aux_vec5, float* restrict aux_vec6, float* restrict aux_vec7,
        float* restrict aux_vec8, jmtx_solver_arguments* args);

/**
 *  Solves the linear problem A x = y for a general matrix A by using the relations used for Bi-CG, but does not
 *  explicitly solve the adjoint problem, instead computing values by computing results of polynomial relations for it.
 *  Stabilized method also computes these indirectly by using a polynomial with a lower condition number, giving better
 *  convergence behaviour.
 *
 *  This version uses incomplete LU decomposition (ILU) of the matrix, which then allows for better convergence
 *  properties. The decomposition must be given to the function.
 *
 *  This version of the funciton does not check if its inputs are valid and just assumes they are.
 *
 *  This version uses OpenMP to solve the problem in parallel using multiple threads.
 *
 * @param mtx system matrix A
 * @param l lower triangular matrix
 * @param u upper triangular matrix
 * @param y solution to the system A x = y
 * @param x the initial guess of x, which receives the final solution
 * @param aux_vec1 auxiliary memory used by the algorithm which needs to be the same size as x and y
 * @param aux_vec2 auxiliary memory used by the algorithm which needs to be the same size as x and y
 * @param aux_vec3 auxiliary memory used by the algorithm which needs to be the same size as x and y
 * @param aux_vec4 auxiliary memory used by the algorithm which needs to be the same size as x and y
 * @param aux_vec5 auxiliary memory used by the algorithm which needs to be the same size as x and y
 * @param aux_vec6 auxiliary memory used by the algorithm which needs to be the same size as x and y
 * @param aux_vec7 auxiliary memory used by the algorithm which needs to be the same size as x and y
 * @param aux_vec8 auxiliary memory used by the algorithm which needs to be the same size as x and y
 * @param args::in_convergence_criterion tolerance to determine if the solution is close enough
 * @param args::in_max_iterations number of iterations to stop at
 * @param args::out_last_error receives the value of the error criterion at the final iteration
 * @param args::out_last_iteration receives the number of the final iteration
 * @param args::opt_error_evolution (optional) pointer to an array of length max_iterations, that receives the error
 * value of each iteration
 * @return JMTX_RESULT_SUCCESS if solution converged, JMTX_RESULT_NOT_CONVERGED if solution did not converge in the
 * given number of iterations
 */
jmtx_result jmtx_solve_iterative_pilubicgstab_crs_parallel(
        const jmtx_matrix_crs* mtx, const jmtx_matrix_crs* l, const jmtx_matrix_crs* u, const float* restrict y,
        float* restrict x, double* restrict aux_vec1, float* restrict aux_vec2, float* restrict aux_vec3,
        float* restrict aux_vec4, float* restrict aux_vec5, float* restrict aux_vec6, float* restrict aux_vec7,
        float* restrict aux_vec8, jmtx_solver_arguments* args);
#endif


#ifdef JMTX_SPARSE_DIAGONAL_COMPRESSED_H
/**
 *  Solves the linear problem A x = y for a general matrix A by using the relations used for Bi-CG, but does not
 *  explicitly solve the adjoint problem, instead computing values by computing results of polynomial relations for it.
 *  Stabilized method also computes these indirectly by using a polynomial with a lower condition number, giving better
 *  convergence behaviour.
 *
 *  This version of the funciton does not check if its inputs are valid and just assumes they are.
 *
 * @param mtx system matrix A
 * @param y solution to the system A x = y
 * @param x the initial guess of x, which receives the final solution
 * @param aux_vec1 auxiliary memory used by the algorithm which needs to be the same size as x and y
 * @param aux_vec2 auxiliary memory used by the algorithm which needs to be the same size as x and y
 * @param aux_vec3 auxiliary memory used by the algorithm which needs to be the same size as x and y
 * @param aux_vec4 auxiliary memory used by the algorithm which needs to be the same size as x and y
 * @param aux_vec5 auxiliary memory used by the algorithm which needs to be the same size as x and y
 * @param aux_vec6 auxiliary memory used by the algorithm which needs to be the same size as x and y
 * @param args::in_convergence_criterion tolerance to determine if the solution is close enough
 * @param args::in_max_iterations number of iterations to stop at
 * @param args::out_last_error receives the value of the error criterion at the final iteration
 * @param args::out_last_iteration receives the number of the final iteration
 * @param args::opt_error_evolution (optional) pointer to an array of length max_iterations, that receives the error
 * value of each iteration
 * @return JMTX_RESULT_SUCCESS if solution converged, JMTX_RESULT_NOT_CONVERGED if solution did not converge in the
 * given number of iterations
 */
jmtx_result jmtx_solve_iterative_bicgstab_cds(
        const jmtx_matrix_cds* mtx, const float* restrict y, float* restrict x, float* restrict aux_vec1,
        float* restrict aux_vec2, float* restrict aux_vec3, float* restrict aux_vec4, float* restrict aux_vec5,
        float* restrict aux_vec6, jmtx_solver_arguments* args);

/**
 *  Solves the linear problem A x = y for a general matrix A by using the relations used for Bi-CG, but does not
 *  explicitly solve the adjoint problem, instead computing values by computing results of polynomial relations for it.
 *  Stabilized method also computes these indirectly by using a polynomial with a lower condition number, giving better
 *  convergence behaviour.
 *
 *  This version of the funciton checks for appropriate matrix type and dimensions, as well as for memory not
 *  overlapping.
 *
 * @param mtx system matrix A
 * @param y solution to the system A x = y
 * @param x the initial guess of x, which receives the final solution
 * @param aux_vec1 auxiliary memory used by the algorithm which needs to be the same size as x and y
 * @param aux_vec2 auxiliary memory used by the algorithm which needs to be the same size as x and y
 * @param aux_vec3 auxiliary memory used by the algorithm which needs to be the same size as x and y
 * @param aux_vec4 auxiliary memory used by the algorithm which needs to be the same size as x and y
 * @param aux_vec5 auxiliary memory used by the algorithm which needs to be the same size as x and y
 * @param aux_vec6 auxiliary memory used by the algorithm which needs to be the same size as x and y
 * @param args::in_convergence_criterion tolerance to determine if the solution is close enough
 * @param args::in_max_iterations number of iterations to stop at
 * @param args::out_last_error receives the value of the error criterion at the final iteration
 * @param args::out_last_iteration receives the number of the final iteration
 * @param args::opt_error_evolution (optional) pointer to an array of length max_iterations, that receives the error
 * value of each iteration
 * @return JMTX_RESULT_SUCCESS if solution converged, JMTX_RESULT_NOT_CONVERGED if solution did not converge in the
 * given number of iterations, other error codes in case of other errors
 */
jmtx_result jmtxs_solve_iterative_bicgstab_cds(
        const jmtx_matrix_cds* mtx, uint32_t n, const float y[restrict static n], float x[restrict n], float aux_vec1[restrict n],
        float aux_vec2[restrict n], float aux_vec3[restrict n], float aux_vec4[restrict n], float aux_vec5[restrict n],
        float aux_vec6[restrict n], jmtx_solver_arguments* args);
#endif



#endif //JMTX_BICGSTAB_ITERATION_H
