// Automatically generated from include/jmtx/float/solvers/bicgstab_iteration.h on Sun Dec 17 15:56:06 2023
//
// Created by jan on 17.6.2022.
//
// Properly done and understood on the 17.12.2023

#ifndef JMTXZ_BICGSTAB_ITERATION_H
#define JMTXZ_BICGSTAB_ITERATION_H
#ifndef JMTX_SOLVER_BASE_H
    #include "../../solver_base.h"
#endif

#ifdef JMTXZ_SPARSE_ROW_COMPRESSED_H
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
jmtx_result jmtxz_solve_iterative_bicgstab_crs(
        const jmtxz_matrix_crs* mtx, const _Complex double* restrict y, _Complex double* restrict x, _Complex double* restrict aux_vec1,
        _Complex double* restrict aux_vec2, _Complex double* restrict aux_vec3, _Complex double* restrict aux_vec4, _Complex double* restrict aux_vec5,
        _Complex double* restrict aux_vec6, jmtxd_solver_arguments* args);

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
jmtx_result jmtxzs_solve_iterative_bicgstab_crs(
        const jmtxz_matrix_crs* mtx, uint32_t n, const _Complex double y[JMTX_ARRAY_ATTRIB(restrict static n)], _Complex double x[JMTX_ARRAY_ATTRIB(restrict n)], _Complex double aux_vec1[JMTX_ARRAY_ATTRIB(restrict n)],
        _Complex double aux_vec2[JMTX_ARRAY_ATTRIB(restrict n)], _Complex double aux_vec3[JMTX_ARRAY_ATTRIB(restrict n)], _Complex double aux_vec4[JMTX_ARRAY_ATTRIB(restrict n)], _Complex double aux_vec5[JMTX_ARRAY_ATTRIB(restrict n)],
        _Complex double aux_vec6[JMTX_ARRAY_ATTRIB(restrict n)], jmtxd_solver_arguments* args);
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
jmtx_result jmtxz_solve_iterative_pilubicgstab_crs(
        const jmtxz_matrix_crs* mtx, const jmtxz_matrix_crs* l, const jmtxz_matrix_crs* u, const _Complex double* restrict y,
        _Complex double* restrict x, _Complex double* restrict aux_vec1, _Complex double* restrict aux_vec2, _Complex double* restrict aux_vec3,
        _Complex double* restrict aux_vec4, _Complex double* restrict aux_vec5, _Complex double* restrict aux_vec6, _Complex double* restrict aux_vec7,
        _Complex double* restrict aux_vec8, jmtxd_solver_arguments* args);

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
jmtx_result jmtxz_solve_iterative_pilubicgstab_crs_parallel(
        const jmtxz_matrix_crs* mtx, const jmtxz_matrix_crs* l, const jmtxz_matrix_crs* u, const _Complex double* restrict y,
        _Complex double* restrict x, _Complex double* restrict aux_vec1, _Complex double* restrict aux_vec2, _Complex double* restrict aux_vec3,
        _Complex double* restrict aux_vec4, _Complex double* restrict aux_vec5, _Complex double* restrict aux_vec6, _Complex double* restrict aux_vec7,
        _Complex double* restrict aux_vec8, jmtxd_solver_arguments* args);
#endif


#ifdef JMTXZ_SPARSE_DIAGONAL_COMPRESSED_H
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
jmtx_result jmtxz_solve_iterative_bicgstab_cds(
        const jmtxz_matrix_cds* mtx, const _Complex double* restrict y, _Complex double* restrict x, _Complex double* restrict aux_vec1,
        _Complex double* restrict aux_vec2, _Complex double* restrict aux_vec3, _Complex double* restrict aux_vec4, _Complex double* restrict aux_vec5,
        _Complex double* restrict aux_vec6, jmtxd_solver_arguments* args);

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
jmtx_result jmtxzs_solve_iterative_bicgstab_cds(
        const jmtxz_matrix_cds* mtx, uint32_t n, const _Complex double y[JMTX_ARRAY_ATTRIB(restrict static n)], _Complex double x[JMTX_ARRAY_ATTRIB(restrict n)], _Complex double aux_vec1[JMTX_ARRAY_ATTRIB(restrict n)],
        _Complex double aux_vec2[JMTX_ARRAY_ATTRIB(restrict n)], _Complex double aux_vec3[JMTX_ARRAY_ATTRIB(restrict n)], _Complex double aux_vec4[JMTX_ARRAY_ATTRIB(restrict n)], _Complex double aux_vec5[JMTX_ARRAY_ATTRIB(restrict n)],
        _Complex double aux_vec6[JMTX_ARRAY_ATTRIB(restrict n)], jmtxd_solver_arguments* args);
#endif



#endif //JMTXZ_BICGSTAB_ITERATION_H
