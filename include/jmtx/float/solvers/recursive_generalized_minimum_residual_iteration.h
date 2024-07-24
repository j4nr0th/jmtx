//
// Created by jan on 1.1.2024.
//

#ifndef JMTX_RECURSIVE_GENERALIZED_MINIMUM_RESIDUAL_ITERATION_H
#define JMTX_RECURSIVE_GENERALIZED_MINIMUM_RESIDUAL_ITERATION_H
#ifndef JMTX_SOLVER_BASE_H
    #include "../../solver_base.h"
#endif


#ifdef JMTX_SPARSE_DIAGONAL_COMPRESSED_H
    #ifndef JMTX_BAND_ROW_MAJOR_H
        #include "../matrices/band_row_major.h"
    #endif
/**
 * Applies Generalized Minimum Residual Recursive method (known as GMRESR) to solve a linear system A x = y.
 * Builds up a set of m orthonormal basis for the Krylov subspace in order to find an optimal search direction, then
 * makes its contribution to residual orthogonal to previous search directions. This makes it in general much better
 * behaved than GMRES with a short restart interval. Previously used search directions are saved up to the specified
 * number and removed from the previously applied directions.
 *
 * @param mtx system matrix A
 * @param y the solution of the system A x = y
 * @param x the solution vector which contains the initial guess of the solution
 * @param m the GMRES restart interval
 * @param l the CGR truncation interval
 * @param r_mtx an m by m upper triangular matrix (lbw = 0, ubw = m - 1) that is to be used in solving the least squares
 * problem
 * @param aux_vec1 auxiliary memory for a vector of m elements
 * @param aux_vec2 auxiliary memory for a vector of m elements
 * @param aux_vec3 auxiliary memory for a vector of m elements
 * @param aux_vec4 auxiliary memory for a vector of m elements
 * @param aux_vec5 auxiliary memory for a vector of m elements
 * @param aux_vec6 auxiliary memory for a vector of the same size as x and y
 * @param aux_vecs1 auxiliary memory for m vectors of the same size as x and y (n by m)
 * @param aux_vecs2 auxiliary memory for l vectors of the same size as x and y (n by m)
 * @param aux_vecs3 auxiliary memory for l vectors of the same size as x and y (n by m)
 * @param args::in_convergence_criterion tolerance to determine if the solution is close enough
 * @param args::in_max_iterations number of iterations to stop at
 * @param args::out_last_error receives the value of the error criterion at the final iteration
 * @param args::out_last_iteration receives the number of the final iteration
 * @param args::opt_error_evolution (optional) pointer to an array of length max_iterations, that receives the error value of each
 * iteration
 * @return JMTX_RESULT_SUCCESS if solution converged, JMTX_RESULT_NOT_CONVERGED if solution did not converge in the
 * given number of iterations, other error codes for other errors
 */
jmtx_result jmtx_solve_iterative_gmresr_cds(const jmtx_matrix_cds* mtx, const float* restrict y, float* restrict x,
                                             uint32_t m, uint32_t l, jmtx_matrix_brm* r_mtx,
                                             float aux_vec1[JMTX_ARRAY_ATTRIB(restrict m)], float aux_vec2[JMTX_ARRAY_ATTRIB(restrict m)],
                                             float aux_vec3[JMTX_ARRAY_ATTRIB(restrict m)], float aux_vec4[JMTX_ARRAY_ATTRIB(restrict m)],
                                             float aux_vec5[JMTX_ARRAY_ATTRIB(restrict m)], float* restrict aux_vec6,
                                             float* restrict aux_vecs1, float* restrict aux_vecs2,
                                             float* restrict aux_vecs3, jmtx_solver_arguments* args);
/**
 * Applies Generalized Minimum Residual Recursive method (known as GMRESR) to solve a linear system A x = y.
 * Builds up a set of m orthonormal basis for the Krylov subspace in order to find an optimal search direction, then
 * makes its contribution to residual orthogonal to previous search directions. This makes it in general much better
 * behaved than GMRES with a short restart interval. Previously used search directions are saved up to the specified
 * number and removed from the previously applied directions.
 *
 * @param mtx system matrix A
 * @param n the size of the system
 * @param y the solution of the system A x = y
 * @param x the solution vector which contains the initial guess of the solution
 * @param m the GMRES restart interval
 * @param l the CGR truncation interval
 * @param r_mtx an m by m upper triangular matrix (lbw = 0, ubw = m - 1) that is to be used in solving the least squares
 * problem
 * @param aux_vec1 auxiliary memory for a vector of m elements
 * @param aux_vec2 auxiliary memory for a vector of m elements
 * @param aux_vec3 auxiliary memory for a vector of m elements
 * @param aux_vec4 auxiliary memory for a vector of m elements
 * @param aux_vec5 auxiliary memory for a vector of m elements
 * @param aux_vec6 auxiliary memory for a vector of the same size as x and y
 * @param aux_vecs1 auxiliary memory for m vectors of the same size as x and y (n by m)
 * @param aux_vecs2 auxiliary memory for l vectors of the same size as x and y (n by l)
 * @param aux_vecs3 auxiliary memory for l vectors of the same size as x and y (n by l)
 * @param args::in_convergence_criterion tolerance to determine if the solution is close enough
 * @param args::in_max_iterations number of iterations to stop at
 * @param args::out_last_error receives the value of the error criterion at the final iteration
 * @param args::out_last_iteration receives the number of the final iteration
 * @param args::opt_error_evolution (optional) pointer to an array of length max_iterations, that receives the error value of each
 * iteration
 * @return JMTX_RESULT_SUCCESS if solution converged, JMTX_RESULT_NOT_CONVERGED if solution did not converge in the
 * given number of iterations, other error codes for other errors
 */
jmtx_result jmtxs_solve_iterative_gmresr_cds(const jmtx_matrix_cds* mtx, uint32_t n, const float y[JMTX_ARRAY_ATTRIB(static restrict n)],
                                             float x[JMTX_ARRAY_ATTRIB(static restrict n)], uint32_t m, uint32_t l, jmtx_matrix_brm* r_mtx,
                                             float aux_vec1[JMTX_ARRAY_ATTRIB(restrict m)], float aux_vec2[JMTX_ARRAY_ATTRIB(restrict m)],
                                             float aux_vec3[JMTX_ARRAY_ATTRIB(restrict m)], float aux_vec4[JMTX_ARRAY_ATTRIB(restrict m)],
                                             float aux_vec5[JMTX_ARRAY_ATTRIB(restrict m)], float aux_vec6[JMTX_ARRAY_ATTRIB(restrict n)],
                                             float aux_vecs1[JMTX_ARRAY_ATTRIB(restrict m * n)], float aux_vecs2[JMTX_ARRAY_ATTRIB(restrict l * n)],
                                             float aux_vecs3[JMTX_ARRAY_ATTRIB(restrict l * n)], jmtx_solver_arguments* args);
#endif

#ifdef JMTX_SPARSE_ROW_COMPRESSED_H
    #ifndef JMTX_BAND_ROW_MAJOR_H
        #include "../matrices/band_row_major.h"
    #endif
/**
 * Applies Generalized Minimum Residual Recursive method (known as GMRESR) to solve a linear system A x = y.
 * Builds up a set of m orthonormal basis for the Krylov subspace in order to find an optimal search direction, then
 * makes its contribution to residual orthogonal to previous search directions. This makes it in general much better
 * behaved than GMRES with a short restart interval. Previously used search directions are saved up to the specified
 * number and removed from the previously applied directions.
 *
 * @param mtx system matrix A
 * @param y the solution of the system A x = y
 * @param x the solution vector which contains the initial guess of the solution
 * @param m the GMRES restart interval
 * @param l the CGR truncation interval
 * @param r_mtx an m by m upper triangular matrix (lbw = 0, ubw = m - 1) that is to be used in solving the least squares
 * problem
 * @param aux_vec1 auxiliary memory for a vector of m elements
 * @param aux_vec2 auxiliary memory for a vector of m elements
 * @param aux_vec3 auxiliary memory for a vector of m elements
 * @param aux_vec4 auxiliary memory for a vector of m elements
 * @param aux_vec5 auxiliary memory for a vector of m elements
 * @param aux_vec6 auxiliary memory for a vector of the same size as x and y
 * @param aux_vecs1 auxiliary memory for m vectors of the same size as x and y (n by m)
 * @param aux_vecs2 auxiliary memory for l vectors of the same size as x and y (n by m)
 * @param aux_vecs3 auxiliary memory for l vectors of the same size as x and y (n by m)
 * @param args::in_convergence_criterion tolerance to determine if the solution is close enough
 * @param args::in_max_iterations number of iterations to stop at
 * @param args::out_last_error receives the value of the error criterion at the final iteration
 * @param args::out_last_iteration receives the number of the final iteration
 * @param args::opt_error_evolution (optional) pointer to an array of length max_iterations, that receives the error value of each
 * iteration
 * @return JMTX_RESULT_SUCCESS if solution converged, JMTX_RESULT_NOT_CONVERGED if solution did not converge in the
 * given number of iterations, other error codes for other errors
 */
jmtx_result jmtx_solve_iterative_gmresr_crs(const jmtx_matrix_crs* mtx, const float* restrict y, float* restrict x,
                                             uint32_t m, uint32_t l, jmtx_matrix_brm* r_mtx,
                                             float aux_vec1[JMTX_ARRAY_ATTRIB(restrict m)], float aux_vec2[JMTX_ARRAY_ATTRIB(restrict m)],
                                             float aux_vec3[JMTX_ARRAY_ATTRIB(restrict m)], float aux_vec4[JMTX_ARRAY_ATTRIB(restrict m)],
                                             float aux_vec5[JMTX_ARRAY_ATTRIB(restrict m)], float* restrict aux_vec6,
                                             float* restrict aux_vecs1, float* restrict aux_vecs2,
                                             float* restrict aux_vecs3, jmtx_solver_arguments* args);
/**
 * Applies Generalized Minimum Residual Recursive method (known as GMRESR) to solve a linear system A x = y.
 * Builds up a set of m orthonormal basis for the Krylov subspace in order to find an optimal search direction, then
 * makes its contribution to residual orthogonal to previous search directions. This makes it in general much better
 * behaved than GMRES with a short restart interval. Previously used search directions are saved up to the specified
 * number and removed from the previously applied directions.
 *
 * @param mtx system matrix A
 * @param n the size of the system
 * @param y the solution of the system A x = y
 * @param x the solution vector which contains the initial guess of the solution
 * @param m the GMRES restart interval
 * @param l the CGR truncation interval
 * @param r_mtx an m by m upper triangular matrix (lbw = 0, ubw = m - 1) that is to be used in solving the least squares
 * problem
 * @param aux_vec1 auxiliary memory for a vector of m elements
 * @param aux_vec2 auxiliary memory for a vector of m elements
 * @param aux_vec3 auxiliary memory for a vector of m elements
 * @param aux_vec4 auxiliary memory for a vector of m elements
 * @param aux_vec5 auxiliary memory for a vector of m elements
 * @param aux_vec6 auxiliary memory for a vector of the same size as x and y
 * @param aux_vecs1 auxiliary memory for m vectors of the same size as x and y (n by m)
 * @param aux_vecs2 auxiliary memory for l vectors of the same size as x and y (n by l)
 * @param aux_vecs3 auxiliary memory for l vectors of the same size as x and y (n by l)
 * @param args::in_convergence_criterion tolerance to determine if the solution is close enough
 * @param args::in_max_iterations number of iterations to stop at
 * @param args::out_last_error receives the value of the error criterion at the final iteration
 * @param args::out_last_iteration receives the number of the final iteration
 * @param args::opt_error_evolution (optional) pointer to an array of length max_iterations, that receives the error value of each
 * iteration
 * @return JMTX_RESULT_SUCCESS if solution converged, JMTX_RESULT_NOT_CONVERGED if solution did not converge in the
 * given number of iterations, other error codes for other errors
 */
jmtx_result jmtxs_solve_iterative_gmresr_crs(const jmtx_matrix_crs* mtx, uint32_t n, const float y[JMTX_ARRAY_ATTRIB(static restrict n)],
                                             float x[JMTX_ARRAY_ATTRIB(static restrict n)], uint32_t m, uint32_t l, jmtx_matrix_brm* r_mtx,
                                             float aux_vec1[JMTX_ARRAY_ATTRIB(restrict m)], float aux_vec2[JMTX_ARRAY_ATTRIB(restrict m)],
                                             float aux_vec3[JMTX_ARRAY_ATTRIB(restrict m)], float aux_vec4[JMTX_ARRAY_ATTRIB(restrict m)],
                                             float aux_vec5[JMTX_ARRAY_ATTRIB(restrict m)], float aux_vec6[JMTX_ARRAY_ATTRIB(restrict n)],
                                             float aux_vecs1[JMTX_ARRAY_ATTRIB(restrict m * n)], float aux_vecs2[JMTX_ARRAY_ATTRIB(restrict l * n)],
                                             float aux_vecs3[JMTX_ARRAY_ATTRIB(restrict l * n)], jmtx_solver_arguments* args);
#endif

#endif //JMTX_RECURSIVE_GENERALIZED_MINIMUM_RESIDUAL_ITERATION_H
