//
// Created by jan on 1.1.2024.
//

#ifndef JMTXZ_GENERALIZED_MINIMUM_RESIDUAL_H
#define JMTXZ_GENERALIZED_MINIMUM_RESIDUAL_H
#ifndef JMTX_SOLVER_BASE_H
    #include "../../solver_base.h"
#endif
#ifndef JMTXZ_BAND_ROW_MAJOR_H
    #include "../matrices/band_row_major.h"
#endif

#ifdef JMTXZ_SPARSE_ROW_COMPRESSED_H

/**
 * Applies Generalized Minimum Residual method with a restart interval of M (known as GMRES(M)). Builds up a set of m
 * orthonormal basis for the Krylov subspace, then solves a least squares problem to minimize the residual using these
 * basis to solve a problem A x = y.
 *
 *
 * @param mtx system matrix A
 * @param y the solution of the system A x = y
 * @param x the solution vector which contains the initial guess of the solution
 * @param m the GMRES restart interval
 * @param r an m by m upper triangular matrix (lbw = 0, ubw = m - 1) that is to be used in solving the least squares
 * problem
 * @param aux_vec1 auxiliary memory for a vector of m elements
 * @param aux_vec2 auxiliary memory for a vector of m elements
 * @param aux_vec3 auxiliary memory for a vector of m elements
 * @param aux_vec4 auxiliary memory for a vector of m elements
 * @param aux_vec5 auxiliary memory for a vector of m elements
 * @param aux_vecs auxiliary memory for m vectors of the same size as x and y (n by m)
 * @param args::in_convergence_criterion tolerance to determine if the solution is close enough
 * @param args::in_max_iterations number of iterations to stop at
 * @param args::out_last_error receives the value of the error criterion at the final iteration
 * @param args::out_last_iteration receives the number of the final iteration
 * @param args::opt_error_evolution (optional) pointer to an array of length max_iterations, that receives the error value of each
 * iteration
 * @return JMTX_RESULT_SUCCESS if solution converged, JMTX_RESULT_NOT_CONVERGED if solution did not converge in the
 * given number of iterations, other error codes for other errors
 */
jmtx_result jmtxz_solve_iterative_gmresm_crs(const jmtxz_matrix_crs* mtx, const _Complex double* restrict y, _Complex double* restrict x,
                                             uint32_t m, jmtxz_matrix_brm* r, _Complex double aux_vec1[restrict m],
                                             _Complex double aux_vec2[restrict m], _Complex double aux_vec3[restrict m],
                                             _Complex double aux_vec4[restrict m], _Complex double aux_vec5[restrict m],
                                             _Complex double* restrict aux_vecs, jmtxd_solver_arguments* args);
#endif


#ifdef JMTXZ_SPARSE_DIAGONAL_COMPRESSED_H

/**
 * Applies Generalized Minimum Residual method with a restart interval of M (known as GMRES(M)). Builds up a set of m
 * orthonormal basis for the Krylov subspace, then solves a least squares problem to minimize the residual using these
 * basis to solve a problem A x = y.
 *
 *
 * @param mtx system matrix A
 * @param y the solution of the system A x = y
 * @param x the solution vector which contains the initial guess of the solution
 * @param m the GMRES restart interval
 * @param r an m by m upper triangular matrix (lbw = 0, ubw = m - 1) that is to be used in solving the least squares
 * problem
 * @param aux_vec1 auxiliary memory for a vector of m elements
 * @param aux_vec2 auxiliary memory for a vector of m elements
 * @param aux_vec3 auxiliary memory for a vector of m elements
 * @param aux_vec4 auxiliary memory for a vector of m elements
 * @param aux_vec5 auxiliary memory for a vector of m elements
 * @param aux_vecs auxiliary memory for m vectors of the same size as x and y (n by m)
 * @param args::in_convergence_criterion tolerance to determine if the solution is close enough
 * @param args::in_max_iterations number of iterations to stop at
 * @param args::out_last_error receives the value of the error criterion at the final iteration
 * @param args::out_last_iteration receives the number of the final iteration
 * @param args::opt_error_evolution (optional) pointer to an array of length max_iterations, that receives the error value of each
 * iteration
 * @return JMTX_RESULT_SUCCESS if solution converged, JMTX_RESULT_NOT_CONVERGED if solution did not converge in the
 * given number of iterations, other error codes for other errors
 */
jmtx_result jmtxz_solve_iterative_gmresm_cds(const jmtxz_matrix_cds* mtx, const _Complex double* restrict y, _Complex double* restrict x,
                                             uint32_t m, jmtxz_matrix_brm* r, _Complex double aux_vec1[restrict m],
                                             _Complex double aux_vec2[restrict m], _Complex double aux_vec3[restrict m],
                                             _Complex double aux_vec4[restrict m], _Complex double aux_vec5[restrict m],
                                             _Complex double* restrict aux_vecs, jmtxd_solver_arguments* args);

/**
 * Applies Generalized Minimum Residual method with a restart interval of M (known as GMRES(M)). Builds up a set of m
 * orthonormal basis for the Krylov subspace, then solves a least squares problem to minimize the residual using these
 * basis to solve a problem A x = y. Using a preconditioner means that the Krylov subspace is not based on the system
 * matrix A, but instead on the preconditioned matrix.
 *
 * Uses Right Preconditioning with the Jacobi iteration, meaning it uses the it actually solves a different system:
 *                                      A D⁻¹ z = y, then solves z = D x
 * This is done in hopes of A D⁻¹ having a lower condition number than A. Right preconditioning retains the exact
 * same residual as the original system A x = y, which may or may not be desired.
 *
 *
 * @param mtx system matrix A
 * @param y the solution of the system A x = y
 * @param x the solution vector which contains the initial guess of the solution
 * @param m the GMRES restart interval
 * @param r an m by m upper triangular matrix (lbw = 0, ubw = m - 1) that is to be used in solving the least squares
 * problem
 * @param aux_vec1 auxiliary memory for a vector of m elements
 * @param aux_vec2 auxiliary memory for a vector of m elements
 * @param aux_vec3 auxiliary memory for a vector of m elements
 * @param aux_vec4 auxiliary memory for a vector of m elements
 * @param aux_vec5 auxiliary memory for a vector of m elements
 * @param aux_vec6 auxiliary memory for a vector of of the same size as x and y
 * @param aux_vec7 auxiliary memory for a vector of of the same size as x and y
 * @param aux_vecs auxiliary memory for m vectors of the same size as x and y (n by m)
 * @param args::in_convergence_criterion tolerance to determine if the solution is close enough
 * @param args::in_max_iterations number of iterations to stop at
 * @param args::out_last_error receives the value of the error criterion at the final iteration
 * @param args::out_last_iteration receives the number of the final iteration
 * @param args::opt_error_evolution (optional) pointer to an array of length max_iterations, that receives the error value of each
 * iteration
 * @return JMTX_RESULT_SUCCESS if solution converged, JMTX_RESULT_NOT_CONVERGED if solution did not converge in the
 * given number of iterations, other error codes for other errors
 */
jmtx_result jmtxz_solve_iterative_gmresm_rpc_jacobi_cds(const jmtxz_matrix_cds* mtx, const _Complex double* restrict y,
                                                        _Complex double* restrict x, uint32_t m, jmtxz_matrix_brm* r,
                                                        _Complex double aux_vec1[restrict m], _Complex double aux_vec2[restrict m],
                                                        _Complex double aux_vec3[restrict m], _Complex double aux_vec4[restrict m],
                                                        _Complex double aux_vec5[restrict m], _Complex double* restrict aux_vec6,
                                                        _Complex double* restrict aux_vec7, _Complex double* restrict aux_vecs,
                                                        jmtxd_solver_arguments* args);

/**
 * Applies Generalized Minimum Residual method with a restart interval of M (known as GMRES(M)). Builds up a set of m
 * orthonormal basis for the Krylov subspace, then solves a least squares problem to minimize the residual using these
 * basis to solve a problem A x = y. Using a preconditioner means that the Krylov subspace is not based on the system
 * matrix A, but instead on the preconditioned matrix.
 *
 * Uses Left Preconditioning with the Jacobi iteration, meaning it uses the it actually solves a different system:
 *                                            D⁻¹ A x = D⁻¹ y
 * This is done in hopes of D⁻¹ A having a lower condition number than A. Left preconditioning has a the consequence of
 * not actually using/minimizing the real residual but instead the residual of the preconditioned system. This may or
 * may not be desired.
 *
 *
 * @param mtx system matrix A
 * @param y the solution of the system A x = y
 * @param x the solution vector which contains the initial guess of the solution
 * @param m the GMRES restart interval
 * @param r an m by m upper triangular matrix (lbw = 0, ubw = m - 1) that is to be used in solving the least squares
 * problem
 * @param aux_vec1 auxiliary memory for a vector of m elements
 * @param aux_vec2 auxiliary memory for a vector of m elements
 * @param aux_vec3 auxiliary memory for a vector of m elements
 * @param aux_vec4 auxiliary memory for a vector of m elements
 * @param aux_vec5 auxiliary memory for a vector of m elements
 * @param aux_vec6 auxiliary memory for a vector of of the same size as x and y
 * @param aux_vecs auxiliary memory for m vectors of the same size as x and y (n by m)
 * @param args::in_convergence_criterion tolerance to determine if the solution is close enough
 * @param args::in_max_iterations number of iterations to stop at
 * @param args::out_last_error receives the value of the error criterion at the final iteration
 * @param args::out_last_iteration receives the number of the final iteration
 * @param args::opt_error_evolution (optional) pointer to an array of length max_iterations, that receives the error value of each
 * iteration
 * @return JMTX_RESULT_SUCCESS if solution converged, JMTX_RESULT_NOT_CONVERGED if solution did not converge in the
 * given number of iterations, other error codes for other errors
 */
jmtx_result jmtxz_solve_iterative_gmresm_lpc_jacobi_cds(const jmtxz_matrix_cds* mtx, const _Complex double* restrict y,
                                                        _Complex double* restrict x, uint32_t m, jmtxz_matrix_brm* r,
                                                        _Complex double aux_vec1[restrict m], _Complex double aux_vec2[restrict m],
                                                        _Complex double aux_vec3[restrict m], _Complex double aux_vec4[restrict m],
                                                        _Complex double aux_vec5[restrict m], _Complex double* restrict aux_vec6,
                                                        _Complex double* restrict aux_vecs,
                                                        jmtxd_solver_arguments* args);

#endif

#endif //JMTXZ_GENERALIZED_MINIMUM_RESIDUAL_H
