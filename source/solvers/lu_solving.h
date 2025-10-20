#ifndef JMTX_LU_SOLVING_H
#define JMTX_LU_SOLVING_H
#include "../solver_base.h"

#ifdef JMTX_BAND_ROW_MAJOR_H

/**
 * Solves a problem L U x = y, where L is a lower triangular matrix with the diagonal equal to 1 and U is an upper
 * triangular matrix.
 * @param l lower triangular matrix with all entries on its main diagonal equal to 1
 * @param u upper triangular matrix
 * @param y memory containing forcing vector
 * @param x memory which receives the solution
 */
void JMTX_NAME_TYPED(solve_direct_lu_brm)(const JMTX_NAME_TYPED(matrix_brm) * l, const JMTX_NAME_TYPED(matrix_brm) * u,
                                          const JMTX_SCALAR_T *restrict y, JMTX_SCALAR_T *restrict x);

/**
 * Solves a problem L U x = y, where L is a lower triangular matrix with the diagonal equal to 1 and U is an upper
 * triangular matrix. This version of the function stores the solution vector x back into the same memory where the
 * forcing vector was.
 * @param l lower triangular matrix with all entries on its main diagonal equal to 1
 * @param u upper triangular matrix
 * @param x memory which contains the forcing vector and receives the solution
 */
void JMTX_NAME_TYPED(solve_direct_lu_brm_inplace)(const JMTX_NAME_TYPED(matrix_brm) * l,
                                                  const JMTX_NAME_TYPED(matrix_brm) * u, JMTX_SCALAR_T *restrict x);

/**
 * Solves the A x = L U x = y problem by computing the residual, then solving for L U e = r for the error e if residual
 * is too large to further refine the solution. This can help eliminate rounding errors, or alternatively can be used
 * with an incomplete LU decomposition (or ILU) to work as an iterative solver. Must be given the matrices LU.
 * @param a system matrix A
 * @param l lower triangular matrix with all entries on its main diagonal equal to 1
 * @param u upper triangular matrix
 * @param y memory containing forcing vector
 * @param x memory which receives the solution
 * @param aux_vec auxiliary memory for a vector of the same size as x and y
 * @param args::in_convergence_criterion tolerance to determine if the solution is close enough
 * @param args::in_max_iterations number of iterations to stop at
 * @param args::out_last_error receives the value of the error criterion at the final iteration
 * @param args::out_last_iteration receives the number of the final iteration
 * @param args::opt_error_evolution (optional) pointer to an array of length max_iterations, that receives the error
 * value of each iteration
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_NOT_CONVERGED if it hasn't reached given stopping criterion,
 * in case of failure it returns the associated error code
 */
jmtx_result JMTX_NAME_TYPED(solve_iterative_lu_brm_refine)(
    const JMTX_NAME_TYPED(matrix_brm) * a, const JMTX_NAME_TYPED(matrix_brm) * l, const JMTX_NAME_TYPED(matrix_brm) * u,
    const JMTX_SCALAR_T y[JMTX_ARRAY_ATTRIB(restrict)], JMTX_SCALAR_T x[JMTX_ARRAY_ATTRIB(restrict)],
    JMTX_SCALAR_T aux_vec[JMTX_ARRAY_ATTRIB(restrict)], JMTX_NAME_TYPED(solver_arguments) * args);

/**
 * Solves the A x = L U x = y problem by computing the residual, then solving for L U e = r for the error e if residual
 * is too large to further refine the solution. This can help eliminate rounding errors, or alternatively can be used
 * with an incomplete LU decomposition (or ILU) to work as an iterative solver. Must be given the matrices LU.
 * This version offers a minor degree of parallelization, where calculation of the residual and applying error
 * correction are done in parallel, but the main bottleneck of dealing with inverting L U e = r is done in series.
 * @param a system matrix A
 * @param l lower triangular matrix with all entries on its main diagonal equal to 1
 * @param u upper triangular matrix
 * @param y memory containing forcing vector
 * @param x memory which receives the solution
 * @param aux_vec auxiliary memory for a vector of the same size as x and y
 * @param args::in_convergence_criterion tolerance to determine if the solution is close enough
 * @param args::in_max_iterations number of iterations to stop at
 * @param args::out_last_error receives the value of the error criterion at the final iteration
 * @param args::out_last_iteration receives the number of the final iteration
 * @param args::opt_error_evolution (optional) pointer to an array of length max_iterations, that receives the error
 * value of each iteration
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_NOT_CONVERGED if it hasn't reached given stopping criterion,
 * in case of failure it returns the associated error code
 */
jmtx_result JMTX_NAME_TYPED(solve_iterative_lu_brm_refine_parallel)(
    const JMTX_NAME_TYPED(matrix_brm) * a, const JMTX_NAME_TYPED(matrix_brm) * l, const JMTX_NAME_TYPED(matrix_brm) * u,
    const JMTX_SCALAR_T y[JMTX_ARRAY_ATTRIB(const restrict)], JMTX_SCALAR_T x[JMTX_ARRAY_ATTRIB(const restrict)],
    JMTX_SCALAR_T aux_vec[JMTX_ARRAY_ATTRIB(const restrict)], JMTX_NAME_TYPED(solver_arguments) * args);
#endif
#ifdef JMTX_SPARSE_ROW_COMPRESSED_H

/**
 * Solves a problem L U x = y, where L is a lower triangular matrix with the diagonal equal to 1 and U is an upper
 * triangular matrix.
 * @param l lower triangular matrix with all entries on its main diagonal equal to 1
 * @param u upper triangular matrix
 * @param y memory containing forcing vector
 * @param x memory which receives the solution
 */
void JMTX_NAME_TYPED(solve_direct_lu_crs)(const JMTX_NAME_TYPED(matrix_crs) * l, const JMTX_NAME_TYPED(matrix_crs) * u,
                                          const JMTX_SCALAR_T *restrict y, JMTX_SCALAR_T *restrict x);

/**
 * Solves a problem L U x = y, where L is a lower triangular matrix with the diagonal equal to 1 and U is an upper
 * triangular matrix. This version of the function stores the solution vector x back into the same memory where the
 * forcing vector was.
 * @param l lower triangular matrix with all entries on its main diagonal equal to 1
 * @param u upper triangular matrix
 * @param x memory which contains the forcing vector and receives the solution
 */
void JMTX_NAME_TYPED(solve_direct_lu_crs_inplace)(const JMTX_NAME_TYPED(matrix_crs) * l,
                                                  const JMTX_NAME_TYPED(matrix_crs) * u, JMTX_SCALAR_T *restrict x);

/**
 * Solves the A x = L U x = y problem by computing the residual, then solving for L U e = r for the error e if residual
 * is too large to further refine the solution. This can help eliminate rounding errors, or alternatively can be used
 * with an incomplete LU decomposition (or ILU) to work as an iterative solver. Must be given the matrices LU.
 * @param a system matrix A
 * @param y memory containing forcing vector
 * @param x memory which receives the solution
 * @param aux_vec auxiliary memory for a vector of the same size as x and y
 * @param args::in_convergence_criterion tolerance to determine if the solution is close enough
 * @param args::in_max_iterations number of iterations to stop at
 * @param args::out_last_error receives the value of the error criterion at the final iteration
 * @param args::out_last_iteration receives the number of the final iteration
 * @param args::opt_error_evolution (optional) pointer to an array of length max_iterations, that receives the error
 * value of each iteration
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_NOT_CONVERGED if it hasn't reached given stopping criterion,
 * in case of failure it returns the associated error code
 */
jmtx_result JMTX_NAME_TYPED(solve_iterative_ilu_crs)(const JMTX_NAME_TYPED(matrix_crs) * mtx,
                                                     const JMTX_SCALAR_T *restrict y, JMTX_SCALAR_T *restrict x,
                                                     JMTX_SCALAR_T *aux_vec, JMTX_NAME_TYPED(solver_arguments) * args,
                                                     const jmtx_allocator_callbacks *allocator_callbacks);

/**
 * Solves the A x = L U x = y problem by computing the residual, then solving for L U e = r for the error e if residual
 * is too large to further refine the solution. This can help eliminate rounding errors, or alternatively can be used
 * with an incomplete LU decomposition (or ILU) to work as an iterative solver. Must be given the matrices LU.
 * @param a system matrix A
 * @param l lower triangular matrix with all entries on its main diagonal equal to 1
 * @param u upper triangular matrix
 * @param y memory containing forcing vector
 * @param x memory which receives the solution
 * @param aux_vec auxiliary memory for a vector of the same size as x and y
 * @param args::in_convergence_criterion tolerance to determine if the solution is close enough
 * @param args::in_max_iterations number of iterations to stop at
 * @param args::out_last_error receives the value of the error criterion at the final iteration
 * @param args::out_last_iteration receives the number of the final iteration
 * @param args::opt_error_evolution (optional) pointer to an array of length max_iterations, that receives the error
 * value of each iteration
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_NOT_CONVERGED if it hasn't reached given stopping criterion,
 * in case of failure it returns the associated error code
 */
jmtx_result JMTX_NAME_TYPED(solve_iterative_ilu_crs_precomputed)(const JMTX_NAME_TYPED(matrix_crs) * mtx,
                                                                 const JMTX_NAME_TYPED(matrix_crs) * l,
                                                                 const JMTX_NAME_TYPED(matrix_crs) * u,
                                                                 const JMTX_SCALAR_T *restrict y,
                                                                 JMTX_SCALAR_T *restrict x, JMTX_SCALAR_T *aux_vec,
                                                                 JMTX_NAME_TYPED(solver_arguments) * args);

/**
 * Solves the A x = L U x = y problem by computing the residual, then solving for L U e = r for the error e if residual
 * is too large to further refine the solution. This can help eliminate rounding errors, or alternatively can be used
 * with an incomplete LU decomposition (or ILU) to work as an iterative solver. Must be given the matrices LU.
 * This version offers a minor degree of parallelization, where calculation of the residual and applying error
 * correction are done in parallel, but the main bottleneck of dealing with inverting L U e = r is done in series.
 * @param a system matrix A
 * @param y memory containing forcing vector
 * @param x memory which receives the solution
 * @param aux_vec auxiliary memory for a vector of the same size as x and y
 * @param args::in_convergence_criterion tolerance to determine if the solution is close enough
 * @param args::in_max_iterations number of iterations to stop at
 * @param args::out_last_error receives the value of the error criterion at the final iteration
 * @param args::out_last_iteration receives the number of the final iteration
 * @param args::opt_error_evolution (optional) pointer to an array of length max_iterations, that receives the error
 * value of each iteration
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_NOT_CONVERGED if it hasn't reached given stopping criterion,
 * in case of failure it returns the associated error code
 */
jmtx_result JMTX_NAME_TYPED(solve_iterative_ilu_crs_parallel)(const JMTX_NAME_TYPED(matrix_crs) * mtx,
                                                              const JMTX_SCALAR_T *restrict y,
                                                              JMTX_SCALAR_T *restrict x,
                                                              JMTX_SCALAR_T *restrict aux_vec,
                                                              JMTX_NAME_TYPED(solver_arguments) * args,
                                                              const jmtx_allocator_callbacks *allocator_callbacks);

/**
 * Solves the A x = L U x = y problem by computing the residual, then solving for L U e = r for the error e if residual
 * is too large to further refine the solution. This can help eliminate rounding errors, or alternatively can be used
 * with an incomplete LU decomposition (or ILU) to work as an iterative solver. Must be given the matrices LU.
 * This version offers a minor degree of parallelization, where calculation of the residual and applying error
 * correction are done in parallel, but the main bottleneck of dealing with inverting L U e = r is done in series.
 * @param a system matrix A
 * @param l lower triangular matrix with all entries on its main diagonal equal to 1
 * @param u upper triangular matrix
 * @param y memory containing forcing vector
 * @param x memory which receives the solution
 * @param aux_vec auxiliary memory for a vector of the same size as x and y
 * @param args::in_convergence_criterion tolerance to determine if the solution is close enough
 * @param args::in_max_iterations number of iterations to stop at
 * @param args::out_last_error receives the value of the error criterion at the final iteration
 * @param args::out_last_iteration receives the number of the final iteration
 * @param args::opt_error_evolution (optional) pointer to an array of length max_iterations, that receives the error
 * value of each iteration
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_NOT_CONVERGED if it hasn't reached given stopping criterion,
 * in case of failure it returns the associated error code
 */
jmtx_result JMTX_NAME_TYPED(solve_iterative_ilu_crs_precomputed_parallel)(
    const JMTX_NAME_TYPED(matrix_crs) * mtx, const JMTX_NAME_TYPED(matrix_crs) * l,
    const JMTX_NAME_TYPED(matrix_crs) * u, const JMTX_SCALAR_T *restrict y, JMTX_SCALAR_T *restrict x,
    JMTX_SCALAR_T *restrict aux_vec, JMTX_NAME_TYPED(solver_arguments) * args);

#endif

#ifdef JMTX_DENSE_ROW_MAJOR_H

/**
 * Solves a problem L U x = y, where L is a lower triangular matrix with the diagonal equal to 1 and U is an upper
 * triangular matrix.
 * @param decomposed LU decomposition of the system to solve
 * @param y memory containing forcing vector
 * @param x memory which receives the solution
 * @param aux_vec extra memory of same size as x and y, used to store intermediate values
 */
void JMTX_NAME_TYPED(solve_direct_lu_drm)(const JMTX_NAME_TYPED(matrix_drm) * decomposed,
                                          const JMTX_SCALAR_T *restrict y, JMTX_SCALAR_T *restrict x,
                                          JMTX_SCALAR_T *restrict aux_vec);

#endif

#endif // JMTX_LU_SOLVING_H
