// Automatically generated from include/jmtx/float/solvers/lu_solving.h on Sun Dec 17 17:59:24 2023
//
// Created by jan on 6.11.2023.
//

#ifndef JMTXZ_LU_SOLVING_H
#define JMTXZ_LU_SOLVING_H
#ifndef JMTX_SOLVER_BASE_H
    #include "../../solver_base.h"
#endif

#ifdef JMTXZ_BAND_ROW_MAJOR_H

/**
 * Solves a problem L U x = y, where L is a lower triangular matrix with the diagonal equal to 1 and U is an upper
 * triangular matrix.
 * @param l lower triangular matrix with all entries on its main diagonal equal to 1
 * @param u upper triangular matrix
 * @param y memory containing forcing vector
 * @param x memory which receives the solution
 */
void jmtxz_lu_solve_brm(const jmtxz_matrix_brm* l, const jmtxz_matrix_brm* u, const _Complex double* restrict y, _Complex double* restrict x);

/**
 * Solves a problem L U x = y, where L is a lower triangular matrix with the diagonal equal to 1 and U is an upper
 * triangular matrix. This version of the function stores the solution vector x back into the same memory where the
 * forcing vector was.
 * @param l lower triangular matrix with all entries on its main diagonal equal to 1
 * @param u upper triangular matrix
 * @param x memory which contains the forcing vector and receives the solution
 */
void jmtxz_lu_solve_inplace_brm(const jmtxz_matrix_brm* l, const jmtxz_matrix_brm* u, _Complex double* restrict x);

/**
 * Solves a problem L U x = y, where L is a lower triangular matrix with the diagonal equal to 1 and U is an upper
 * triangular matrix.
 * @param l lower triangular matrix with all entries on its main diagonal equal to 1
 * @param u upper triangular matrix
 * @param n dimensions of L, U, y, and x
 * @param y memory containing forcing vector
 * @param x memory which receives the solution
 * @returns JMTX_RESULT_SUCCESS if successful, otherwise an error code indicating error in the input parameters
 */
jmtx_result jmtxzs_lu_solve_brm(const jmtxz_matrix_brm* l, const jmtxz_matrix_brm* u, uint32_t n,
                               const _Complex double y[static restrict n], _Complex double x[static restrict n]);

/**
 * Solves a problem L U x = y, where L is a lower triangular matrix with the diagonal equal to 1 and U is an upper
 * triangular matrix. This version of the function stores the solution vector x back into the same memory where the
 * forcing vector was.
 * @param l lower triangular matrix with all entries on its main diagonal equal to 1
 * @param u upper triangular matrix
 * @param n dimensions of L, U, y, and x
 * @param x memory which receives the solution
 * @returns JMTX_RESULT_SUCCESS if successful, otherwise an error code indicating error in the input parameters
 */
jmtx_result jmtxzs_lu_solve_inplace_brm(const jmtxz_matrix_brm* l, const jmtxz_matrix_brm* u, uint32_t n, _Complex double x[static n]);

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
 * @param args::opt_error_evolution (optional) pointer to an array of length max_iterations, that receives the error value of each
 * iteration
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_NOT_CONVERGED if it hasn't reached given stopping criterion,
 * in case of failure it returns the associated error code
 */
jmtx_result jmtxz_lu_solve_iterative_bmr(const jmtxz_matrix_brm* a, const jmtxz_matrix_brm* l, const jmtxz_matrix_brm* u,
                                        const _Complex double y[restrict], _Complex double x[restrict],
                                        _Complex double aux_vec[restrict], jmtxd_solver_arguments* args);


/**
 * Solves the A x = L U x = y problem by computing the residual, then solving for L U e = r for the error e if residual
 * is too large to further refine the solution. This can help eliminate rounding errors, or alternatively can be used
 * with an incomplete LU decomposition (or ILU) to work as an iterative solver. Must be given the matrices LU.
 * @param a system matrix A
 * @param l lower triangular matrix with all entries on its main diagonal equal to 1
 * @param u upper triangular matrix
 * @param n dimensions of L, U, y, and x
 * @param y memory containing forcing vector
 * @param x memory which receives the solution
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
jmtx_result jmtxzs_lu_solve_iterative_bmr(const jmtxz_matrix_brm* a, const jmtxz_matrix_brm* l, const jmtxz_matrix_brm* u,
                                        uint32_t n, const _Complex double y[restrict static n], _Complex double x[restrict n],
                                        _Complex double aux_vec[restrict n], jmtxd_solver_arguments* args);

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
 * @param args::opt_error_evolution (optional) pointer to an array of length max_iterations, that receives the error value of each
 * iteration
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_NOT_CONVERGED if it hasn't reached given stopping criterion,
 * in case of failure it returns the associated error code
 */
jmtx_result jmtxz_lu_solve_iterative_bmr_parallel(const jmtxz_matrix_brm* a, const jmtxz_matrix_brm* l,
                                                 const jmtxz_matrix_brm* u,  const _Complex double y[const restrict],
                                                 _Complex double x[const restrict], _Complex double aux_vec[const restrict],
                                                 jmtxd_solver_arguments* args);
#endif
#ifdef JMTXZ_SPARSE_ROW_COMPRESSED_H

/**
 * Solves a problem L U x = y, where L is a lower triangular matrix with the diagonal equal to 1 and U is an upper
 * triangular matrix.
 * @param l lower triangular matrix with all entries on its main diagonal equal to 1
 * @param u upper triangular matrix
 * @param y memory containing forcing vector
 * @param x memory which receives the solution
 */
void jmtxz_lu_solve_crs(const jmtxz_matrix_crs* l, const jmtxz_matrix_crs* u, const _Complex double* restrict y, _Complex double* restrict x);

/**
 * Solves a problem L U x = y, where L is a lower triangular matrix with the diagonal equal to 1 and U is an upper
 * triangular matrix. This version of the function stores the solution vector x back into the same memory where the
 * forcing vector was.
 * @param l lower triangular matrix with all entries on its main diagonal equal to 1
 * @param u upper triangular matrix
 * @param x memory which contains the forcing vector and receives the solution
 */
void jmtxz_lu_solve_inplace_crs(const jmtxz_matrix_crs* l, const jmtxz_matrix_crs* u, _Complex double* restrict x);

/**
 * Solves a problem L U x = y, where L is a lower triangular matrix with the diagonal equal to 1 and U is an upper
 * triangular matrix.
 * @param l lower triangular matrix with all entries on its main diagonal equal to 1
 * @param u upper triangular matrix
 * @param n dimensions of L, U, y, and x
 * @param y memory containing forcing vector
 * @param x memory which receives the solution
 * @returns JMTX_RESULT_SUCCESS if successful, otherwise an error code indicating error in the input parameters
 */
jmtx_result jmtxzs_lu_solve_crs(const jmtxz_matrix_crs* l, const jmtxz_matrix_crs* u, uint32_t n,
                               const _Complex double y[static restrict n], _Complex double x[restrict n]);

/**
 * Solves a problem L U x = y, where L is a lower triangular matrix with the diagonal equal to 1 and U is an upper
 * triangular matrix. This version of the function stores the solution vector x back into the same memory where the
 * forcing vector was.
 * @param l lower triangular matrix with all entries on its main diagonal equal to 1
 * @param u upper triangular matrix
 * @param n dimensions of L, U, y, and x
 * @param x memory which receives the solution
 * @returns JMTX_RESULT_SUCCESS if successful, otherwise an error code indicating error in the input parameters
 */
jmtx_result jmtxzs_lu_solve_inplace_crs(const jmtxz_matrix_crs* l, const jmtxz_matrix_crs* u, uint32_t n, _Complex double x[static n]);

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
 * @param args::opt_error_evolution (optional) pointer to an array of length max_iterations, that receives the error value of each
 * iteration
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_NOT_CONVERGED if it hasn't reached given stopping criterion,
 * in case of failure it returns the associated error code
 */
jmtx_result jmtxz_incomplete_lu_decomposition_solve_crs(
        const jmtxz_matrix_crs* mtx, const _Complex double* restrict y, _Complex double* restrict x, _Complex double* aux_vec, jmtxd_solver_arguments* args,
        const jmtx_allocator_callbacks* allocator_callbacks);

/**
 * Solves the A x = L U x = y problem by computing the residual, then solving for L U e = r for the error e if residual
 * is too large to further refine the solution. This can help eliminate rounding errors, or alternatively can be used
 * with an incomplete LU decomposition (or ILU) to work as an iterative solver. Must be given the matrices LU.
 * @param a system matrix A
 * @param n dimensions of L, U, y, and x
 * @param y memory containing forcing vector
 * @param x memory which receives the solution
 * @param aux_vec auxiliary memory for a vector of the same size as x and y
 * @param args::in_convergence_criterion tolerance to determine if the solution is close enough
 * @param args::in_max_iterations number of iterations to stop at
 * @param args::out_last_error receives the value of the error criterion at the final iteration
 * @param args::out_last_iteration receives the number of the final iteration
 * @param args::opt_error_evolution (optional) pointer to an array of length max_iterations, that receives the error value of each
 * iteration
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_NOT_CONVERGED if it hasn't reached given stopping criterion,
 * in case of failure it returns the associated error code
 */
jmtx_result jmtxzs_incomplete_lu_decomposition_solve_crs(
        const jmtxz_matrix_crs* mtx, uint32_t n, const _Complex double y[restrict static n], _Complex double x[restrict n],
        _Complex double aux_vec[restrict n], jmtxd_solver_arguments* args, const jmtx_allocator_callbacks* allocator_callbacks);

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
 * @param args::opt_error_evolution (optional) pointer to an array of length max_iterations, that receives the error value of each
 * iteration
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_NOT_CONVERGED if it hasn't reached given stopping criterion,
 * in case of failure it returns the associated error code
 */
jmtx_result jmtxz_incomplete_lu_decomposition_solve_precomputed_crs(
        const jmtxz_matrix_crs* mtx, const jmtxz_matrix_crs* l, const jmtxz_matrix_crs* u, const _Complex double* restrict y, _Complex double* restrict x,
        _Complex double* aux_vec, jmtxd_solver_arguments* args);

/**
 * Solves the A x = L U x = y problem by computing the residual, then solving for L U e = r for the error e if residual
 * is too large to further refine the solution. This can help eliminate rounding errors, or alternatively can be used
 * with an incomplete LU decomposition (or ILU) to work as an iterative solver. Must be given the matrices LU.
 * @param a system matrix A
 * @param l lower triangular matrix with all entries on its main diagonal equal to 1
 * @param u upper triangular matrix
 * @param n dimensions of L, U, y, and x
 * @param y memory containing forcing vector
 * @param x memory which receives the solution
 * @param aux_vec auxiliary memory for a vector of the same size as x and y
 * @param args::in_convergence_criterion tolerance to determine if the solution is close enough
 * @param args::in_max_iterations number of iterations to stop at
 * @param args::out_last_error receives the value of the error criterion at the final iteration
 * @param args::out_last_iteration receives the number of the final iteration
 * @param args::opt_error_evolution (optional) pointer to an array of length max_iterations, that receives the error value of each
 * iteration
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_NOT_CONVERGED if it hasn't reached given stopping criterion,
 * in case of failure it returns the associated error code
 */
jmtx_result jmtxzs_incomplete_lu_decomposition_solve_precomputed_crs(
        const jmtxz_matrix_crs* mtx, const jmtxz_matrix_crs* l, const jmtxz_matrix_crs* u, uint32_t n, const _Complex double y[restrict static n],
        _Complex double x[restrict n], _Complex double aux_vec[restrict n], jmtxd_solver_arguments* args);

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
 * @param args::opt_error_evolution (optional) pointer to an array of length max_iterations, that receives the error value of each
 * iteration
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_NOT_CONVERGED if it hasn't reached given stopping criterion,
 * in case of failure it returns the associated error code
 */
jmtx_result jmtxz_incomplete_lu_decomposition_solve_crs_parallel(
        const jmtxz_matrix_crs* mtx, const _Complex double* restrict y, _Complex double* restrict x, _Complex double* restrict aux_vec, jmtxd_solver_arguments* args,
        const jmtx_allocator_callbacks* allocator_callbacks);

/**
 * Solves the A x = L U x = y problem by computing the residual, then solving for L U e = r for the error e if residual
 * is too large to further refine the solution. This can help eliminate rounding errors, or alternatively can be used
 * with an incomplete LU decomposition (or ILU) to work as an iterative solver. Must be given the matrices LU.
 * This version offers a minor degree of parallelization, where calculation of the residual and applying error
 * correction are done in parallel, but the main bottleneck of dealing with inverting L U e = r is done in series.
 * @param a system matrix A
 * @param n dimensions of L, U, y, and x
 * @param y memory containing forcing vector
 * @param x memory which receives the solution
 * @param aux_vec auxiliary memory for a vector of the same size as x and y
 * @param args::in_convergence_criterion tolerance to determine if the solution is close enough
 * @param args::in_max_iterations number of iterations to stop at
 * @param args::out_last_error receives the value of the error criterion at the final iteration
 * @param args::out_last_iteration receives the number of the final iteration
 * @param args::opt_error_evolution (optional) pointer to an array of length max_iterations, that receives the error value of each
 * iteration
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_NOT_CONVERGED if it hasn't reached given stopping criterion,
 * in case of failure it returns the associated error code
 */
jmtx_result jmtxzs_incomplete_lu_decomposition_solve_crs_parallel(
        const jmtxz_matrix_crs* mtx, uint32_t n, const _Complex double y[restrict static n],
        _Complex double x[restrict n], _Complex double aux_vec[restrict n], jmtxd_solver_arguments* args,
        const jmtx_allocator_callbacks* allocator_callbacks);

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
 * @param args::opt_error_evolution (optional) pointer to an array of length max_iterations, that receives the error value of each
 * iteration
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_NOT_CONVERGED if it hasn't reached given stopping criterion,
 * in case of failure it returns the associated error code
 */
jmtx_result jmtxz_incomplete_lu_decomposition_solve_precomputed_crs_parallel(
        const jmtxz_matrix_crs* mtx, const jmtxz_matrix_crs* l, const jmtxz_matrix_crs* u, const _Complex double* restrict y, _Complex double* restrict x,
        _Complex double* restrict aux_vec, jmtxd_solver_arguments* args);

/**
 * Solves the A x = L U x = y problem by computing the residual, then solving for L U e = r for the error e if residual
 * is too large to further refine the solution. This can help eliminate rounding errors, or alternatively can be used
 * with an incomplete LU decomposition (or ILU) to work as an iterative solver. Must be given the matrices LU.
 * This version offers a minor degree of parallelization, where calculation of the residual and applying error
 * correction are done in parallel, but the main bottleneck of dealing with inverting L U e = r is done in series.
 * @param a system matrix A
 * @param l lower triangular matrix with all entries on its main diagonal equal to 1
 * @param u upper triangular matrix
 * @param n dimensions of L, U, y, and x
 * @param y memory containing forcing vector
 * @param x memory which receives the solution
 * @param aux_vec auxiliary memory for a vector of the same size as x and y
 * @param args::in_convergence_criterion tolerance to determine if the solution is close enough
 * @param args::in_max_iterations number of iterations to stop at
 * @param args::out_last_error receives the value of the error criterion at the final iteration
 * @param args::out_last_iteration receives the number of the final iteration
 * @param args::opt_error_evolution (optional) pointer to an array of length max_iterations, that receives the error value of each
 * iteration
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_NOT_CONVERGED if it hasn't reached given stopping criterion,
 * in case of failure it returns the associated error code
 */
jmtx_result jmtxzs_incomplete_lu_decomposition_solve_precomputed_crs_parallel(
        const jmtxz_matrix_crs* mtx, const jmtxz_matrix_crs* l, const jmtxz_matrix_crs* u, uint32_t n,
        const _Complex double y[restrict static n], _Complex double x[restrict n], _Complex double aux_vec[restrict n], jmtxd_solver_arguments* args);

#endif

#endif //JMTXZ_LU_SOLVING_H
