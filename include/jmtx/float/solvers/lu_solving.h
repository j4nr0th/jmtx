//
// Created by jan on 6.11.2023.
//

#ifndef JMTX_LU_SOLVING_H
#define JMTX_LU_SOLVING_H
#ifndef JMTX_SOLVER_BASE_H
    #include "../../solver_base.h"
#endif

#ifdef JMTX_BAND_ROW_MAJOR_H

/**
 * Solves a problem L U x = y, where L is a lower triangular matrix with the diagonal equal to 1 and U is an upper
 * triangular matrix.
 * @param l lower triangular matrix with all entries on its main diagonal equal to 1
 * @param u upper triangular matrix
 * @param y memory containing forcing vector
 * @param x memory which receives the solution
 */
void jmtx_solve_direct_lu_brm(const jmtx_matrix_brm* l, const jmtx_matrix_brm* u, const float* restrict y, float* restrict x);

/**
 * Solves a problem L U x = y, where L is a lower triangular matrix with the diagonal equal to 1 and U is an upper
 * triangular matrix. This version of the function stores the solution vector x back into the same memory where the
 * forcing vector was.
 * @param l lower triangular matrix with all entries on its main diagonal equal to 1
 * @param u upper triangular matrix
 * @param x memory which contains the forcing vector and receives the solution
 */
void jmtx_solve_direct_lu_brm_inplace(const jmtx_matrix_brm* l, const jmtx_matrix_brm* u, float* restrict x);

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
jmtx_result jmtxs_solve_direct_lu_brm(const jmtx_matrix_brm* l, const jmtx_matrix_brm* u, uint32_t n,
                               const float y[JMTX_ARRAY_ATTRIB(static restrict n)], float x[JMTX_ARRAY_ATTRIB(static restrict n)]);

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
jmtx_result jmtxs_solve_direct_lu_brm_inplace(const jmtx_matrix_brm* l, const jmtx_matrix_brm* u, uint32_t n, float x[JMTX_ARRAY_ATTRIB(static n)]);

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
jmtx_result jmtx_solve_iterative_lu_brm_refine(const jmtx_matrix_brm* a, const jmtx_matrix_brm* l, const jmtx_matrix_brm* u,
                                        const float y[JMTX_ARRAY_ATTRIB(restrict)], float x[JMTX_ARRAY_ATTRIB(restrict)],
                                        float aux_vec[JMTX_ARRAY_ATTRIB(restrict)], jmtx_solver_arguments* args);


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
jmtx_result jmtxs_solve_iterative_lu_brm_refine(const jmtx_matrix_brm* a, const jmtx_matrix_brm* l, const jmtx_matrix_brm* u,
                                        uint32_t n, const float y[JMTX_ARRAY_ATTRIB(restrict static n)], float x[JMTX_ARRAY_ATTRIB(restrict n)],
                                        float aux_vec[JMTX_ARRAY_ATTRIB(restrict n)], jmtx_solver_arguments* args);

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
jmtx_result jmtx_solve_iterative_lu_brm_refine_parallel(const jmtx_matrix_brm* a, const jmtx_matrix_brm* l,
                                                 const jmtx_matrix_brm* u,  const float y[JMTX_ARRAY_ATTRIB(const restrict)],
                                                 float x[JMTX_ARRAY_ATTRIB(const restrict)], float aux_vec[JMTX_ARRAY_ATTRIB(const restrict)],
                                                 jmtx_solver_arguments* args);
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
void jmtx_solve_direct_lu_crs(const jmtx_matrix_crs* l, const jmtx_matrix_crs* u, const float* restrict y, float* restrict x);

/**
 * Solves a problem L U x = y, where L is a lower triangular matrix with the diagonal equal to 1 and U is an upper
 * triangular matrix. This version of the function stores the solution vector x back into the same memory where the
 * forcing vector was.
 * @param l lower triangular matrix with all entries on its main diagonal equal to 1
 * @param u upper triangular matrix
 * @param x memory which contains the forcing vector and receives the solution
 */
void jmtx_solve_direct_lu_crs_inplace(const jmtx_matrix_crs* l, const jmtx_matrix_crs* u, float* restrict x);

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
jmtx_result jmtxs_solve_direct_lu_crs(const jmtx_matrix_crs* l, const jmtx_matrix_crs* u, uint32_t n,
                               const float y[JMTX_ARRAY_ATTRIB(static restrict n)], float x[JMTX_ARRAY_ATTRIB(restrict n)]);

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
jmtx_result jmtxs_solve_direct_lu_crs_inplace(const jmtx_matrix_crs* l, const jmtx_matrix_crs* u, uint32_t n, float x[JMTX_ARRAY_ATTRIB(static n)]);

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
jmtx_result jmtx_solve_iterative_ilu_crs(
        const jmtx_matrix_crs* mtx, const float* restrict y, float* restrict x, float* aux_vec, jmtx_solver_arguments* args,
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
jmtx_result jmtxs_solve_iterative_ilu_crs(
        const jmtx_matrix_crs* mtx, uint32_t n, const float y[JMTX_ARRAY_ATTRIB(restrict static n)], float x[JMTX_ARRAY_ATTRIB(restrict n)],
        float aux_vec[JMTX_ARRAY_ATTRIB(restrict n)], jmtx_solver_arguments* args, const jmtx_allocator_callbacks* allocator_callbacks);

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
jmtx_result jmtx_solve_iterative_ilu_crs_precomputed(
        const jmtx_matrix_crs* mtx, const jmtx_matrix_crs* l, const jmtx_matrix_crs* u, const float* restrict y, float* restrict x,
        float* aux_vec, jmtx_solver_arguments* args);

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
jmtx_result jmtxs_solve_iterative_ilu_crs_precomputed(
        const jmtx_matrix_crs* mtx, const jmtx_matrix_crs* l, const jmtx_matrix_crs* u, uint32_t n, const float y[JMTX_ARRAY_ATTRIB(restrict static n)],
        float x[JMTX_ARRAY_ATTRIB(restrict n)], float aux_vec[JMTX_ARRAY_ATTRIB(restrict n)], jmtx_solver_arguments* args);

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
jmtx_result jmtx_solve_iterative_ilu_crs_parallel(
        const jmtx_matrix_crs* mtx, const float* restrict y, float* restrict x, float* restrict aux_vec, jmtx_solver_arguments* args,
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
jmtx_result jmtxs_solve_iterative_ilu_crs_parallel(
        const jmtx_matrix_crs* mtx, uint32_t n, const float y[JMTX_ARRAY_ATTRIB(restrict static n)],
        float x[JMTX_ARRAY_ATTRIB(restrict n)], float aux_vec[JMTX_ARRAY_ATTRIB(restrict n)], jmtx_solver_arguments* args,
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
jmtx_result jmtx_solve_iterative_ilu_crs_precomputed_parallel(
        const jmtx_matrix_crs* mtx, const jmtx_matrix_crs* l, const jmtx_matrix_crs* u, const float* restrict y, float* restrict x,
        float* restrict aux_vec, jmtx_solver_arguments* args);

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
jmtx_result jmtxs_solve_iterative_ilu_crs_precomputed_parallel(
        const jmtx_matrix_crs* mtx, const jmtx_matrix_crs* l, const jmtx_matrix_crs* u, uint32_t n,
        const float y[JMTX_ARRAY_ATTRIB(restrict static n)], float x[JMTX_ARRAY_ATTRIB(restrict n)], float aux_vec[JMTX_ARRAY_ATTRIB(restrict n)], jmtx_solver_arguments* args);

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
void jmtx_solve_direct_lu_drm(const jmtx_matrix_drm* decomposed, const float* restrict y, float* restrict x, float* restrict aux_vec);


#endif
#endif //JMTX_LU_SOLVING_H
