//
// Created by jan on 2.11.2023.
//

#ifndef JMTX_INCOMPLETE_LU_DECOMPOSITION_H
#define JMTX_INCOMPLETE_LU_DECOMPOSITION_H

#include "../matrices/sparse_row_compressed.h"
#include "../matrices/sparse_column_compressed.h"
#include "solver_base.h"

/**
 * Uses relations for LU decomposition to compute an approximate decomposition with L' and U' such that the matrix
 * L'U' has the same sparsity as the starting matrix. This decomposition can be used as a preconditioner or directly.
 *
 * For a symmetric SPD matrix, an incomplete Cholesky factorization is used instead, which exploits the symmetry of the
 * matrix to give the decomposition in the form of C' C'^T = A, where C' has the same sparsity pattern as the top of
 * the matrix A.
 *
 * @param a matrix to decompose
 * @param p_l pointer which receives the resulting L' CRS matrix
 * @param p_u pointer which receives the resulting U' CCS matrix
 * @param allocator_callbacks Pointer to allocators to use for allocating L', U', and auxiliary memory. If NULL, malloc
 * and free are used.
 * @return JMTX_RESULT_SUCCESS if successfully converged in to tolerance in under max iterations,
 * JMTX_RESULT_NOT_CONVERGED if convergence was not achieved in number of specified iterations,
 * other jmtx_result values on other failures.
 */
jmtx_result jmtx_incomplete_lu_crs(
        const jmtx_matrix_crs* a, jmtx_matrix_crs** p_l, jmtx_matrix_ccs** p_u,
        const jmtx_allocator_callbacks* allocator_callbacks);

/**
 * Uses relations for LU decomposition to compute an approximate decomposition with L' and U' such that the matrix
 * L'U' has the same sparsity as the starting matrix. This decomposition can be used as a preconditioner or directly.
 *
 * For a symmetric SPD matrix, an incomplete Cholesky factorization is used instead, which exploits the symmetry of the
 * matrix to give the decomposition in the form of C' C'^T = A, where C' has the same sparsity pattern as the top of
 * the matrix A.
 *
 * Parallel version of the function, which is parallelized. Best case it performs as serial version, behaving like Gauss
 * Seidel, worst case, it becomes the same as if point Jacobi iterations were done instead.
 *
 * @param a matrix to decompose
 * @param p_l pointer which receives the resulting L' CRS matrix
 * @param p_u pointer which receives the resulting U' CCS matrix
 * @param allocator_callbacks Pointer to allocators to use for allocating L', U', and auxiliary memory. Does not need to
 * be thread-safe. If NULL, malloc and free are used.
 * @param args::in_convergence_criterion Stopping criterion to find if more iterations are needed. This is done based on
 * the fraction between the change of an element and its magnitude. When the largest value of that fraction falls bellow
 * this value, iterations stop.
 * @param args::in_max_iterations Maximum number of iterations to perform.
 * @param args::out_last_error Receives the largest value of the stopping criterion on the last iteration.
 * @param args::out_last_iteration Receives the number of the last iteration.
 * @param args::opt_error_evolution Receives the evolution of error.
 * @return JMTX_RESULT_SUCCESS if successfully converged in to tolerance in under max iterations,
 * JMTX_RESULT_NOT_CONVERGED if convergence was not achieved in number of specified iterations,
 * other jmtx_result values on other failures.
 */
jmtx_result jmtx_incomplete_lu_crs_parallel(
        jmtx_matrix_crs* a, jmtx_matrix_crs** p_l, jmtx_matrix_ccs** p_u,
        const jmtx_allocator_callbacks* allocator_callbacks, jmtx_solver_arguments* args);


jmtx_result jmtx_incomplete_lu_crs_old(
        const jmtx_matrix_crs* a, jmtx_matrix_crs** p_l, jmtx_matrix_ccs** p_u, float convergence, uint32_t max_iterations,
        float* final_max_change, uint32_t* p_last_iteration, const jmtx_allocator_callbacks* allocator_callbacks);


#endif //JMTX_INCOMPLETE_LU_DECOMPOSITION_H
