//
// Created by jan on 15.6.2022.
//

#ifndef JMTX_JACOBI_POINT_ITERATION_H
#define JMTX_JACOBI_POINT_ITERATION_H
#include "../matrices/sparse_row_compressed.h"
#include "../matrices/sparse_column_compressed.h"

/*
 * Jacobi Point Iteration is an iterative method for solving the system Ax = y. It works by splitting the matrix A into
 * matrices D (diagonal), L (lower triangular with zero diagonal), and U (upper triangular with zero diagonal), such
 * that A = D + L + U. The equation is then expressed as x_(n+1) = D^{-1} (y - (L + U)x_(n)). This means that all the
 * functions in this file require that the diagonals of the matrices are non-zero.
 * Note: None of these functions do any safety checks. They are built to be fast and straight forward, not to hold your
 * hand. This will likely be a pain in the ass, but as long as they are tested well enough, it should be fine.
 */

/**
 * Uses Jacobi point iteration (also known as Jacobi method: https://en.wikipedia.org/wiki/Jacobi_method)
 * to solve the linear system Ax = y
 *
 * @param mtx pointer to the memory where matrix A is stored as a compressed row sparse matrix
 * @param y pointer to the memory where the vector y is stored
 * @param x pointer to the memory where the solution vector x will be stored
 * @param convergence_dif when the largest value of change per iteration for an element in x is less than this,
 * the solution is considered found
 * @param n_max_iter maximum number of iterations to perform before giving up
 * @param p_iter pointer which if not null receives the number of iterations that were performed
 * @param p_error pointer which receives the evolution of error (may be null)
 * @param p_final_error pointer which receives the final error
 * @param allocator_callbacks allocator callbacks used for internal memory allocations (may be null)
 * @return zero if successful
 */
jmtx_result jmtx_jacobi_crs(
        const jmtx_matrix_crs* mtx, const jmtx_scalar_t* y, jmtx_scalar_t* x, jmtx_scalar_t convergence_dif,
        uint32_t n_max_iter, uint32_t* p_iter, jmtx_scalar_t* p_error, jmtx_scalar_t* p_final_error,
        const jmtx_allocator_callbacks* allocator_callbacks);

/**
 * Uses Jacobi point iteration (also known as Jacobi method: https://en.wikipedia.org/wiki/Jacobi_method)
 * to solve the linear system Ax = y. Uses a relaxation factor ω for the following relation:
 *
 *  x(n + 1) = ω D⁻¹ (y - A x(n)) + x(n)
 *
 *  This may allow for better convergence
 *
 * @param mtx pointer to the memory where matrix A is stored as a compressed row sparse matrix
 * @param y pointer to the memory where the vector y is stored
 * @param x pointer to the memory where the solution vector x will be stored
 * @param convergence_dif when the largest value of change per iteration for an element in x is less than this,
 * the solution is considered found
 * @param n_max_iter maximum number of iterations to perform before giving up
 * @param p_iter pointer which if not null receives the number of iterations that were performed
 * @param p_error pointer which receives the evolution of error (may be null)
 * @param p_final_error pointer which receives the final error
 * @param allocator_callbacks allocator callbacks used for internal memory allocations (may be null)
 * @return zero if successful
 */
jmtx_result jmtx_jacobi_relaxed_crs(
        const jmtx_matrix_crs* mtx, const jmtx_scalar_t* y, jmtx_scalar_t* x, jmtx_scalar_t relaxation_factor,
        jmtx_scalar_t convergence_dif, uint32_t n_max_iter, uint32_t* p_iter, jmtx_scalar_t* p_error,
        jmtx_scalar_t* p_final_error, const jmtx_allocator_callbacks* allocator_callbacks);

/**
 * @warning Currently broken for God knows why. Not in a hurry to fix it, since it is inferior
 *
 * Uses Jacobi point iteration (also known as Jacobi method: https://en.wikipedia.org/wiki/Jacobi_method)
 * to solve the linear system Ax = y
 *
 * @param mtx pointer to the memory where matrix A is stored as a compressed row sparse matrix
 * @param y pointer to the memory where the vector y is stored
 * @param x pointer to the memory where the solution vector x will be stored
 * @param convergence_dif when the largest value of change per iteration for an element in x is less than this,
 * the solution is considered found
 * @param n_max_iter maximum number of iterations to perform before giving up
 * @param p_iter pointer which if not null receives the number of iterations that were performed
 * @param p_error pointer which receives the evolution of error (may be null)
 * @param p_final_error pointer which receives the final error
 * @param allocator_callbacks allocator callbacks used for internal memory allocations (may be null)
 * @param n_thrds the number of threads to use for this (if left as 0, the default number is selected)
 * @return zero if successful
 */
#ifdef __GNUC__
__attribute__((warning("BROKEN, DO NOT USE (YET)")))
#endif
jmtx_result jmtx_jacobi_crs_mt(
        const jmtx_matrix_crs* mtx, const jmtx_scalar_t* y, jmtx_scalar_t* x, jmtx_scalar_t convergence_dif,
        uint32_t n_max_iter, uint32_t* p_iter, jmtx_scalar_t* p_error, jmtx_scalar_t* p_final_error,
        const jmtx_allocator_callbacks* allocator_callbacks, uint32_t n_thrds);


#endif //JMTX_JACOBI_POINT_ITERATION_H
