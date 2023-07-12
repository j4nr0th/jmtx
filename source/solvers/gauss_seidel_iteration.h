//
// Created by jan on 16.6.2022.
//

#ifndef JMTX_GAUSS_SEIDEL_ITERATION_H
#define JMTX_GAUSS_SEIDEL_ITERATION_H
#include "../matrices/sparse_row_compressed.h"

/*
 * Gauss-Seidel is an iterative method for solving the system Ax = y. It works by splitting the matrix A into
 * matrices D (diagonal), L (lower triangular with zero diagonal), and U (upper triangular with zero diagonal), such
 * that A = D + L + U. The equation is then expressed as x_(n+1) = (D + L)^{-1} (y - U x_(n)). This means that all the
 * functions in this file require that the diagonals of the matrices are non-zero.
 * Note: None of these functions do any safety checks. They are built to be fast and straight forward, not to hold your
 * hand. This will likely be a pain in the ass, but as long as they are tested well enough, it should be fine.
 */

/**
 * Uses Gauss-Seidel (https://en.wikipedia.org/wiki/Gauss%E2%80%93Seidel_method)
 * to solve the linear system Ax = y
 * @param mtx pointer to the memory where matrix A is stored as a compressed row sparse matrix
 * @param y pointer to the memory where the vector y is stored
 * @param x pointer to the memory where the solution vector x will be stored
 * @param convergence_dif when the largest value of change per iteration for an element in x is less than this,
 * the solution is considered found
 * @param n_max_iter maximum number of iterations to perform before giving up
 * @param p_iter pointer which if not null receives the number of iterations that were performed
 * @param p_final_error pointer which receives the final error (may be null)
 * @param p_error pointer to array of error evolution (may be null)
 * @param allocator_callbacks pointer to allocation callbacks used internally by the function (may be null)
 * @return zero if successful
 */
jmtx_result jmtx_gauss_seidel_crs(
        const jmtx_matrix_crs* mtx, const jmtx_scalar_t* y, jmtx_scalar_t* x, jmtx_scalar_t convergence_dif,
        uint32_t n_max_iter, uint32_t* p_iter, jmtx_scalar_t* p_final_error, jmtx_scalar_t* p_error,
        const jmtx_allocator_callbacks* allocator_callbacks);

/**
 * Uses Gauss-Seidel (https://en.wikipedia.org/wiki/Gauss%E2%80%93Seidel_method) to solve the linear system Ax = y.
 * For working in parallel, the Gauss-Seidel line variation is employed. This gives slightly lower convergence, but
 * allows for multiple threads to work in parallel
 * @param mtx pointer to the memory where matrix A is stored as a compressed row sparse matrix
 * @param y pointer to the memory where the vector y is stored
 * @param x pointer to the memory where the solution vector x will be stored
 * @param convergence_dif when the largest value of change per iteration for an element in x is less than this,
 * the solution is considered found
 * @param n_max_iter maximum number of iterations to perform before giving up
 * @param p_iter pointer which if not null receives the number of iterations that were performed
 * @param p_final_error pointer which receives the final error (may be null)
 * @param p_error pointer which receives the error evolution (may be null)
 * @param allocator_callbacks pointer to allocation callbacks used internally by the function (may be null)
 * @param n_thrds the number of threads to use for this (if left as 0, the default number is selected)
 * @return zero if successful
 */
jmtx_result jmtx_gauss_seidel_crs_mt(
        const jmtx_matrix_crs* mtx, const jmtx_scalar_t* y, jmtx_scalar_t* x, jmtx_scalar_t convergence_dif,
        uint32_t n_max_iter, uint32_t* p_iter, jmtx_scalar_t* p_final_error, jmtx_scalar_t* p_error,
        const jmtx_allocator_callbacks* allocator_callbacks, uint32_t n_thrds);

#endif //JMTX_GAUSS_SEIDEL_ITERATION_H
