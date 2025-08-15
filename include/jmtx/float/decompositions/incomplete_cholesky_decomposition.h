//
// Created by jan on 2.11.2023.
//

#ifndef JMTX_INCOMPLETE_CHOLESKY_DECOMPOSITION_H
#define JMTX_INCOMPLETE_CHOLESKY_DECOMPOSITION_H

#ifndef JMTX_SOLVER_BASE_H
#include "../../solver_base.h"
#endif

#ifdef JMTX_SPARSE_ROW_COMPRESSED_H
/**
 * Uses relations for Cholesky decomposition to compute an approximate decomposition with C' such that the matrix
 * C'C'^T has the same sparsity as the starting matrix. This decomposition can be used as a preconditioner or directly.
 * As with full Cholesky decomposition, the matrix to be decomposed must be SPD. C' is computed so that it has the same
 * sparsity pattern as the matrix A. C' is an lower triangular matrix.
 *
 * @param a matrix to decompose
 * @param p_c pointer which receives the resulting C' CRS matrix
 * @param allocator_callbacks Pointer to allocators to use for allocating C', and auxiliary memory. If NULL, malloc
 * and free are used.
 * @return JMTX_RESULT_SUCCESS if successfully converged in to tolerance in under max iterations,
 * JMTX_RESULT_NOT_CONVERGED if convergence was not achieved in number of specified iterations,
 * other jmtx_result values on other failures.
 */
jmtx_result jmtx_decompose_icho_crs(const jmtx_matrix_crs *a, jmtx_matrix_crs **p_c,
                                    const jmtx_allocator_callbacks *allocator_callbacks);
#endif

#ifdef JMTX_SPARSE_DIAGONAL_COMPRESSED_H
/**
 * Uses relations for Cholesky decomposition to compute an approximate decomposition with C' such that the matrix
 * C'C'^T has the same sparsity as the starting matrix. This decomposition can be used as a preconditioner or directly.
 * As with full Cholesky decomposition, the matrix to be decomposed must be SPD. C' is computed so that it has the same
 * sparsity pattern as the matrix A. C' is an lower triangular matrix.
 *
 * @param a matrix to decompose
 * @param p_c pointer which receives the resulting C' CRS matrix
 * @param allocator_callbacks Pointer to allocators to use for allocating C', and auxiliary memory. If NULL, malloc
 * and free are used.
 * @return JMTX_RESULT_SUCCESS if successfully converged in to tolerance in under max iterations,
 * JMTX_RESULT_NOT_CONVERGED if convergence was not achieved in number of specified iterations,
 * other jmtx_result values on other failures.
 */
jmtx_result jmtx_decompose_icho_cds(const jmtx_matrix_cds *a, jmtx_matrix_cds **p_c,
                                    const jmtx_allocator_callbacks *allocator_callbacks);
#endif

#endif // JMTX_INCOMPLETE_CHOLESKY_DECOMPOSITION_H
