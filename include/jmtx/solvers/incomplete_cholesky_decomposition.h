//
// Created by jan on 2.11.2023.
//

#ifndef JMTX_INCOMPLETE_CHOLESKY_DECOMPOSITION_H
#define JMTX_INCOMPLETE_CHOLESKY_DECOMPOSITION_H

#include "../matrices/sparse_row_compressed.h"
#include "solver_base.h"

/**
 * Uses relations for Cholesky decomposition to compute an approximate decomposition with C' such that the matrix
 * C'C'^T has the same sparsity as the starting matrix. This decomposition can be used as a preconditioner or directly.
 * As with full Cholesky decomposition, the matrix to be decomposed must be SPD. C' is computed so that it has the same
 * sparsity pattern as the matrix A. C' is an lower triangular matrix.
 *
 * @param a matrix to decompose
 * @param p_c pointer which receives the resulting C' CRS matrix
 * @param allocator_callbacks Pointer to allocators to use for allocating L', U', and auxiliary memory. If NULL, malloc
 * and free are used.
 * @return JMTX_RESULT_SUCCESS if successfully converged in to tolerance in under max iterations,
 * JMTX_RESULT_NOT_CONVERGED if convergence was not achieved in number of specified iterations,
 * other jmtx_result values on other failures.
 */
jmtx_result jmtx_incomplete_cholensk_crs(
        const jmtx_matrix_crs* a, jmtx_matrix_crs** p_c, const jmtx_allocator_callbacks* allocator_callbacks);


#endif //JMTX_INCOMPLETE_CHOLESKY_DECOMPOSITION_H