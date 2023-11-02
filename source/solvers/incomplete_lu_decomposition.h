//
// Created by jan on 2.11.2023.
//

#ifndef JMTX_INCOMPLETE_LU_DECOMPOSITION_H
#define JMTX_INCOMPLETE_LU_DECOMPOSITION_H

#include "../matrices/sparse_row_compressed.h"
#include "../matrices/sparse_column_compressed.h"

/**
 * Uses relations for LU decomposition to compute an approximate decomposition with L' and U' such that the matrix
 * L' + U' has the same sparsity as the starting matrix. This decomposition can be used as a preconditioner or directly.
 *
 * For a symmetric SPD matrix, an incomplete Cholesky factorization is used instead, which exploits the symmetry of the
 * matrix to give the decomposition in the form of C'^T C' = A, where C' has the same sparsity pattern as the top of
 * the matrix A.
 */

jmtx_result jmtx_incomplete_lu_crs(
        jmtx_matrix_crs* a, jmtx_matrix_crs** p_l, jmtx_matrix_ccs** p_u, float convergence, uint32_t max_iterations,
        float* final_max_change, uint32_t* p_last_iteration, const jmtx_allocator_callbacks* allocator_callbacks);



#endif //JMTX_INCOMPLETE_LU_DECOMPOSITION_H
