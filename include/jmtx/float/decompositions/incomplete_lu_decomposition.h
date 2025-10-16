//
// Created by jan on 2.11.2023.
//

#ifndef JMTXF_INCOMPLETE_LU_DECOMPOSITION_H
#define JMTXF_INCOMPLETE_LU_DECOMPOSITION_H

#ifndef JMTXF_SPARSE_ROW_COMPRESSED_H
#    include "../matrices/sparse_row_compressed.h"
#endif
#ifndef JMTXF_SPARSE_COLUMN_COMPRESSED_H
#    include "../matrices/sparse_column_compressed.h"
#endif
#ifndef JMTX_SOLVER_BASE_H
#    include "../../solver_base.h"
#endif

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
 * @param allocator_callbacks Pointer to allocators to use for allocating L', U', and auxiliary memory.
 * If NULL, malloc and free are used.
 * @return JMTX_RESULT_SUCCESS if successfully converged in to tolerance in under max iterations,
 * JMTX_RESULT_NOT_CONVERGED if convergence was not achieved in number of specified iterations,
 * other jmtx_result values on other failures.
 */
jmtx_result jmtxf_decompose_ilu_crs(const jmtxf_matrix_crs *a, jmtxf_matrix_crs **p_l, jmtxf_matrix_ccs **p_u,
                                    const jmtx_allocator_callbacks *allocator_callbacks);

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
 * @param allocator_callbacks Pointer to allocators to use for allocating L', U', and auxiliary memory. If NULL, malloc,
 * and free are used.
 * @return JMTX_RESULT_SUCCESS if successfully converged in to tolerance in under max iterations,
 * JMTX_RESULT_NOT_CONVERGED if convergence was not achieved in number of specified iterations,
 * other jmtx_result values on other failures.
 */
jmtx_result jmtxfs_decompose_ilu_crs(const jmtxf_matrix_crs *a, jmtxf_matrix_crs **p_l, jmtxf_matrix_ccs **p_u,
                                     const jmtx_allocator_callbacks *allocator_callbacks);

#endif // JMTXF_INCOMPLETE_LU_DECOMPOSITION_H
