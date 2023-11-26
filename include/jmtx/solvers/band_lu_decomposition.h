//
// Created by jan on 24.11.2023.
//

#ifndef JMTX_BAND_LU_DECOMPOSITION_H
#define JMTX_BAND_LU_DECOMPOSITION_H

#ifndef JMTX_BAND_ROW_MAJOR_H
    #include "../matrices/band_row_major.h"
#endif
/**
 * Uses relations for LU decomposition to compute the full decomposition for the A, such that LU = A. For banded
 * matrices L + U has the same bandwidth as A, so for band matrices, this dramatically reduces costs.
 *
 * For a symmetric SPD matrix, Cholesky factorization should be used instead, which exploits the symmetry of the
 * matrix to give the decomposition in the form of C C^T = A, where C' has the same sparsity pattern as the top of
 * the matrix A.
 *
 * @param a matrix to decompose
 * @param p_l pointer which receives the resulting L BRM matrix
 * @param p_u pointer which receives the resulting U BRM matrix
 * @param allocator_callbacks Pointer to allocators to use for allocating L, U, and auxiliary memory. If NULL, malloc
 * and free are used.
 * @return JMTX_RESULT_SUCCESS if successfully
 */
jmtx_result jmtx_band_lu_decomposition_brm(const jmtx_matrix_brm* a, jmtx_matrix_brm** p_l, jmtx_matrix_brm** p_u,
                                           const jmtx_allocator_callbacks* allocator_callbacks);

#endif //JMTX_BAND_LU_DECOMPOSITION_H
