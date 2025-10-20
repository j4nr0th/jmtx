#ifndef JMTX_DENSE_LU_DECOMPOSITION_H
#define JMTX_DENSE_LU_DECOMPOSITION_H
#include "../matrices/dense_row_major.h"

/**
 * Decomposes a matrix into a lower triangular matrix L and upper triangular matrix U, storing the result in the
 * "decomposed" parameter. The diagonal of the lower triangular matrix is 1 by definition, so it is not stored.
 * Since this decomposition does not need pivoting, it means that decomposition will break down if there's a zero on
 * the diagonal or poorly conditioned.
 * @param mtx square matrix to decompose
 * @param decomposed square matrix which receives the decomposition
 */
void JMTX_NAME_TYPED(decompose_lu_drm)(JMTX_NAME_TYPED(matrix_drm) * mtx, JMTX_NAME_TYPED(matrix_drm) * decomposed);

/**
 * Decomposes a matrix into a lower triangular matrix L and upper triangular matrix U, storing the result in the
 * "decomposed" parameter. The diagonal of the lower triangular matrix is 1 by definition, so it is not stored.
 * This version uses partial pivoting, so the decomposition should not be permuted further in any way, or its
 * permutations committed. Advantage of using partial pivoting is that the decomposition won't break down if there's
 * a zero on a diagonal as well as be less sensitive to rounding errors.
 * @param mtx square matrix to decompose
 * @param decomposed square matrix which receives the decomposition
 * @return JMTX_RESULT_SUCCESS on success, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result JMTX_NAME_TYPED(decompose_lu_pivot_drm)(JMTX_NAME_TYPED(matrix_drm) * mtx,
                                                    JMTX_NAME_TYPED(matrix_drm) * decomposed);

#endif // JMTX_DENSE_LU_DECOMPOSITION_H
