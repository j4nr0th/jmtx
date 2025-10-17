#ifndef JMTX_CHOLESKY_SOLVING_H
#define JMTX_CHOLESKY_SOLVING_H
#include "../solver_base.h"

#ifdef JMTX_SPARSE_ROW_COMPRESSED_H

/**
 * Solves a problem A x = C C^T x = y, where C is a lower triangular matrix.
 * @param c lower triangular matrix C in the CCS format
 * @param ct transpose of the matrix C in the CRS format
 * @param y memory containing forcing vector
 * @param x memory which receives the solution
 */
void JMTX_NAME_TYPED(solve_direct_cholesky_crs)(const JMTX_NAME_TYPED(matrix_crs) * c,
                                                const JMTX_NAME_TYPED(matrix_crs) * ct, const JMTX_SCALAR_T *restrict y,
                                                JMTX_SCALAR_T *restrict x);

/**
 * Solves a problem A x = C C^T x = y, where C is a lower triangular matrix. This version of the function stores the
 * solution vector x back into the same memory where the forcing vector was.
 * @param c lower triangular matrix C in the CCS format
 * @param ct transpose of the matrix C in the CRS format
 * @param x memory which contains the forcing vector and receives the solution
 */
void JMTX_NAME_TYPED(solve_direct_cholesky_crs_inplace)(const JMTX_NAME_TYPED(matrix_crs) * c,
                                                        const JMTX_NAME_TYPED(matrix_crs) * ct,
                                                        JMTX_SCALAR_T *restrict x);

#endif

#endif // JMTX_CHOLESKY_SOLVING_H
