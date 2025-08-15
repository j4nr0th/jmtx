// Automatically generated from include/jmtx/float/solvers/cholesky_solving.h on Sun Dec 17 20:13:54 2023
//
// Created by jan on 6.11.2023.
//

#ifndef JMTXZ_CHOLESKY_SOLVING_H
#define JMTXZ_CHOLESKY_SOLVING_H
#ifndef JMTX_SOLVER_BASE_H
#include "../../solver_base.h"
#endif

#ifdef JMTXZ_SPARSE_ROW_COMPRESSED_H
/**
 * Solves a problem A x = C C^T x = y, where C is a lower triangular matrix.
 * @param c lower triangular matrix C in the CCS format
 * @param ct transpose of the matrix C in the CRS format
 * @param y memory containing forcing vector
 * @param x memory which receives the solution
 */
void jmtxz_solve_direct_cholesky_crs(const jmtxz_matrix_crs *c, const jmtxz_matrix_crs *ct,
                                     const _Complex double *restrict y, _Complex double *restrict x);

/**
 * Solves a problem A x = C C^T x = y, where C is a lower triangular matrix. This version of the function stores the
 * solution vector x back into the same memory where the forcing vector was.
 * @param c lower triangular matrix C in the CCS format
 * @param ct transpose of the matrix C in the CRS format
 * @param x memory which contains the forcing vector and receives the solution
 */
void jmtxz_solve_direct_cholesky_crs_inplace(const jmtxz_matrix_crs *c, const jmtxz_matrix_crs *ct,
                                             _Complex double *restrict x);

/**
 * Solves a problem A x = C C^T x = y, where C is a lower triangular matrix.
 * @param c lower triangular matrix C in the CCS format
 * @param ct transpose of the matrix C in the CRS format
 * @param y memory containing forcing vector
 * @param x memory which receives the solution
 * @returns JMTX_RESULT_SUCCESS if successful, otherwise an error code indicating error in the input parameters
 */
jmtx_result jmtxzs_solve_direct_cholesky_crs(const jmtxz_matrix_crs *c, const jmtxz_matrix_crs *ct, uint32_t n,
                                             const _Complex double y[JMTX_ARRAY_ATTRIB(static restrict n)],
                                             _Complex double x[JMTX_ARRAY_ATTRIB(restrict n)]);

/**
 * Solves a problem A x = C C^T x = y, where C is a lower triangular matrix. This version of the function stores the
 * solution vector x back into the same memory where the forcing vector was.
 * @param c lower triangular matrix C in the CCS format
 * @param ct transpose of the matrix C in the CRS format
 * @param x memory which contains the forcing vector and receives the solution
 * @returns JMTX_RESULT_SUCCESS if successful, otherwise an error code indicating error in the input parameters
 */
jmtx_result jmtxzs_solve_direct_cholesky_crs_inplace(const jmtxz_matrix_crs *c, const jmtxz_matrix_crs *ct, uint32_t n,
                                                     _Complex double x[JMTX_ARRAY_ATTRIB(static n)]);

#endif

#endif // JMTXZ_CHOLESKY_SOLVING_H
