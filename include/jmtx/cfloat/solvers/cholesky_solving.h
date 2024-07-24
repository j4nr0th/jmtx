// Automatically generated from include/jmtx/float/solvers/cholesky_solving.h on Sun Dec 17 20:14:02 2023
//
// Created by jan on 6.11.2023.
//

#ifndef JMTXC_CHOLESKY_SOLVING_H
#define JMTXC_CHOLESKY_SOLVING_H
#ifndef JMTX_SOLVER_BASE_H
    #include "../../solver_base.h"
#endif

#ifdef JMTXC_SPARSE_ROW_COMPRESSED_H
/**
 * Solves a problem A x = C C^T x = y, where C is a lower triangular matrix.
 * @param c lower triangular matrix C in the CCS format
 * @param ct transpose of the matrix C in the CRS format
 * @param y memory containing forcing vector
 * @param x memory which receives the solution
 */
void jmtxc_solve_direct_cholesky_crs(const jmtxc_matrix_crs* c, const jmtxc_matrix_crs* ct, const _Complex float* restrict y, _Complex float* restrict x);

/**
 * Solves a problem A x = C C^T x = y, where C is a lower triangular matrix. This version of the function stores the
 * solution vector x back into the same memory where the forcing vector was.
 * @param c lower triangular matrix C in the CCS format
 * @param ct transpose of the matrix C in the CRS format
 * @param x memory which contains the forcing vector and receives the solution
 */
void jmtxc_solve_direct_cholesky_crs_inplace(const jmtxc_matrix_crs* c, const jmtxc_matrix_crs* ct, _Complex float* restrict x);

/**
 * Solves a problem A x = C C^T x = y, where C is a lower triangular matrix.
 * @param c lower triangular matrix C in the CCS format
 * @param ct transpose of the matrix C in the CRS format
 * @param y memory containing forcing vector
 * @param x memory which receives the solution
 * @returns JMTX_RESULT_SUCCESS if successful, otherwise an error code indicating error in the input parameters
 */
jmtx_result jmtxcs_solve_direct_cholesky_crs(const jmtxc_matrix_crs* c, const jmtxc_matrix_crs* ct, uint32_t n,
                                 const _Complex float y[JMTX_ARRAY_ATTRIB(static restrict n)], _Complex float x[JMTX_ARRAY_ATTRIB(restrict n)]);

/**
 * Solves a problem A x = C C^T x = y, where C is a lower triangular matrix. This version of the function stores the
 * solution vector x back into the same memory where the forcing vector was.
 * @param c lower triangular matrix C in the CCS format
 * @param ct transpose of the matrix C in the CRS format
 * @param x memory which contains the forcing vector and receives the solution
 * @returns JMTX_RESULT_SUCCESS if successful, otherwise an error code indicating error in the input parameters
 */
jmtx_result jmtxcs_solve_direct_cholesky_crs_inplace(const jmtxc_matrix_crs* c, const jmtxc_matrix_crs* ct, uint32_t n,
                                         _Complex float x[JMTX_ARRAY_ATTRIB(static n)]);

#endif

#endif //JMTXC_CHOLESKY_SOLVING_H
