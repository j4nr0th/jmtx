//
// Created by jan on 6.11.2023.
//

#ifndef JMTX_CHOLESKY_SOLVING_H
#define JMTX_CHOLESKY_SOLVING_H
#ifndef JMTX_SOLVER_BASE_H
    #include "../../solver_base.h"
#endif

#ifdef JMTX_SPARSE_ROW_COMPRESSED_H
/**
 * Solves a problem A x = C C^T x = y, where C is a lower triangular matrix.
 * @param c lower triangular matrix C in the CCS format
 * @param ct transpose of the matrix C in the CRS format
 * @param y memory containing forcing vector
 * @param x memory which receives the solution
 */
void jmtx_solve_direct_cholesky_crs(const jmtx_matrix_crs* c, const jmtx_matrix_crs* ct, const float* restrict y, float* restrict x);

/**
 * Solves a problem A x = C C^T x = y, where C is a lower triangular matrix. This version of the function stores the
 * solution vector x back into the same memory where the forcing vector was.
 * @param c lower triangular matrix C in the CCS format
 * @param ct transpose of the matrix C in the CRS format
 * @param x memory which contains the forcing vector and receives the solution
 */
void jmtx_solve_direct_cholesky_crs_inplace(const jmtx_matrix_crs* c, const jmtx_matrix_crs* ct, float* restrict x);

/**
 * Solves a problem A x = C C^T x = y, where C is a lower triangular matrix.
 * @param c lower triangular matrix C in the CCS format
 * @param ct transpose of the matrix C in the CRS format
 * @param y memory containing forcing vector
 * @param x memory which receives the solution
 * @returns JMTX_RESULT_SUCCESS if successful, otherwise an error code indicating error in the input parameters
 */
jmtx_result jmtxs_solve_direct_cholesky_crs(const jmtx_matrix_crs* c, const jmtx_matrix_crs* ct, uint32_t n,
                                 const float y[JMTX_ARRAY_ATTRIB(static restrict n)], float x[JMTX_ARRAY_ATTRIB(restrict n)]);

/**
 * Solves a problem A x = C C^T x = y, where C is a lower triangular matrix. This version of the function stores the
 * solution vector x back into the same memory where the forcing vector was.
 * @param c lower triangular matrix C in the CCS format
 * @param ct transpose of the matrix C in the CRS format
 * @param x memory which contains the forcing vector and receives the solution
 * @returns JMTX_RESULT_SUCCESS if successful, otherwise an error code indicating error in the input parameters
 */
jmtx_result jmtxs_solve_direct_cholesky_crs_inplace(const jmtx_matrix_crs* c, const jmtx_matrix_crs* ct, uint32_t n,
                                         float x[JMTX_ARRAY_ATTRIB(static n)]);

#endif

#endif //JMTX_CHOLESKY_SOLVING_H
