#ifndef JMTX_CCS_CONVERSION_H
#define JMTX_CCS_CONVERSION_H

/**
 * Creates a new CCS matrix with single precision from a CCS matrix with JMTX_SCALAR_T precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result JMTX_NAME_CONVERSION(matrix_ccs)(JMTX_NAME_OUT(matrix_ccs) * *p_mtx, const JMTX_NAME_IN(matrix_ccs) * in,
                                             const jmtx_allocator_callbacks *allocator_callbacks);
/**
 * Creates a new CCS matrix with single precision from a CCS matrix with JMTX_SCALAR_T precision. Requires no memory
 * allocation by reusing the memory of the initial matrix. Can not fail if the input matrix is valid.
 * @param in matrix which to convert (will be invalid if function succeeds)
 * @return converted matrix
 */
JMTX_NAME_OUT(matrix_ccs) * JMTX_NAME_CONVERSION(matrix_ccs_inplace)(JMTX_NAME_IN(matrix_ccs) * in);

#endif // JMTX_CCS_CONVERSION_H
