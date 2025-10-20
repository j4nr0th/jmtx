#ifndef JMTX_CDS_CONVERSION_H
#define JMTX_CDS_CONVERSION_H

/**
 * Creates a new CDS matrix with single precision from a CDS matrix with JMTX_SCALAR_T precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result JMTX_NAME_CONVERSION(matrix_cds)(JMTX_NAME_OUT(matrix_cds) * *p_mtx, const JMTX_NAME_IN(matrix_cds) * in,
                                             const jmtx_allocator_callbacks *allocator_callbacks);
/**
 * Creates a new CDS matrix with single precision from a CDS matrix with JMTX_SCALAR_T precision. Requires no memory
 * allocation by reusing the memory of the initial matrix. Can not fail if the input matrix is valid.
 * @param in matrix which to convert (will be invalid if function succeeds)
 * @return converted matrix
 */
JMTX_NAME_OUT(matrix_cds) * JMTX_NAME_CONVERSION(matrix_cds_inplace)(JMTX_NAME_IN(matrix_cds) * in);

#endif // JMTX_CDS_CONVERSION_H
