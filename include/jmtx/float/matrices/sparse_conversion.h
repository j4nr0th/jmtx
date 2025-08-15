//
// Created by jan on 2.11.2023.
//

#ifndef JMTX_SPARSE_CONVERSION_H
#define JMTX_SPARSE_CONVERSION_H

#if defined(JMTX_SPARSE_ROW_COMPRESSED_H) && defined(JMTX_SPARSE_COLUMN_COMPRESSED_H)
/**
 * Converts a CRS matrix into the CCS format. Input matrix remains untouched.
 * @param in CRS matrix to convert
 * @param p_out pointer which receives the converted matrix
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on allocation failure
 */
jmtx_result jmtx_convert_crs_to_ccs(const jmtx_matrix_crs *in, jmtx_matrix_ccs **p_out,
                                    const jmtx_allocator_callbacks *allocator_callbacks);

/**
 * Converts a CCS matrix into the CRS format. Input matrix remains untouched.
 * @param in CCS matrix to convert
 * @param p_out pointer which receives the converted matrix
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on allocation failure
 */
jmtx_result jmtx_convert_ccs_to_crs(const jmtx_matrix_ccs *in, jmtx_matrix_crs **p_out,
                                    const jmtx_allocator_callbacks *allocator_callbacks);
/**
 * Changes the type of a CRS matrix to CCS, which also transposes it.
 * @param in CRS matrix to change the type
 * @return the same pointer as passed to the function, now as CCS matrix
 */
jmtx_matrix_ccs *jmtx_convert_crs_to_ccs_inplace_transpose(jmtx_matrix_crs *in);

/**
 * Changes the type of a CCS matrix to CRS, which also transposes it.
 * @param in CCS matrix to change the type
 * @return the same pointer as passed to the function, now as CCS matrix
 */
jmtx_matrix_crs *jmtx_convert_ccs_to_crs_inplace_transpose(jmtx_matrix_ccs *in);
#endif

#if defined(JMTX_SPARSE_ROW_COMPRESSED_H) && defined(JMTX_SPARSE_DIAGONAL_COMPRESSED_H)
/**
 * Converts a CDS matrix into the CRS format. Input matrix remains untouched.
 * @param in CDS matrix to convert
 * @param p_out pointer which receives the converted matrix
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on allocation failure
 */
jmtx_result jmtx_convert_cds_to_crs(const jmtx_matrix_cds *in, jmtx_matrix_crs **p_out,
                                    const jmtx_allocator_callbacks *allocator_callbacks);

/**
 * Converts a CRS matrix into the CRS format. Input matrix remains untouched.
 * @param in CRS matrix to convert
 * @param p_out pointer which receives the converted matrix
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on allocation failure
 */
jmtx_result jmtx_convert_crs_to_cds(const jmtx_matrix_crs *in, jmtx_matrix_cds **p_out,
                                    const jmtx_allocator_callbacks *allocator_callbacks);
#endif

#if defined(JMTX_SPARSE_COLUMN_COMPRESSED_H) && defined(JMTX_SPARSE_DIAGONAL_COMPRESSED_H)
/**
 * Converts a CDS matrix into the CCS format. Input matrix remains untouched.
 * @param in CDS matrix to convert
 * @param p_out pointer which receives the converted matrix
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on allocation failure
 */
jmtx_result jmtx_convert_cds_to_ccs(const jmtx_matrix_cds *in, jmtx_matrix_ccs **p_out,
                                    const jmtx_allocator_callbacks *allocator_callbacks);

/**
 * Converts a CCS matrix into the CDS format. Input matrix remains untouched.
 * @param in CCS matrix to convert
 * @param p_out pointer which receives the converted matrix
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on allocation failure
 */
jmtx_result jmtx_convert_ccs_to_cds(const jmtx_matrix_ccs *in, jmtx_matrix_cds **p_out,
                                    const jmtx_allocator_callbacks *allocator_callbacks);
#endif

#if defined(JMTX_BAND_ROW_MAJOR_H) && defined(JMTX_SPARSE_ROW_COMPRESSED_H)
/**
 * Converts a BRM matrix into the CRS format. Input matrix remains untouched.
 * @param in BRM matrix to convert
 * @param p_out pointer which receives the converted matrix
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on allocation failure
 */
jmtx_result jmtx_convert_brm_to_crs(const jmtx_matrix_brm *in, jmtx_matrix_crs **p_out,
                                    const jmtx_allocator_callbacks *allocator_callbacks);

/**
 * Converts a CRS matrix into the BRM format. Input matrix remains untouched.
 * @param in CRS matrix to convert
 * @param p_out pointer which receives the converted matrix
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on allocation failure
 */
jmtx_result jmtx_convert_crs_to_brm(const jmtx_matrix_crs *in, jmtx_matrix_brm **p_out,
                                    const jmtx_allocator_callbacks *allocator_callbacks);
#endif

#if defined(JMTX_BAND_ROW_MAJOR_H) && defined(JMTX_SPARSE_COLUMN_COMPRESSED_H)
/**
 * Converts a BRM matrix into the CCS format. Input matrix remains untouched.
 * @param in BRM matrix to convert
 * @param p_out pointer which receives the converted matrix
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on allocation failure
 */
jmtx_result jmtx_convert_brm_to_ccs(const jmtx_matrix_brm *in, jmtx_matrix_ccs **p_out,
                                    const jmtx_allocator_callbacks *allocator_callbacks);

/**
 * Converts a CCS matrix into the BRM format. Input matrix remains untouched.
 * @param in CCS matrix to convert
 * @param p_out pointer which receives the converted matrix
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on allocation failure
 */
jmtx_result jmtx_convert_ccs_to_brm(const jmtx_matrix_ccs *in, jmtx_matrix_brm **p_out,
                                    const jmtx_allocator_callbacks *allocator_callbacks);
#endif

#if defined(JMTX_BAND_ROW_MAJOR_H) && defined(JMTX_SPARSE_DIAGONAL_COMPRESSED_H)
/**
 * Converts a BRM matrix into the CDS format. Input matrix remains untouched.
 * @param in BRM matrix to convert
 * @param p_out pointer which receives the converted matrix
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on allocation failure
 */
jmtx_result jmtx_convert_brm_to_cds(const jmtx_matrix_brm *in, jmtx_matrix_cds **p_out,
                                    const jmtx_allocator_callbacks *allocator_callbacks);

/**
 * Converts a CDS matrix into the BRM format. Input matrix remains untouched.
 * @param in CDS matrix to convert
 * @param p_out pointer which receives the converted matrix
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on allocation failure
 */
jmtx_result jmtx_convert_cds_to_brm(const jmtx_matrix_cds *in, jmtx_matrix_brm **p_out,
                                    const jmtx_allocator_callbacks *allocator_callbacks);
#endif

#endif // JMTX_SPARSE_CONVERSION_H
