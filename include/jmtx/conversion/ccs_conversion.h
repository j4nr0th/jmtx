//
// Created by jan on 1.12.2023.
//

#ifndef JMTX_CCS_CONVERSION_H
#define JMTX_CCS_CONVERSION_H



/***********************************************************************************************************************
 *                                                                                                                     *
 *                               Current possibilities for type conversion                                             *
 *                                                                                                                     *
 *                                      FLOAT <---------> DOUBLE                                                       *
 *                                        ^                 ^                                                          *
 *                                        |                 |                                                          *
 *                                        v                 v                                                          *
 *                                  COMPLEX FLOAT <-> COMPLEX DOUBLE                                                   *
 *                                                                                                                     *
 **********************************************************************************************************************/



/***********************************************************************************************************************
 *                                                                                                                     *
 *                                          FLOAT <-> DOUBLE                                                           *
 *                                                                                                                     *
 **********************************************************************************************************************/
#if defined(JMTX_SPARSE_COLUMN_COMPRESSED_H) && defined(JMTXD_SPARSE_COLUMN_COMPRESSED_H)
/**
 * Creates a new CCS matrix with single precision from a CCS matrix with double precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtx_matrix_ccs_from_double(jmtx_matrix_ccs** p_mtx, const jmtxd_matrix_ccs* in,
                                        const jmtx_allocator_callbacks* allocator_callbacks);
/**
 * Creates a new CCS matrix with single precision from a CCS matrix with double precision. Requires no memory
 * allocation by reusing the memory of the initial matrix. Can not fail if the input matrix is valid.
 * @param in matrix which to convert (will be invalid if function succeeds)
 * @return converted matrix
 */
jmtx_matrix_ccs* jmtx_matrix_ccs_from_double_inplace(jmtxd_matrix_ccs* in);

/**
 * Creates a new CCS matrix with single precision from a CCS matrix with double precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxs_matrix_ccs_from_double(jmtx_matrix_ccs** p_mtx, const jmtxd_matrix_ccs* in,
                                        const jmtx_allocator_callbacks* allocator_callbacks);

/**
 * Creates a new CCS matrix with double precision from a CCS matrix with single precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxd_matrix_ccs_from_float(jmtxd_matrix_ccs** p_mtx, const jmtx_matrix_ccs* in,
                                        const jmtx_allocator_callbacks* allocator_callbacks);

/**
 * Creates a new CCS matrix with double precision from a CCS matrix with single precision. Only one memory reallocation
 * may be needed. Can not fail if the input matrix is valid.
 * @param in matrix which to convert
 * @return converted matrix, or NULL in case of allocation failure
 */
jmtxd_matrix_ccs* jmtxd_matrix_ccs_from_float_inplace(jmtx_matrix_ccs* in);

/**
 * Creates a new CCS matrix with double precision from a CCS matrix with single precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxds_matrix_ccs_from_float(jmtxd_matrix_ccs** p_mtx, const jmtx_matrix_ccs* in,
                                        const jmtx_allocator_callbacks* allocator_callbacks);
#endif

/***********************************************************************************************************************
 *                                                                                                                     *
 *                                          FLOAT <-> COMPLEX FLOAT                                                    *
 *                                                                                                                     *
 **********************************************************************************************************************/
#if defined(JMTX_SPARSE_COLUMN_COMPRESSED_H) && defined(JMTXC_SPARSE_COLUMN_COMPRESSED_H) && !defined(JMTX_MSVC)
/**
 * Creates a new CCS matrix with single precision from the real part of a complex CCS matrix with single precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtx_matrix_ccs_from_cfloat_real(jmtx_matrix_ccs** p_mtx, const jmtxc_matrix_ccs* in,
                                        const jmtx_allocator_callbacks* allocator_callbacks);
/**
 * Creates a new CCS matrix with single precision from the imaginary part of a complex CCS matrix with single precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtx_matrix_ccs_from_cfloat_imag(jmtx_matrix_ccs** p_mtx, const jmtxc_matrix_ccs* in,
                                             const jmtx_allocator_callbacks* allocator_callbacks);

/**
 * Creates a new CCS matrix with single precision from the real part of a complex CCS matrix with single precision.
 * Requires no memory allocation by reusing the memory of the initial matrix. Can not fail if the input matrix is valid.
 * @param in matrix which to convert (will be invalid if function succeeds)
 * @return converted matrix
 */
jmtx_matrix_ccs* jmtx_matrix_ccs_from_cfloat_real_inplace(jmtxc_matrix_ccs* in);

/**
 * Creates a new CCS matrix with single precision from the imaginary part of a complex CCS matrix with single precision.
 * Requires no memory allocation by reusing the memory of the initial matrix. Can not fail if the input matrix is valid.
 * @param in matrix which to convert (will be invalid if function succeeds)
 * @return converted matrix
 */
jmtx_matrix_ccs* jmtx_matrix_ccs_from_cfloat_imag_inplace(jmtxc_matrix_ccs* in);

/**
 * Creates a new CCS matrix with single precision from the real part of a complex CCS matrix with single precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxs_matrix_ccs_from_cfloat_real(jmtx_matrix_ccs** p_mtx, const jmtxc_matrix_ccs* in,
                                        const jmtx_allocator_callbacks* allocator_callbacks);

/**
 * Creates a new CCS matrix with single precision from the imaginary part of a complex CCS matrix with single precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxs_matrix_ccs_from_cfloat_imag(jmtx_matrix_ccs** p_mtx, const jmtxc_matrix_ccs* in,
                                              const jmtx_allocator_callbacks* allocator_callbacks);

/**
 * Creates a new complex CCS matrix with single precision from a CCS matrix with single precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in_real matrix to use as the real component of the new matrix (if NULL real part is zero)
 * @param in_imag matrix to use as the real component of the new matrix (if NULL imaginary part is zero)
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxc_matrix_ccs_from_float(jmtxc_matrix_ccs** p_mtx, const jmtx_matrix_ccs* in_real,
                                        const jmtx_matrix_ccs* in_imag,
                                        const jmtx_allocator_callbacks* allocator_callbacks);

/**
 * Creates a new complex CCS matrix with single precision from a CCS matrix with single precision as its real part.
 * Only one memory reallocation.
 * may be needed. Can not fail if the input matrix is valid.
 * @param in matrix which to convert
 * @return converted matrix, or NULL in case of allocation failure
 */
jmtxc_matrix_ccs* jmtxc_matrix_ccs_from_float_real_inplace(jmtx_matrix_ccs* in);

/**
 * Creates a new CCS matrix with single precision from a CCS matrix with single precision as its imaginary part.
 * Only one memory reallocation.
 * may be needed. Can not fail if the input matrix is valid.
 * @param in matrix which to convert
 * @return converted matrix, or NULL in case of allocation failure
 */
jmtxc_matrix_ccs* jmtxc_matrix_ccs_from_float_imag_inplace(jmtx_matrix_ccs* in);

/**
 * Creates a new complex CCS matrix with single precision from a CCS matrix with single precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in_real matrix to use as the real component of the new matrix (if NULL real part is zero)
 * @param in_imag matrix to use as the real component of the new matrix (if NULL imaginary part is zero)
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxcs_matrix_ccs_from_float(jmtxc_matrix_ccs** p_mtx, const jmtx_matrix_ccs* in_real,
                                        const jmtx_matrix_ccs* in_imag,
                                        const jmtx_allocator_callbacks* allocator_callbacks);
#endif

/***********************************************************************************************************************
 *                                                                                                                     *
 *                                          DOUBLE <-> COMPLEX DOUBLE                                                  *
 *                                                                                                                     *
 **********************************************************************************************************************/
#if defined(JMTXD_SPARSE_COLUMN_COMPRESSED_H) && defined(JMTXZ_SPARSE_COLUMN_COMPRESSED_H) && !defined(JMTX_MSVC)
/**
 * Creates a new CCS matrix with double precision from the real part of a complex CCS matrix with double precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxd_matrix_ccs_from_cdouble_real(jmtxd_matrix_ccs** p_mtx, const jmtxz_matrix_ccs* in,
                                        const jmtx_allocator_callbacks* allocator_callbacks);
/**
 * Creates a new CCS matrix with double precision from the imaginary part of a complex CCS matrix with double precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxd_matrix_ccs_from_cdouble_imag(jmtxd_matrix_ccs** p_mtx, const jmtxz_matrix_ccs* in,
                                             const jmtx_allocator_callbacks* allocator_callbacks);

/**
 * Creates a new CCS matrix with double precision from the real part of a complex CCS matrix with double precision.
 * Requires no memory allocation by reusing the memory of the initial matrix. Can not fail if the input matrix is valid.
 * @param in matrix which to convert (will be invalid if function succeeds)
 * @return converted matrix
 */
jmtxd_matrix_ccs* jmtxd_matrix_ccs_from_cdouble_real_inplace(jmtxz_matrix_ccs* in);

/**
 * Creates a new CCS matrix with double precision from the imaginary part of a complex CCS matrix with double precision.
 * Requires no memory allocation by reusing the memory of the initial matrix. Can not fail if the input matrix is valid.
 * @param in matrix which to convert (will be invalid if function succeeds)
 * @return converted matrix
 */
jmtxd_matrix_ccs* jmtxd_matrix_ccs_from_cdouble_imag_inplace(jmtxz_matrix_ccs* in);

/**
 * Creates a new CCS matrix with double precision from the real part of a complex CCS matrix with double precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxds_matrix_ccs_from_cdouble_real(jmtxd_matrix_ccs** p_mtx, const jmtxz_matrix_ccs* in,
                                        const jmtx_allocator_callbacks* allocator_callbacks);

/**
 * Creates a new CCS matrix with double precision from the imaginary part of a complex CCS matrix with double precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxds_matrix_ccs_from_cdouble_imag(jmtxd_matrix_ccs** p_mtx, const jmtxz_matrix_ccs* in,
                                              const jmtx_allocator_callbacks* allocator_callbacks);

/**
 * Creates a new complex CCS matrix with double precision from a CCS matrix with double precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in_real matrix to use as the real component of the new matrix (if NULL real part is zero)
 * @param in_imag matrix to use as the real component of the new matrix (if NULL imaginary part is zero)
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxz_matrix_ccs_from_double(jmtxz_matrix_ccs** p_mtx, const jmtxd_matrix_ccs* in_real,
                                        const jmtxd_matrix_ccs* in_imag,
                                        const jmtx_allocator_callbacks* allocator_callbacks);

/**
 * Creates a new complex CCS matrix with double precision from a CCS matrix with double precision as its real part.
 * Only one memory reallocation.
 * may be needed. Can not fail if the input matrix is valid.
 * @param in matrix which to convert
 * @return converted matrix, or NULL in case of allocation failure
 */
jmtxz_matrix_ccs* jmtxz_matrix_ccs_from_double_real_inplace(jmtxd_matrix_ccs* in);

/**
 * Creates a new CCS matrix with double precision from a CCS matrix with double precision as its imaginary part.
 * Only one memory reallocation.
 * may be needed. Can not fail if the input matrix is valid.
 * @param in matrix which to convert
 * @return converted matrix, or NULL in case of allocation failure
 */
jmtxz_matrix_ccs* jmtxz_matrix_ccs_from_double_imag_inplace(jmtxd_matrix_ccs* in);

/**
 * Creates a new complex CCS matrix with double precision from a CCS matrix with double precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in_real matrix to use as the real component of the new matrix (if NULL real part is zero)
 * @param in_imag matrix to use as the real component of the new matrix (if NULL imaginary part is zero)
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxzs_matrix_ccs_from_double(jmtxz_matrix_ccs** p_mtx, const jmtxd_matrix_ccs* in_real,
                                        const jmtxd_matrix_ccs* in_imag,
                                        const jmtx_allocator_callbacks* allocator_callbacks);
#endif

/***********************************************************************************************************************
 *                                                                                                                     *
 *                                 COMPLEX FLOAT <->  COMPLEX DOUBLE                                                   *
 *                                                                                                                     *
 **********************************************************************************************************************/
#if defined(JMTXC_SPARSE_COLUMN_COMPRESSED_H) && defined(JMTXZ_SPARSE_COLUMN_COMPRESSED_H) && !defined(JMTX_MSVC)
/**
 * Creates a new complex CCS matrix with single precision from a complex CCS matrix with double precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxc_matrix_ccs_from_cdouble(jmtxc_matrix_ccs** p_mtx, const jmtxz_matrix_ccs* in,
                                        const jmtx_allocator_callbacks* allocator_callbacks);
/**
 * Creates a new complex CCS matrix with single precision from a complex CCS matrix with double precision. Requires no memory
 * allocation by reusing the memory of the initial matrix. Can not fail if the input matrix is valid.
 * @param in matrix which to convert (will be invalid if function succeeds)
 * @return converted matrix
 */
jmtxc_matrix_ccs* jmtxc_matrix_ccs_from_cdouble_inplace(jmtxz_matrix_ccs* in);

/**
 * Creates a new complex CCS matrix with single precision from a complex CCS matrix with double precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxcs_matrix_ccs_from_cdouble(jmtxc_matrix_ccs** p_mtx, const jmtxz_matrix_ccs* in,
                                        const jmtx_allocator_callbacks* allocator_callbacks);

/**
 * Creates a new complex CCS matrix with double precision from a complex CCS matrix with single precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxz_matrix_ccs_from_cfloat(jmtxz_matrix_ccs** p_mtx, const jmtxc_matrix_ccs* in,
                                        const jmtx_allocator_callbacks* allocator_callbacks);

/**
 * Creates a new complex CCS matrix with double precision from a complex CCS matrix with single precision. Only one memory reallocation
 * may be needed. Can not fail if the input matrix is valid.
 * @param in matrix which to convert
 * @return converted matrix, or NULL in case of allocation failure
 */
jmtxz_matrix_ccs* jmtxz_matrix_ccs_from_cfloat_inplace(jmtxc_matrix_ccs* in);

/**
 * Creates a new complex CCS matrix with double precision from a complex CCS matrix with single precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxzs_matrix_ccs_from_cfloat(jmtxz_matrix_ccs** p_mtx, const jmtxc_matrix_ccs* in,
                                        const jmtx_allocator_callbacks* allocator_callbacks);
#endif



#endif //JMTX_CCS_CONVERSION_H
