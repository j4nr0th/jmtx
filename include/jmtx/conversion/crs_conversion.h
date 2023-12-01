//
// Created by jan on 1.12.2023.
//

#ifndef JMTX_CRS_CONVERSION_H
#define JMTX_CRS_CONVERSION_H
#ifndef JMTX_SPARSE_ROW_COMPRESSED_H
    #include "../float/matrices/sparse_row_compressed.h"
#endif
#ifndef JMTXD_SPARSE_ROW_COMPRESSED_H
    #include "../double/matrices/sparse_row_compressed.h"
#endif

/**
 * Creates a new CRS matrix with single precision from a CRS matrix with double precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtx_matrix_crs_from_double(jmtx_matrix_crs** p_mtx, const jmtxd_matrix_crs* in,
                                        const jmtx_allocator_callbacks* allocator_callbacks);

/**
 * Creates a new CRS matrix with single precision from a CRS matrix with double precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxs_matrix_crs_from_double(jmtx_matrix_crs** p_mtx, const jmtxd_matrix_crs* in,
                                        const jmtx_allocator_callbacks* allocator_callbacks);

/**
 * Creates a new CRS matrix with double precision from a CRS matrix with single precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxd_matrix_crs_from_float(jmtxd_matrix_crs** p_mtx, const jmtx_matrix_crs* in,
                                        const jmtx_allocator_callbacks* allocator_callbacks);

/**
 * Creates a new CRS matrix with double precision from a CRS matrix with single precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxds_matrix_crs_from_float(jmtxd_matrix_crs** p_mtx, const jmtx_matrix_crs* in,
                                        const jmtx_allocator_callbacks* allocator_callbacks);

#endif //JMTX_CRS_CONVERSION_H
