//
// Created by jan on 2.11.2023.
//

#ifndef JMTXC_SPARSE_MULTIPLICATION_SAFE_H
#define JMTXC_SPARSE_MULTIPLICATION_SAFE_H
#ifndef JMTX_COMMON_H
#include "../../common.h"
#endif

#if defined(JMTXC_SPARSE_ROW_COMPRESSED_H) && defined(JMTXC_SPARSE_COLUMN_COMPRESSED_H)
/**
 * Multiplies CRS and CCS matrix together and saves the result into a CRS matrix
 * @param a CRS matrix
 * @param b CCS matrix
 * @param p_out pointer which receives the resulting CRS matrix
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result jmtxcs_multiply_matrix_crs(const jmtxc_matrix_crs *a, const jmtxc_matrix_ccs *b, jmtxc_matrix_crs **p_out,
                                       const jmtx_allocator_callbacks *allocator_callbacks);

/**
 * Multiplies CRS and CCS matrix together and saves the result into a CCS matrix
 * @param a CRS matrix
 * @param b CCS matrix
 * @param p_out pointer which receives the resulting CCS matrix
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result jmtxcs_multiply_matrix_ccs(const jmtxc_matrix_crs *a, const jmtxc_matrix_ccs *b, jmtxc_matrix_ccs **p_out,
                                       const jmtx_allocator_callbacks *allocator_callbacks);
#endif

#ifdef JMTXC_BAND_ROW_MAJOR_H
/**
 * Multiplies two BRM matrices together and produces a BRM matrix with the result of the matrix multiplication.
 * @param a BRM matrix
 * @param b BRM matrix
 * @param p_out pointer which receives the resulting BRM matrix
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result jmtxcs_multiply_matrix_brm(const jmtxc_matrix_brm *a, const jmtxc_matrix_brm *b, jmtxc_matrix_brm **p_out,
                                       const jmtx_allocator_callbacks *allocator_callbacks);
#endif

#ifdef JMTXC_SPARSE_DIAGONAL_COMPRESSED_H
/**
 * Multiplies two CDS matrices together and produces a CDS matrix with the result of the matrix multiplication.
 * @param a CDS matrix
 * @param b CDS matrix
 * @param p_out pointer which receives the resulting CDS matrix
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result jmtxcs_multiply_matrix_cds(const jmtxc_matrix_cds *a, const jmtxc_matrix_cds *b, jmtxc_matrix_cds **p_out,
                                       const jmtx_allocator_callbacks *allocator_callbacks);
#endif

#endif // JMTXC_SPARSE_MULTIPLICATION_SAFE_H
