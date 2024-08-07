// Automatically generated from include/jmtx/float/matrices/sparse_multiplication.h on Fri Dec  1 06:43:05 2023
//
// Created by jan on 2.11.2023.
//

#ifndef JMTXD_SPARSE_MULTIPLICATION_H
#define JMTXD_SPARSE_MULTIPLICATION_H
#ifndef JMTX_COMMON_H
    #include "../../common.h"
#endif

#if defined(JMTXD_SPARSE_COLUMN_COMPRESSED_H) && defined(JMTXD_SPARSE_ROW_COMPRESSED_H)
/**
 * Multiplies CRS and CCS matrix together and saves the result into a CRS matrix
 * @param a CRS matrix
 * @param b CCS matrix
 * @param p_out pointer which receives the resulting CRS matrix
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on allocation failure
 */
jmtx_result jmtxd_multiply_matrix_crs(const jmtxd_matrix_crs* a, const jmtxd_matrix_ccs* b, jmtxd_matrix_crs** p_out,
                                     const jmtx_allocator_callbacks* allocator_callbacks);

/**
 * Multiplies CRS and CCS matrix together and saves the result into a CCS matrix
 * @param a CRS matrix
 * @param b CCS matrix
 * @param p_out pointer which receives the resulting CCS matrix
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on allocation failure
 */
jmtx_result jmtxd_multiply_matrix_ccs(const jmtxd_matrix_crs* a, const jmtxd_matrix_ccs* b, jmtxd_matrix_ccs** p_out,
                                     const jmtx_allocator_callbacks* allocator_callbacks);
#endif

/**
 * Computes the inner product of two sparse vectors
 * @param n_a number of non-zero entries in the first vector
 * @param i_a sorted array of indices of non-zero entries in the first vector
 * @param v_a values of non-zero entries of the first vector
 * @param n_b number of non-zero entries in the second vector
 * @param i_b sorted array of indices of non-zero entries in the second vector
 * @param v_b values of non-zero entries of the second vector
 * @return inner product of the two vectors
 */
double jmtxd_multiply_matrix_sparse_vectors(uint32_t n_a, const uint32_t i_a[JMTX_ARRAY_ATTRIB(static n_a)], const double v_a[JMTX_ARRAY_ATTRIB(static n_a)],
                                          uint32_t n_b, const uint32_t i_b[JMTX_ARRAY_ATTRIB(static n_b)], const double v_b[JMTX_ARRAY_ATTRIB(static n_b)]);
/**
 * Computes the inner product of two sparse vectors, but stops once it reaches a maximum value of the non-zero entry
 * indices (useful when inner product should be done for the first max_a/max_b components)
 * @param max_a maximum non-zero index in the first vector that can be reached
 * @param max_b maximum non-zero index in the second vector that can be reached
 * @param n_a number of non-zero entries in the first vector
 * @param i_a sorted array of indices of non-zero entries in the first vector
 * @param v_a values of non-zero entries of the first vector
 * @param n_b number of non-zero entries in the second vector
 * @param i_b sorted array of indices of non-zero entries in the second vector
 * @param v_b values of non-zero entries of the second vector
 * @return inner product of the two vectors
 */
double jmtxd_multiply_matrix_sparse_vectors_limit(uint32_t max_a, uint32_t max_b, uint32_t n_a,
                                                const uint32_t i_a[JMTX_ARRAY_ATTRIB(static n_a)], const double v_a[JMTX_ARRAY_ATTRIB(static max_a)],
                                                uint32_t n_b, const uint32_t i_b[JMTX_ARRAY_ATTRIB(static n_b)],
                                                const double v_b[JMTX_ARRAY_ATTRIB(static max_b)]);

#ifdef JMTXD_BAND_ROW_MAJOR_H
/**
 * Multiplies two BRM matrices together and produces a BRM matrix with the result of the matrix multiplication.
 * @param a BRM matrix
 * @param b BRM matrix
 * @param p_out pointer which receives the resulting BRM matrix
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on allocation failure
 */
jmtx_result jmtxd_multiply_matrix_brm(const jmtxd_matrix_brm* a, const jmtxd_matrix_brm* b, jmtxd_matrix_brm** p_out,
                                     const jmtx_allocator_callbacks* allocator_callbacks);
#endif

#ifdef JMTXD_SPARSE_DIAGONAL_COMPRESSED_H
/**
 * Multiplies two CDS matrices together and produces a CDS matrix with the result of the matrix multiplication.
 * @param a CDS matrix
 * @param b CDS matrix
 * @param p_out pointer which receives the resulting CDS matrix
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on allocation failure
 */
jmtx_result jmtxd_multiply_matrix_cds(const jmtxd_matrix_cds* a, const jmtxd_matrix_cds* b, jmtxd_matrix_cds** p_out,
                                     const jmtx_allocator_callbacks* allocator_callbacks);
#endif

#endif //JMTXD_SPARSE_MULTIPLICATION_H
