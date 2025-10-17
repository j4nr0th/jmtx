// Automatically generated from include/jmtx/float/matrices/sparse_multiplication.h on Fri Dec  1 06:43:05 2023
//
// Created by jan on 2.11.2023.
//

#ifndef JMTX_SPARSE_MULTIPLICATION_H
#define JMTX_SPARSE_MULTIPLICATION_H
#include "../common.h"

#if defined(JMTX_SPARSE_COLUMN_COMPRESSED_H) && defined(JMTX_SPARSE_ROW_COMPRESSED_H)
/**
 * Multiplies CRS and CCS matrix together and saves the result into a CRS matrix
 * @param a CRS matrix
 * @param b CCS matrix
 * @param p_out pointer which receives the resulting CRS matrix
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on allocation failure
 */
jmtx_result JMTX_NAME_TYPED(multiply_matrix_crs)(const JMTX_NAME_TYPED(matrix_crs) * a,
                                                 const JMTX_NAME_TYPED(matrix_ccs) * b,
                                                 JMTX_NAME_TYPED(matrix_crs) * *p_out,
                                                 const jmtx_allocator_callbacks *allocator_callbacks);

/**
 * Multiplies CRS and CCS matrix together and saves the result into a CCS matrix
 * @param a CRS matrix
 * @param b CCS matrix
 * @param p_out pointer which receives the resulting CCS matrix
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on allocation failure
 */
jmtx_result JMTX_NAME_TYPED(multiply_matrix_ccs)(const JMTX_NAME_TYPED(matrix_crs) * a,
                                                 const JMTX_NAME_TYPED(matrix_ccs) * b,
                                                 JMTX_NAME_TYPED(matrix_ccs) * *p_out,
                                                 const jmtx_allocator_callbacks *allocator_callbacks);
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
JMTX_SCALAR_T JMTX_NAME_TYPED(multiply_matrix_sparse_vectors)(JMTX_INDEX_T n_a,
                                                              const JMTX_INDEX_T i_a[JMTX_ARRAY_ATTRIB(static n_a)],
                                                              const JMTX_SCALAR_T v_a[JMTX_ARRAY_ATTRIB(static n_a)],
                                                              JMTX_INDEX_T n_b,
                                                              const JMTX_INDEX_T i_b[JMTX_ARRAY_ATTRIB(static n_b)],
                                                              const JMTX_SCALAR_T v_b[JMTX_ARRAY_ATTRIB(static n_b)]);
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
JMTX_SCALAR_T JMTX_NAME_TYPED(multiply_matrix_sparse_vectors_limit)(
    JMTX_INDEX_T max_a, JMTX_INDEX_T max_b, JMTX_INDEX_T n_a, const JMTX_INDEX_T i_a[JMTX_ARRAY_ATTRIB(static n_a)],
    const JMTX_SCALAR_T v_a[JMTX_ARRAY_ATTRIB(static max_a)], JMTX_INDEX_T n_b,
    const JMTX_INDEX_T i_b[JMTX_ARRAY_ATTRIB(static n_b)], const JMTX_SCALAR_T v_b[JMTX_ARRAY_ATTRIB(static max_b)]);

#ifdef JMTX_BAND_ROW_MAJOR_H
/**
 * Multiplies two BRM matrices together and produces a BRM matrix with the result of the matrix multiplication.
 * @param a BRM matrix
 * @param b BRM matrix
 * @param p_out pointer which receives the resulting BRM matrix
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on allocation failure
 */
jmtx_result JMTX_NAME_TYPED(multiply_matrix_brm)(const JMTX_NAME_TYPED(matrix_brm) * a,
                                                 const JMTX_NAME_TYPED(matrix_brm) * b,
                                                 JMTX_NAME_TYPED(matrix_brm) * *p_out,
                                                 const jmtx_allocator_callbacks *allocator_callbacks);
#endif

#ifdef JMTX_SPARSE_DIAGONAL_COMPRESSED_H
/**
 * Multiplies two CDS matrices together and produces a CDS matrix with the result of the matrix multiplication.
 * @param a CDS matrix
 * @param b CDS matrix
 * @param p_out pointer which receives the resulting CDS matrix
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on allocation failure
 */
jmtx_result JMTX_NAME_TYPED(multiply_matrix_cds)(const JMTX_NAME_TYPED(matrix_cds) * a,
                                                 const JMTX_NAME_TYPED(matrix_cds) * b,
                                                 JMTX_NAME_TYPED(matrix_cds) * *p_out,
                                                 const jmtx_allocator_callbacks *allocator_callbacks);
#endif

#endif // JMTX_SPARSE_MULTIPLICATION_H
