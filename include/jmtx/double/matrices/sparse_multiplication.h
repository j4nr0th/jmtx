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
jmtx_result jmtxd_matrix_multiply_crs(const jmtxd_matrix_crs* a, const jmtxd_matrix_ccs* b, jmtxd_matrix_crs** p_out,
                                     const jmtx_allocator_callbacks* allocator_callbacks);

jmtx_result jmtxd_matrix_multiply_ccs(const jmtxd_matrix_crs* a, const jmtxd_matrix_ccs* b, jmtxd_matrix_ccs** p_out,
                                     const jmtx_allocator_callbacks* allocator_callbacks);
#endif

double jmtxd_matrix_multiply_sparse_vectors(uint32_t n_a, const uint32_t i_a[static n_a], const double v_a[static n_a],
                                          uint32_t n_b, const uint32_t i_b[static n_b], const double v_b[static n_b]);

double jmtxd_matrix_multiply_sparse_vectors_limit(uint32_t max_a, uint32_t max_b, uint32_t n_a,
                                                const uint32_t i_a[static n_a], const double v_a[static max_a],
                                                uint32_t n_b, const uint32_t i_b[static n_b],
                                                const double v_b[static max_b]);

#ifdef JMTXD_BAND_ROW_MAJOR_H
jmtx_result jmtxd_matrix_multiply_brm(const jmtxd_matrix_brm* a, const jmtxd_matrix_brm* b, jmtxd_matrix_brm** p_out,
                                     const jmtx_allocator_callbacks* allocator_callbacks);
#endif

#ifdef JMTXD_SPARSE_DIAGONAL_COMPRESSED_H
jmtx_result jmtxd_matrix_multiply_cds(const jmtxd_matrix_cds* a, const jmtxd_matrix_cds* b, jmtxd_matrix_cds** p_out,
                                     const jmtx_allocator_callbacks* allocator_callbacks);
#endif

#endif //JMTXD_SPARSE_MULTIPLICATION_H
