// Automatically generated from include/jmtx/float/matrices/sparse_multiplication.h on Fri Dec  1 17:35:57 2023
//
// Created by jan on 2.11.2023.
//

#ifndef JMTXC_SPARSE_MULTIPLICATION_H
#define JMTXC_SPARSE_MULTIPLICATION_H
#ifndef JMTXC_SPARSE_ROW_COMPRESSED_H
    #include "sparse_row_compressed.h"
#endif
#ifndef JMTXC_SPARSE_COLUMN_COMPRESSED_H
    #include "sparse_column_compressed.h"
#endif
#ifndef JMTXC_BAND_ROW_MAJOR_H
    #include "band_row_major.h"
#endif

jmtx_result jmtxc_matrix_multiply_crs(const jmtxc_matrix_crs* a, const jmtxc_matrix_ccs* b, jmtxc_matrix_crs** p_out,
                                     const jmtx_allocator_callbacks* allocator_callbacks);

jmtx_result jmtxc_matrix_multiply_ccs(const jmtxc_matrix_crs* a, const jmtxc_matrix_ccs* b, jmtxc_matrix_ccs** p_out,
                                     const jmtx_allocator_callbacks* allocator_callbacks);

_Complex float jmtxc_matrix_multiply_sparse_vectors(uint32_t n_a, const uint32_t i_a[static n_a], const _Complex float v_a[static n_a],
                                          uint32_t n_b, const uint32_t i_b[static n_b], const _Complex float v_b[static n_b]);

_Complex float jmtxc_matrix_multiply_sparse_vectors_limit(uint32_t max_a, uint32_t max_b, uint32_t n_a,
                                                const uint32_t i_a[static n_a], const _Complex float v_a[static max_a],
                                                uint32_t n_b, const uint32_t i_b[static n_b],
                                                const _Complex float v_b[static max_b]);

jmtx_result jmtxc_matrix_multiply_brm(const jmtxc_matrix_brm* a, const jmtxc_matrix_brm* b, jmtxc_matrix_brm** p_out,
                                     const jmtx_allocator_callbacks* allocator_callbacks);

#endif //JMTXC_SPARSE_MULTIPLICATION_H
