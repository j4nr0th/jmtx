//
// Created by jan on 2.11.2023.
//

#ifndef JMTX_SPARSE_MULTIPLICATION_H
#define JMTX_SPARSE_MULTIPLICATION_H
#ifndef JMTX_SPARSE_ROW_COMPRESSED_H
    #include "sparse_row_compressed.h"
#endif
#ifndef JMTX_SPARSE_COLUMN_COMPRESSED_H
    #include "sparse_column_compressed.h"
#endif
#ifndef JMTX_BAND_ROW_MAJOR_H
    #include "band_row_major.h"
#endif

jmtx_result jmtx_matrix_multiply_crs(const jmtx_matrix_crs* a, const jmtx_matrix_ccs* b, jmtx_matrix_crs** p_out,
                                     const jmtx_allocator_callbacks* allocator_callbacks);

jmtx_result jmtx_matrix_multiply_ccs(const jmtx_matrix_crs* a, const jmtx_matrix_ccs* b, jmtx_matrix_ccs** p_out,
                                     const jmtx_allocator_callbacks* allocator_callbacks);

float jmtx_matrix_multiply_sparse_vectors(uint32_t n_a, const uint32_t i_a[static n_a], const float v_a[static n_a],
                                          uint32_t n_b, const uint32_t i_b[static n_b], const float v_b[static n_b]);

float jmtx_matrix_multiply_sparse_vectors_limit(uint32_t max_a, uint32_t max_b, uint32_t n_a,
                                                const uint32_t i_a[static n_a], const float v_a[static max_a],
                                                uint32_t n_b, const uint32_t i_b[static n_b],
                                                const float v_b[static max_b]);

jmtx_result jmtx_matrix_multiply_brm(const jmtx_matrix_brm* a, const jmtx_matrix_brm* b, jmtx_matrix_brm** p_out,
                                     const jmtx_allocator_callbacks* allocator_callbacks);

#endif //JMTX_SPARSE_MULTIPLICATION_H
