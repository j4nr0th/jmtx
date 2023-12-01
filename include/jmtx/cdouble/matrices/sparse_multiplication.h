// Automatically generated from include/jmtx/cfloat/matrices/sparse_multiplication.h on Fri Dec  1 18:48:13 2023
// Automatically generated from include/jmtx/cdouble/matrices/sparse_multiplication.h on Fri Dec  1 17:35:57 2023
//
// Created by jan on 2.11.2023.
//

#ifndef JMTXZ_SPARSE_MULTIPLICATION_H
#define JMTXZ_SPARSE_MULTIPLICATION_H
#ifndef JMTXZ_SPARSE_ROW_COMPRESSED_H
    #include "sparse_row_compressed.h"
#endif
#ifndef JMTXZ_SPARSE_COLUMN_COMPRESSED_H
    #include "sparse_column_compressed.h"
#endif
#ifndef JMTXZ_BAND_ROW_MAJOR_H
    #include "band_row_major.h"
#endif

jmtx_result jmtxz_matrix_multiply_crs(const jmtxz_matrix_crs* a, const jmtxz_matrix_ccs* b, jmtxz_matrix_crs** p_out,
                                     const jmtx_allocator_callbacks* allocator_callbacks);

jmtx_result jmtxz_matrix_multiply_ccs(const jmtxz_matrix_crs* a, const jmtxz_matrix_ccs* b, jmtxz_matrix_ccs** p_out,
                                     const jmtx_allocator_callbacks* allocator_callbacks);

_Complex double jmtxz_matrix_multiply_sparse_vectors(uint32_t n_a, const uint32_t i_a[static n_a], const _Complex double v_a[static n_a],
                                          uint32_t n_b, const uint32_t i_b[static n_b], const _Complex double v_b[static n_b]);

_Complex double jmtxz_matrix_multiply_sparse_vectors_limit(uint32_t max_a, uint32_t max_b, uint32_t n_a,
                                                const uint32_t i_a[static n_a], const _Complex double v_a[static max_a],
                                                uint32_t n_b, const uint32_t i_b[static n_b],
                                                const _Complex double v_b[static max_b]);

jmtx_result jmtxz_matrix_multiply_brm(const jmtxz_matrix_brm* a, const jmtxz_matrix_brm* b, jmtxz_matrix_brm** p_out,
                                     const jmtx_allocator_callbacks* allocator_callbacks);

#endif //JMTXZ_SPARSE_MULTIPLICATION_H
