//
// Created by jan on 2.11.2023.
//

#ifndef JMTX_SPARSE_MULTIPLICATION_H
#define JMTX_SPARSE_MULTIPLICATION_H
#include "sparse_row_compressed.h"
#include "sparse_column_compressed.h"


jmtx_result jmtx_matrix_multiply_crs(const jmtx_matrix_crs* a, const jmtx_matrix_ccs* b, jmtx_matrix_crs** p_out,
                                     const jmtx_allocator_callbacks* allocator_callbacks);

jmtx_result jmtx_matrix_multiply_ccs(const jmtx_matrix_crs* a, const jmtx_matrix_ccs* b, jmtx_matrix_ccs** p_out,
                                     const jmtx_allocator_callbacks* allocator_callbacks);

#endif //JMTX_SPARSE_MULTIPLICATION_H
