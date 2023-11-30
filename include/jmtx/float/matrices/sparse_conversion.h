//
// Created by jan on 2.11.2023.
//

#ifndef JMTX_SPARSE_CONVERSION_H
#define JMTX_SPARSE_CONVERSION_H
#ifndef JMTX_SPARSE_ROW_COMPRESSED_H
    #include "sparse_row_compressed.h"
#endif
#ifndef JMTX_SPARSE_COLUMN_COMPRESSED_H
    #include "sparse_column_compressed.h"
#endif

jmtx_result jmtx_convert_crs_to_ccs(const jmtx_matrix_crs* in, jmtx_matrix_ccs** p_out, const jmtx_allocator_callbacks* allocator_callbacks);

jmtx_result jmtx_convert_ccs_to_crs(const jmtx_matrix_ccs* in, jmtx_matrix_crs** p_out, const jmtx_allocator_callbacks* allocator_callbacks);

jmtx_result jmtx_convert_crs_to_ccs_inplace_transpose(jmtx_matrix_crs* in, jmtx_matrix_ccs** p_out);

jmtx_result jmtx_convert_ccs_to_crs_inplace_transpose(jmtx_matrix_ccs* in, jmtx_matrix_crs** p_out);


#endif //JMTX_SPARSE_CONVERSION_H
