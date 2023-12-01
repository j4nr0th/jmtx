// Automatically generated from include/jmtx/float/matrices/sparse_conversion.h on Fri Dec  1 17:35:57 2023
//
// Created by jan on 2.11.2023.
//

#ifndef JMTXC_SPARSE_CONVERSION_H
#define JMTXC_SPARSE_CONVERSION_H
#ifndef JMTXC_SPARSE_ROW_COMPRESSED_H
    #include "sparse_row_compressed.h"
#endif
#ifndef JMTXC_SPARSE_COLUMN_COMPRESSED_H
    #include "sparse_column_compressed.h"
#endif
#ifndef JMTXC_SPARSE_DIAGONAL_COMPRESSED_H
    #include "sparse_diagonal_compressed.h"
#endif
#ifndef JMTXC_BAND_ROW_MAJOR_H
    #include "band_row_major.h"
#endif

jmtx_result jmtxc_convert_crs_to_ccs(const jmtxc_matrix_crs* in, jmtxc_matrix_ccs** p_out,
                                    const jmtx_allocator_callbacks* allocator_callbacks);

jmtx_result jmtxc_convert_ccs_to_crs(const jmtxc_matrix_ccs* in, jmtxc_matrix_crs** p_out,
                                    const jmtx_allocator_callbacks* allocator_callbacks);

jmtx_result jmtxc_convert_crs_to_ccs_inplace_transpose(jmtxc_matrix_crs* in, jmtxc_matrix_ccs** p_out);

jmtx_result jmtxc_convert_ccs_to_crs_inplace_transpose(jmtxc_matrix_ccs* in, jmtxc_matrix_crs** p_out);

jmtx_result jmtxc_convert_cds_to_crs(const jmtxc_matrix_cds* in, jmtxc_matrix_crs** p_out,
                                    const jmtx_allocator_callbacks* allocator_callbacks);

jmtx_result jmtxc_convert_cds_to_ccs(const jmtxc_matrix_cds* in, jmtxc_matrix_ccs** p_out,
                                    const jmtx_allocator_callbacks* allocator_callbacks);

jmtx_result jmtxc_convert_crs_to_cds(const jmtxc_matrix_crs* in, jmtxc_matrix_cds** p_out,
                                    const jmtx_allocator_callbacks* allocator_callbacks);

jmtx_result jmtxc_convert_ccs_to_cds(const jmtxc_matrix_ccs* in, jmtxc_matrix_cds** p_out,
                                    const jmtx_allocator_callbacks* allocator_callbacks);

jmtx_result jmtxc_convert_brm_to_crs(const jmtxc_matrix_brm* in, jmtxc_matrix_crs** p_out,
                                    const jmtx_allocator_callbacks* allocator_callbacks);

jmtx_result jmtxc_convert_brm_to_ccs(const jmtxc_matrix_brm* in, jmtxc_matrix_ccs** p_out,
                                    const jmtx_allocator_callbacks* allocator_callbacks);

jmtx_result jmtxc_convert_brm_to_cds(const jmtxc_matrix_brm* in, jmtxc_matrix_cds** p_out,
                                    const jmtx_allocator_callbacks* allocator_callbacks);

jmtx_result jmtxc_convert_crs_to_brm(const jmtxc_matrix_crs* in, jmtxc_matrix_brm** p_out,
                                    const jmtx_allocator_callbacks* allocator_callbacks);

jmtx_result jmtxc_convert_ccs_to_brm(const jmtxc_matrix_ccs* in, jmtxc_matrix_brm** p_out,
                                    const jmtx_allocator_callbacks* allocator_callbacks);

jmtx_result jmtxc_convert_cds_to_brm(const jmtxc_matrix_cds* in, jmtxc_matrix_brm** p_out,
                                    const jmtx_allocator_callbacks* allocator_callbacks);


#endif //JMTXC_SPARSE_CONVERSION_H
