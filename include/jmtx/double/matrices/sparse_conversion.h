// Automatically generated from include/jmtx/float/matrices/sparse_conversion.h on Fri Dec  1 06:43:05 2023
//
// Created by jan on 2.11.2023.
//

#ifndef JMTXD_SPARSE_CONVERSION_H
#define JMTXD_SPARSE_CONVERSION_H
#ifndef JMTXD_SPARSE_ROW_COMPRESSED_H
    #include "sparse_row_compressed.h"
#endif
#ifndef JMTXD_SPARSE_COLUMN_COMPRESSED_H
    #include "sparse_column_compressed.h"
#endif
#ifndef JMTXD_SPARSE_DIAGONAL_COMPRESSED_H
    #include "sparse_diagonal_compressed.h"
#endif
#ifndef JMTXD_BAND_ROW_MAJOR_H
    #include "band_row_major.h"
#endif

jmtx_result jmtxd_convert_crs_to_ccs(const jmtxd_matrix_crs* in, jmtxd_matrix_ccs** p_out, const jmtx_allocator_callbacks* allocator_callbacks);

jmtx_result jmtxd_convert_ccs_to_crs(const jmtxd_matrix_ccs* in, jmtxd_matrix_crs** p_out, const jmtx_allocator_callbacks* allocator_callbacks);

jmtx_result jmtxd_convert_crs_to_ccs_inplace_transpose(jmtxd_matrix_crs* in, jmtxd_matrix_ccs** p_out);

jmtx_result jmtxd_convert_ccs_to_crs_inplace_transpose(jmtxd_matrix_ccs* in, jmtxd_matrix_crs** p_out);

jmtx_result jmtxd_convert_cds_to_crs(const jmtxd_matrix_cds* in, jmtxd_matrix_crs** p_out,
                                    const jmtx_allocator_callbacks* allocator_callbacks);

jmtx_result jmtxd_convert_cds_to_ccs(const jmtxd_matrix_cds* in, jmtxd_matrix_ccs** p_out,
                                    const jmtx_allocator_callbacks* allocator_callbacks);

jmtx_result jmtxd_convert_crs_to_cds(const jmtxd_matrix_crs* in, jmtxd_matrix_cds** p_out,
                                    const jmtx_allocator_callbacks* allocator_callbacks);

jmtx_result jmtxd_convert_ccs_to_cds(const jmtxd_matrix_ccs* in, jmtxd_matrix_cds** p_out,
                                    const jmtx_allocator_callbacks* allocator_callbacks);

jmtx_result jmtxd_convert_brm_to_crs(const jmtxd_matrix_brm* in, jmtxd_matrix_crs** p_out,
                                    const jmtx_allocator_callbacks* allocator_callbacks);

jmtx_result jmtxd_convert_brm_to_ccs(const jmtxd_matrix_brm* in, jmtxd_matrix_ccs** p_out,
                                    const jmtx_allocator_callbacks* allocator_callbacks);

jmtx_result jmtxd_convert_brm_to_cds(const jmtxd_matrix_brm* in, jmtxd_matrix_cds** p_out,
                                    const jmtx_allocator_callbacks* allocator_callbacks);

jmtx_result jmtxd_convert_crs_to_brm(const jmtxd_matrix_crs* in, jmtxd_matrix_brm** p_out,
                                    const jmtx_allocator_callbacks* allocator_callbacks);

jmtx_result jmtxd_convert_ccs_to_brm(const jmtxd_matrix_ccs* in, jmtxd_matrix_brm** p_out,
                                    const jmtx_allocator_callbacks* allocator_callbacks);

jmtx_result jmtxd_convert_cds_to_brm(const jmtxd_matrix_cds* in, jmtxd_matrix_brm** p_out,
                                    const jmtx_allocator_callbacks* allocator_callbacks);


#endif //JMTXD_SPARSE_CONVERSION_H
