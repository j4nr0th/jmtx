// Automatically generated from include/jmtx/cfloat/matrices/sparse_conversion.h on Fri Dec  1 18:48:13 2023
// Automatically generated from include/jmtx/cdouble/matrices/sparse_conversion.h on Fri Dec  1 17:35:57 2023
//
// Created by jan on 2.11.2023.
//

#ifndef JMTXZ_SPARSE_CONVERSION_H
#define JMTXZ_SPARSE_CONVERSION_H
#ifndef JMTXZ_SPARSE_ROW_COMPRESSED_H
    #include "sparse_row_compressed.h"
#endif
#ifndef JMTXZ_SPARSE_COLUMN_COMPRESSED_H
    #include "sparse_column_compressed.h"
#endif
#ifndef JMTXZ_SPARSE_DIAGONAL_COMPRESSED_H
    #include "sparse_diagonal_compressed.h"
#endif
#ifndef JMTXZ_BAND_ROW_MAJOR_H
    #include "band_row_major.h"
#endif

jmtx_result jmtxz_convert_crs_to_ccs(const jmtxz_matrix_crs* in, jmtxz_matrix_ccs** p_out,
                                    const jmtx_allocator_callbacks* allocator_callbacks);

jmtx_result jmtxz_convert_ccs_to_crs(const jmtxz_matrix_ccs* in, jmtxz_matrix_crs** p_out,
                                    const jmtx_allocator_callbacks* allocator_callbacks);

jmtx_result jmtxz_convert_crs_to_ccs_inplace_transpose(jmtxz_matrix_crs* in, jmtxz_matrix_ccs** p_out);

jmtx_result jmtxz_convert_ccs_to_crs_inplace_transpose(jmtxz_matrix_ccs* in, jmtxz_matrix_crs** p_out);

jmtx_result jmtxz_convert_cds_to_crs(const jmtxz_matrix_cds* in, jmtxz_matrix_crs** p_out,
                                    const jmtx_allocator_callbacks* allocator_callbacks);

jmtx_result jmtxz_convert_cds_to_ccs(const jmtxz_matrix_cds* in, jmtxz_matrix_ccs** p_out,
                                    const jmtx_allocator_callbacks* allocator_callbacks);

jmtx_result jmtxz_convert_crs_to_cds(const jmtxz_matrix_crs* in, jmtxz_matrix_cds** p_out,
                                    const jmtx_allocator_callbacks* allocator_callbacks);

jmtx_result jmtxz_convert_ccs_to_cds(const jmtxz_matrix_ccs* in, jmtxz_matrix_cds** p_out,
                                    const jmtx_allocator_callbacks* allocator_callbacks);

jmtx_result jmtxz_convert_brm_to_crs(const jmtxz_matrix_brm* in, jmtxz_matrix_crs** p_out,
                                    const jmtx_allocator_callbacks* allocator_callbacks);

jmtx_result jmtxz_convert_brm_to_ccs(const jmtxz_matrix_brm* in, jmtxz_matrix_ccs** p_out,
                                    const jmtx_allocator_callbacks* allocator_callbacks);

jmtx_result jmtxz_convert_brm_to_cds(const jmtxz_matrix_brm* in, jmtxz_matrix_cds** p_out,
                                    const jmtx_allocator_callbacks* allocator_callbacks);

jmtx_result jmtxz_convert_crs_to_brm(const jmtxz_matrix_crs* in, jmtxz_matrix_brm** p_out,
                                    const jmtx_allocator_callbacks* allocator_callbacks);

jmtx_result jmtxz_convert_ccs_to_brm(const jmtxz_matrix_ccs* in, jmtxz_matrix_brm** p_out,
                                    const jmtx_allocator_callbacks* allocator_callbacks);

jmtx_result jmtxz_convert_cds_to_brm(const jmtxz_matrix_cds* in, jmtxz_matrix_brm** p_out,
                                    const jmtx_allocator_callbacks* allocator_callbacks);


#endif //JMTXZ_SPARSE_CONVERSION_H
