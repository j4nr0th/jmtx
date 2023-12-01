// Automatically generated from source/cfloat/matrices/basic_io.h on Fri Dec  1 18:48:06 2023
// Automatically generated from source/cdouble/matrices/basic_io.h on Fri Dec  1 17:36:03 2023
//
// Created by jan on 23.11.2023.
//

#ifndef JMTXZ_BASIC_IO_H
#define JMTXZ_BASIC_IO_H
#ifndef JMTXZ_MATRIX_BASE_H
    #include "../../../include/jmtx/matrix_base.h"
#endif
#ifdef JMTXZ_SPARSE_ROW_COMPRESSED_H

jmtx_result jmtxz_matrix_crs_to_file(const jmtxz_matrix_crs* mtx, const char* filename);

jmtx_result jmtxz_matrix_crs_to_file_explicit(const jmtxz_matrix_crs* mtx, const char* filename);

//jmtx_result jmtxz_matrix_crs_from_file(jmtxz_matrix_crs** mtx, const char* filename, jmtx_allocator_callbacks* allocator_callbacks);

#endif


#ifdef JMTXZ_SPARSE_COLUMN_COMPRESSED_H

jmtx_result jmtxz_matrix_ccs_to_file(const jmtxz_matrix_ccs* mtx, const char* filename);

jmtx_result jmtxz_matrix_ccs_to_file_explicit(const jmtxz_matrix_ccs* mtx, const char* filename);

#endif

#endif //JMTXZ_BASIC_IO_H
