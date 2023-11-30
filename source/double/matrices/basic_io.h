// Automatically generated from source/float/matrices/basic_io.h on Thu Nov 30 20:33:08 2023
//
// Created by jan on 23.11.2023.
//

#ifndef JMTX_BASIC_IO_H
#define JMTX_BASIC_IO_H
#ifndef JMTX_MATRIX_BASE_H
    #include "../../../include/jmtx/matrix_base.h"
#endif
#ifdef JMTX_SPARSE_ROW_COMPRESSED_H

jmtx_result jmtxd_matrix_crs_to_file(const jmtxd_matrix_crs* mtx, const char* filename);

jmtx_result jmtxd_matrix_crs_to_file_explicit(const jmtxd_matrix_crs* mtx, const char* filename);

//jmtx_result jmtxd_matrix_crs_from_file(jmtxd_matrix_crs** mtx, const char* filename, jmtx_allocator_callbacks* allocator_callbacks);

#endif


#ifdef JMTX_SPARSE_COLUMN_COMPRESSED_H

jmtx_result jmtxd_matrix_ccs_to_file(const jmtxd_matrix_ccs* mtx, const char* filename);

jmtx_result jmtxd_matrix_ccs_to_file_explicit(const jmtxd_matrix_ccs* mtx, const char* filename);

#endif

#endif //JMTX_BASIC_IO_H
