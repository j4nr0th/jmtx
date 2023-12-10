// Automatically generated from source/float/matrices/basic_io.h on Fri Dec  1 17:36:03 2023
//
// Created by jan on 23.11.2023.
//

#ifndef JMTXC_BASIC_IO_H
#define JMTXC_BASIC_IO_H
#ifndef JMTX_MATRIX_BASE_H
    #include "../../../include/jmtx/matrix_base.h"
#endif
#ifdef JMTXC_SPARSE_ROW_COMPRESSED_H

jmtx_result jmtxc_matrix_crs_to_file(const jmtxc_matrix_crs* mtx, const char* filename);

jmtx_result jmtxc_matrix_crs_to_file_explicit(const jmtxc_matrix_crs* mtx, const char* filename);

//jmtx_result jmtxc_matrix_crs_from_file(jmtxc_matrix_crs** mtx, const char* filename, jmtx_allocator_callbacks* allocator_callbacks);

#endif


#ifdef JMTXC_SPARSE_COLUMN_COMPRESSED_H

jmtx_result jmtxc_matrix_ccs_to_file(const jmtxc_matrix_ccs* mtx, const char* filename);

jmtx_result jmtxc_matrix_ccs_to_file_explicit(const jmtxc_matrix_ccs* mtx, const char* filename);

#endif

#endif //JMTXC_BASIC_IO_H
