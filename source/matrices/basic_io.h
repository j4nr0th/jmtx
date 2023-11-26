//
// Created by jan on 23.11.2023.
//

#ifndef JMTX_BASIC_IO_H
#define JMTX_BASIC_IO_H
#ifndef JMTX_MATRIX_BASE_H
    #include "../../include/jmtx/matrices/matrix_base.h"
#endif
#ifdef JMTX_SPARSE_ROW_COMPRESSED_H

jmtx_result jmtx_matrix_crs_to_file(const jmtx_matrix_crs* mtx, const char* filename);

jmtx_result jmtx_matrix_crs_to_file_explicit(const jmtx_matrix_crs* mtx, const char* filename);

//jmtx_result jmtx_matrix_crs_from_file(jmtx_matrix_crs** mtx, const char* filename, jmtx_allocator_callbacks* allocator_callbacks);

#endif


#ifdef JMTX_SPARSE_COLUMN_COMPRESSED_H

jmtx_result jmtx_matrix_ccs_to_file(const jmtx_matrix_ccs* mtx, const char* filename);

jmtx_result jmtx_matrix_ccs_to_file_explicit(const jmtx_matrix_ccs* mtx, const char* filename);

#endif

#endif //JMTX_BASIC_IO_H
