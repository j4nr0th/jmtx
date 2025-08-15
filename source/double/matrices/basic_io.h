// Automatically generated from source/float/matrices/basic_io.h on Fri Dec  1 06:43:01 2023
//
// Created by jan on 23.11.2023.
//

#ifndef JMTXD_BASIC_IO_H
#define JMTXD_BASIC_IO_H
#ifndef JMTX_MATRIX_BASE_H
#    include "../../../include/jmtx/matrix_base.h"
#endif
#ifdef JMTXD_SPARSE_ROW_COMPRESSED_H

jmtx_result jmtxd_matrix_crs_to_file(const jmtxd_matrix_crs *mtx, const char *filename);

jmtx_result jmtxd_matrix_crs_to_file_explicit(const jmtxd_matrix_crs *mtx, const char *filename);

// jmtx_result jmtxd_matrix_crs_from_file(jmtxd_matrix_crs** mtx, const char* filename, jmtx_allocator_callbacks*
// allocator_callbacks);

#endif

#ifdef JMTXD_SPARSE_COLUMN_COMPRESSED_H

jmtx_result jmtxd_matrix_ccs_to_file(const jmtxd_matrix_ccs *mtx, const char *filename);

jmtx_result jmtxd_matrix_ccs_to_file_explicit(const jmtxd_matrix_ccs *mtx, const char *filename);

#endif

#endif // JMTXD_BASIC_IO_H
