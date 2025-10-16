//
// Created by jan on 23.11.2023.
//

#ifndef JMTX_BASIC_IO_H
#define JMTX_BASIC_IO_H
#ifndef JMTX_MATRIX_BASE_H
#    include "../../../include/jmtx/matrix_base.h"
#endif
#ifdef JMTXF_SPARSE_ROW_COMPRESSED_H

jmtx_result jmtxf_matrix_crs_to_file(const jmtxf_matrix_crs *mtx, const char *filename);

jmtx_result jmtxf_matrix_crs_to_file_explicit(const jmtxf_matrix_crs *mtx, const char *filename);

// jmtx_result jmtxf_matrix_crs_from_file(jmtxf_matrix_crs** mtx, const char* filename, jmtx_allocator_callbacks*
// allocator_callbacks);

#endif

#ifdef JMTXF_SPARSE_COLUMN_COMPRESSED_H

jmtx_result jmtxf_matrix_ccs_to_file(const jmtxf_matrix_ccs *mtx, const char *filename);

jmtx_result jmtxf_matrix_ccs_to_file_explicit(const jmtxf_matrix_ccs *mtx, const char *filename);

#endif

#endif // JMTX_BASIC_IO_H
