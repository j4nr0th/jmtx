#ifndef JMTX_BASIC_IO_H
#define JMTX_BASIC_IO_H
#include "../matrix_base.h"

jmtx_result JMTX_NAME_TYPED(matrix_crs_to_file)(const JMTX_NAME_TYPED(matrix_crs) * mtx, const char *filename);

jmtx_result JMTX_NAME_TYPED(matrix_crs_to_file_explicit)(const JMTX_NAME_TYPED(matrix_crs) * mtx, const char *filename);

jmtx_result JMTX_NAME_TYPED(matrix_ccs_to_file)(const JMTX_NAME_TYPED(matrix_ccs) * mtx, const char *filename);

jmtx_result JMTX_NAME_TYPED(matrix_ccs_to_file_explicit)(const JMTX_NAME_TYPED(matrix_ccs) * mtx, const char *filename);

#endif // JMTX_BASIC_IO_H
