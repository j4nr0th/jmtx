#ifndef JMTX_TEST_COMMON_H
#define JMTX_TEST_COMMON_H
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include "../source/common.h"
#include INCLUDE_HEADER
#include "matrices/band_row_major.h"
#include "matrices/dense_row_major.h"
#include "matrices/sparse_column_compressed.h"
#include "matrices/sparse_diagonal_compressed.h"
#include "matrices/sparse_row_compressed.h"

#ifdef JMTX_TYPE_DEFINED_FLOAT

void jmtxf_print_crs_matrix(const jmtxf_matrix_crs *mtx);
void jmtxf_print_ccs_matrix(const jmtxf_matrix_ccs *mtx);
void jmtxf_print_brm_matrix(const jmtxf_matrix_brm *mtx);
void jmtxf_print_cds_matrix(const jmtxf_matrix_cds *mtx);
void jmtxf_print_drm_matrix(const jmtxf_matrix_drm *mtx);

#endif // JMTX_TYPE_DEFINED_FLOAT

#ifdef JMTX_TYPE_DEFINED_DOUBLE

void jmtxd_print_crs_matrix(const jmtxd_matrix_crs *mtx);
void jmtxd_print_ccs_matrix(const jmtxd_matrix_ccs *mtx);
void jmtxd_print_brm_matrix(const jmtxd_matrix_brm *mtx);
void jmtxd_print_cds_matrix(const jmtxd_matrix_cds *mtx);
void jmtxd_print_drm_matrix(const jmtxd_matrix_drm *mtx);

#endif // JMTX_TYPE_DEFINED_DOUBLE

#ifdef JMTX_TYPE_DEFINED_CFLOAT

void jmtxc_print_crs_matrix(const jmtxc_matrix_crs *mtx);
void jmtxc_print_ccs_matrix(const jmtxc_matrix_ccs *mtx);
void jmtxc_print_brm_matrix(const jmtxc_matrix_brm *mtx);
void jmtxc_print_cds_matrix(const jmtxc_matrix_cds *mtx);
void jmtxc_print_drm_matrix(const jmtxc_matrix_drm *mtx);

#endif // JMTX_TYPE_DEFINED_CFLOAT

#ifdef JMTX_TYPE_DEFINED_CDOUBLE

void jmtxz_print_crs_matrix(const jmtxz_matrix_crs *mtx);
void jmtxz_print_ccs_matrix(const jmtxz_matrix_ccs *mtx);
void jmtxz_print_brm_matrix(const jmtxz_matrix_brm *mtx);
void jmtxz_print_cds_matrix(const jmtxz_matrix_cds *mtx);
void jmtxz_print_drm_matrix(const jmtxz_matrix_drm *mtx);

#endif // JMTX_TYPE_DEFINED_CDOUBLE

void jmtxf_print_vec(unsigned n, const float x[JMTX_ARRAY_ATTRIB(static n)]);
void jmtxd_print_vec(unsigned n, const double x[JMTX_ARRAY_ATTRIB(static n)]);
void jmtxc_print_vec(unsigned n, const _Complex float x[JMTX_ARRAY_ATTRIB(static n)]);
void jmtxz_print_vec(unsigned n, const _Complex double x[JMTX_ARRAY_ATTRIB(static n)]);

#ifndef NDEBUG
#    ifdef __GNUC__
#        define DBG_BREAK __builtin_trap()
#    endif
#    ifdef _MSC_BUILD
#        define DBG_BREAK __debugbreak()
#    endif
#else
#    define DBG_BREAK (void)0
#endif

#ifndef ASSERT
#    define ASSERT(x)                                                                                                  \
        (void)(!(x) ? (fprintf(stderr, "Failed assertion \"" #x "\" on line %u, file %s\n", __LINE__, __FILE__),       \
                       DBG_BREAK, exit(EXIT_FAILURE), 1)                                                               \
                    : 0)
#endif // ASSERT
#ifndef MATRIX_TEST_CALL
#    define MATRIX_TEST_CALL(x) printf("Called:\t%s -> %s\n", #x, jmtx_result_to_str((mtx_res = (x))))
#endif // MATRIX_TEST_CALL
#endif // JMTX_TEST_COMMON_H
