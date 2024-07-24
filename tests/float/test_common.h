//
// Created by jan on 15.6.2022.
//

#ifndef JMTX_TEST_COMMON_H
#define JMTX_TEST_COMMON_H
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#ifndef JMTX_COMMON_H
#include "../../include/jmtx/common.h"
#endif//!JMTX_COMMON_H

struct jmtx_matrix_crs_struct;
struct jmtx_matrix_ccs_struct;
struct jmtx_matrix_brm_struct;
struct jmtx_matrix_cds_struct;
struct jmtx_matrix_drm_struct;
void print_crs_matrix(const struct jmtx_matrix_crs_struct* mtx);
void print_ccs_matrix(const struct jmtx_matrix_ccs_struct* mtx);
void print_brm_matrix(const struct jmtx_matrix_brm_struct* mtx);
void print_cds_matrix(const struct jmtx_matrix_cds_struct* mtx);
void print_drm_matrix(const struct jmtx_matrix_drm_struct* mtx);

void print_vec(unsigned n, const float x[JMTX_ARRAY_ATTRIB(static n)]);

#ifndef NDEBUG
#   ifdef __GNUC__
#       define DBG_BREAK __builtin_trap()
#   endif
#   ifdef _MSC_BUILD
#       define DBG_BREAK __debugbreak()
#   endif
#else
#define DBG_BREAK (void)0
#endif

#ifndef ASSERT
    #define ASSERT(x) if (!(x)) {fprintf(stderr, "Failed assertion \"" #x "\" on line %u, file %s\n", __LINE__, __FILE__); DBG_BREAK; exit(EXIT_FAILURE);} (void)0
#endif //ASSERT
#ifndef MATRIX_TEST_CALL
    #define MATRIX_TEST_CALL(x) printf("Called:\t%s -> %s\n", #x, jmtx_result_to_str((mtx_res = (x))))
#endif //MATRIX_TEST_CALL

#endif //JMTX_TEST_COMMON_H
