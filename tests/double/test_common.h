// Automatically generated from tests/float/test_common.h on Fri Dec  1 06:43:09 2023
//
// Created by jan on 15.6.2022.
//

#ifndef JMTXD_TEST_COMMON_H
#define JMTXD_TEST_COMMON_H
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>

typedef struct jmtxd_matrix_crs_struct jmtxd_matrix_crs;
typedef struct jmtxd_matrix_ccs_struct jmtxd_matrix_ccs;
typedef struct jmtxd_matrix_brm_struct jmtxd_matrix_brm;
typedef struct jmtxd_matrix_cds_struct jmtxd_matrix_cds;
void print_crsd_matrix(const jmtxd_matrix_crs* mtx);
void print_ccsd_matrix(const jmtxd_matrix_ccs* mtx);
void print_brmd_matrix(const jmtxd_matrix_brm* mtx);
void print_cdsd_matrix(const jmtxd_matrix_cds* mtx);

#ifndef NDEBUG
#   ifdef __GNUC__
#       define DBG_BREAK __builtin_trap()
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
#endif //JMTXD_TEST_COMMON_H
