// Automatically generated from tests/cfloat/test_common.h on Fri Dec  1 18:48:10 2023
// Automatically generated from tests/cdouble/test_common.h on Fri Dec  1 17:35:45 2023
//
// Created by jan on 15.6.2022.
//

#ifndef JMTXZ_TEST_COMMON_H
#define JMTXZ_TEST_COMMON_H
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>

struct jmtxz_matrix_crs_struct;
struct jmtxz_matrix_ccs_struct;
struct jmtxz_matrix_brm_struct;
struct jmtxz_matrix_cds_struct;
void print_crs_matrix(const struct jmtxz_matrix_crs_struct* mtx);
void print_ccs_matrix(const struct jmtxz_matrix_ccs_struct* mtx);
void print_brm_matrix(const struct jmtxz_matrix_brm_struct* mtx);
void print_cds_matrix(const struct jmtxz_matrix_cds_struct* mtx);

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

#endif //JMTXZ_TEST_COMMON_H
