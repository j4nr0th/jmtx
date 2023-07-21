//
// Created by jan on 15.6.2022.
//

#ifndef JMTX_TEST_COMMON_H
#define JMTX_TEST_COMMON_H
#include <stdint.h>
#include <stdlib.h>

typedef union struct_fRNG_state fRNG;

union struct_fRNG_state
{
    uint64_t state[4];
    struct
    {
        uint64_t s0, s1, s2, s3;
    };
};

fRNG fRNG_create(uint64_t seed);

float fRNG_float(fRNG* p_rng);

double fRNG_double(fRNG* p_rng);

long fRNG_long(fRNG* p_rng);

float fRNG_float_range(fRNG* p_rng, float min, float max);

double fRNG_double_range(fRNG* p_rng, double min, double max);

typedef struct jmtx_matrix_crs_struct jmtx_matrix_crs;
typedef struct jmtx_matrix_ccs_struct jmtx_matrix_ccs;
void print_crs_matrix(const jmtx_matrix_crs* mtx);
void print_ccs_matrix(const jmtx_matrix_ccs* mtx);

#ifndef NDEBUG
#   ifdef __GNUC__
#       define DBG_BREAK __builtin_trap()
#   endif
#else
#define DBG_BREAK (void)0
#endif

#define ASSERT(x) if (!(x)) {fprintf(stderr, "Failed assertion \"" #x "\" on line %u, file %s\n", __LINE__, __FILE__); DBG_BREAK; exit(EXIT_FAILURE);} (void)0
#define MATRIX_TEST_CALL(x) printf("Called:\t"#x" -> %s\n", jmtx_result_to_str((mtx_res = (x))))

#endif //JMTX_TEST_COMMON_H
