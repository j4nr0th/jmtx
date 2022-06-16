//
// Created by jan on 15.6.2022.
//

#ifndef MTXLIB_COMMON_H
#define MTXLIB_COMMON_H
#include <stdint.h>

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

typedef struct struct_CRS_Matrix CrsMatrix;
typedef struct struct_CCS_Matrix CcsMatrix;
void print_crs_matrix(const CrsMatrix* mtx);
void print_ccs_matrix(const CcsMatrix* mtx);

#endif //MTXLIB_COMMON_H
