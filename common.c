//
// Created by jan on 15.6.2022.
//

#include "common.h"
#include "sparse_row_compressed.h"
#include "sparse_column_compressed.h"
#include <stdio.h>


fRNG fRNG_create(uint64_t seed)
{
    uint64_t tmp = (seed += 0x9E3779B97f4A7C15);
    tmp = (tmp ^ (tmp >> 30)) * 0xBF58476D1CE4E5B9;
    tmp = (tmp ^ (tmp >> 27)) * 0x94D049BB133111EB;
    tmp = (tmp ^ (tmp >> 31));

    fRNG rng = {0};
    rng.s0 = tmp;

    tmp = (seed += 0x9E3779B97f4A7C15);
    tmp = (tmp ^ (tmp >> 30)) * 0xBF58476D1CE4E5B9;
    tmp = (tmp ^ (tmp >> 27)) * 0x94D049BB133111EB;
    tmp = (tmp ^ (tmp >> 31));

    rng.s1 = tmp;

    tmp = (seed += 0x9E3779B97f4A7C15);
    tmp = (tmp ^ (tmp >> 30)) * 0xBF58476D1CE4E5B9;
    tmp = (tmp ^ (tmp >> 27)) * 0x94D049BB133111EB;
    tmp = (tmp ^ (tmp >> 31));

    rng.s2 = tmp;

    tmp = (seed += 0x9E3779B97f4A7C15);
    tmp = (tmp ^ (tmp >> 30)) * 0xBF58476D1CE4E5B9;
    tmp = (tmp ^ (tmp >> 27)) * 0x94D049BB133111EB;
    tmp = (tmp ^ (tmp >> 31));

    rng.s3 = tmp;
    return rng;
}

static inline uint64_t roll_left(uint64_t x, uint64_t count)
{
    return (x << count) | (x >> (count - (sizeof(x) << 3)));
}

float fRNG_float(fRNG* p_rng)
{
    const uint64_t res = roll_left(p_rng->s1 * 5, 7) * 9;
    const uint64_t tmp = p_rng->s1 << 17;

    p_rng->s2 ^= p_rng->s0;
    p_rng->s3 ^= p_rng->s1;
    p_rng->s1 ^= p_rng->s2;
    p_rng->s0 ^= p_rng->s3;

    p_rng->s2 ^= tmp;
    p_rng->s3 = roll_left(p_rng->s3, 45);

    return *(float*)&res;
}

double fRNG_double(fRNG* p_rng)
{
    const uint64_t res = roll_left(p_rng->s1 * 5, 7) * 9;
    const uint64_t tmp = p_rng->s1 << 17;

    p_rng->s2 ^= p_rng->s0;
    p_rng->s3 ^= p_rng->s1;
    p_rng->s1 ^= p_rng->s2;
    p_rng->s0 ^= p_rng->s3;

    p_rng->s2 ^= tmp;
    p_rng->s3 = roll_left(p_rng->s3, 45);

    return *(double*)&res;
}

long fRNG_long(fRNG* p_rng)
{
    const uint64_t res = roll_left(p_rng->s1 * 5, 7) * 9;
    const uint64_t tmp = p_rng->s1 << 17;

    p_rng->s2 ^= p_rng->s0;
    p_rng->s3 ^= p_rng->s1;
    p_rng->s1 ^= p_rng->s2;
    p_rng->s0 ^= p_rng->s3;

    p_rng->s2 ^= tmp;
    p_rng->s3 = roll_left(p_rng->s3, 45);

    return *(long*)res;
}

float fRNG_float_range(fRNG* p_rng, float min, float max)
{
    uint64_t res = roll_left(p_rng->s1 * 5, 7) * 9;
    const uint64_t tmp = p_rng->s1 << 17;

    p_rng->s2 ^= p_rng->s0;
    p_rng->s3 ^= p_rng->s1;
    p_rng->s1 ^= p_rng->s2;
    p_rng->s0 ^= p_rng->s3;

    p_rng->s2 ^= tmp;
    p_rng->s3 = roll_left(p_rng->s3, 45);

    return (max - min) * ((float)(res >> 11) / (float)(11llu << 53)) + min;
}

double fRNG_double_range(fRNG* p_rng, double min, double max)
{
    uint64_t res = roll_left(p_rng->s1 * 5, 7) * 9;
    const uint64_t tmp = p_rng->s1 << 17;

    p_rng->s2 ^= p_rng->s0;
    p_rng->s3 ^= p_rng->s1;
    p_rng->s1 ^= p_rng->s2;
    p_rng->s0 ^= p_rng->s3;

    p_rng->s2 ^= tmp;
    p_rng->s3 = roll_left(p_rng->s3, 45);

    return (max - min) * ((double)(res >> 11) / (double)(UINT64_MAX >> 11)) + min;
}

void print_crs_matrix(const CrsMatrix* mtx)
{
    printf("elements:");
    for (uint i = 0, l = 0; i < mtx->n_elements; ++i)
    {
        if (l < mtx->rows && mtx->elements_before[l + 1] <= i + 1)
        {
            l += 1;
        }
        printf(" (%u, %u, %g)", l, mtx->indices[i + 1], mtx->elements[i + 1]);
    }
    printf("\nelement offsets:");
    for (uint i = 0; i < mtx->rows; ++i)
    {
        printf(" %u,", mtx->elements_before[i]);
    }
    printf(" %u", mtx->elements_before[mtx->rows]);
    printf("\nMatrix:\n[");

    for (uint i = 0; i < mtx->rows; ++i)
    {
        if (i != 0)
            printf(" [");
        else
            printf("[");
        scalar_t x;
        for (uint j = 0; j < mtx->columns - 1; ++j)
        {
            matrix_crs_get_element(mtx, i, j, &x);
            printf("%g ", x);
        }
        matrix_crs_get_element(mtx, i, mtx->columns - 1, &x);
        printf("%g] - %u", x, mtx->elements_before[i + 1] - mtx->elements_before[i]);
        if (i != mtx->rows - 1) printf("\n");
    }

    printf("]\n");
}

void print_ccs_matrix(const CcsMatrix* mtx)
{
    printf("elements:");
    for (uint i = 0, l = 0; i < mtx->n_elements; ++i)
    {
        if (l < mtx->rows && mtx->elements_before[l + 1] <= i + 1)
        {
            l += 1;
        }
        printf(" (%u, %u, %g)", l, mtx->indices[i + 1], mtx->elements[i + 1]);
    }
    printf("\nelement offsets:");
    for (uint i = 0; i < mtx->rows; ++i)
    {
        printf(" %u,", mtx->elements_before[i]);
    }
    printf(" %u", mtx->elements_before[mtx->rows]);
    printf("\nMatrix:\n[");

    for (uint i = 0; i < mtx->rows; ++i)
    {
        if (i != 0)
            printf(" [");
        else
            printf("[");
        scalar_t x;
        for (uint j = 0; j < mtx->columns - 1; ++j)
        {
            matrix_ccs_get_element(mtx, i, j, &x);
            printf("%g ", x);
        }
        matrix_ccs_get_element(mtx, i, mtx->columns - 1, &x);
        printf("%g] - %u", x, mtx->elements_before[i + 1] - mtx->elements_before[i]);
        if (i != mtx->rows - 1) printf("\n");
    }

    printf("]\n");
}
