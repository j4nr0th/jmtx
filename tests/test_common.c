//
// Created by jan on 15.6.2022.
//

#include "test_common.h"
#include "../source/matrices/sparse_row_compressed.h"
#include "../source/matrices/sparse_column_compressed.h"
#include <stdio.h>
#include "../source/matrices/sparse_row_compressed_internal.h"


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

void print_crs_matrix(const jmtx_matrix_crs* mtx)
{
    printf("Entries:");
    for (uint32_t i = 0, l = 0; i < mtx->n_entries; ++i)
    {
        while (l < mtx->base.rows && mtx->end_of_row_offsets[l] <= i)
        {
            l += 1;
        }
        printf(" (%u, %u, %g)", l, mtx->indices[i], mtx->values[i]);
    }
    printf("\nEnd of row offsets:");
    for (uint32_t i = 0; i < mtx->base.rows; ++i)
    {
        printf(" %u,", mtx->end_of_row_offsets[i]);
    }
    printf("\nMatrix:\n[\n");

    for (uint32_t i = 0; i < mtx->base.rows; ++i)
    {
        printf("\t[");
        float x;
        for (uint32_t j = 0; j < mtx->base.cols; ++j)
        {
            jmtx_matrix_crs_get_element(mtx, i, j, &x);
            printf("%g ", x);
        }
        printf("] - %u", mtx->end_of_row_offsets[i] - (i ? mtx->end_of_row_offsets[i - 1] : 0));
        printf("\n");
    }

    printf("]\n");
}

void print_ccs_matrix(const jmtx_matrix_ccs* mtx)
{
    printf("elements:");
    for (uint32_t i = 0, l = 0; i < mtx->n_elements; ++i)
    {
        if (l < mtx->base.rows && mtx->elements_before[l + 1] <= i + 1)
        {
            l += 1;
        }
        printf(" (%u, %u, %g)", l, mtx->indices[i + 1], mtx->elements[i + 1]);
    }
    printf("\nelement offsets:");
    for (uint32_t i = 0; i < mtx->base.rows; ++i)
    {
        printf(" %u,", mtx->elements_before[i]);
    }
    printf(" %u", mtx->elements_before[mtx->base.rows]);
    printf("\nMatrix:\n[");

    for (uint32_t i = 0; i < mtx->base.rows; ++i)
    {
        if (i != 0)
            printf(" [");
        else
            printf("[");
        float x;
        for (uint32_t j = 0; j < mtx->base.cols - 1; ++j)
        {
            jmtx_matrix_ccs_get_element(mtx, i, j, &x);
            printf("%g ", x);
        }
        jmtx_matrix_ccs_get_element(mtx, i, mtx->base.cols - 1, &x);
        printf("%g] - %u", x, mtx->elements_before[i + 1] - mtx->elements_before[i]);
        if (i != mtx->base.rows - 1) printf("\n");
    }

    printf("]\n");
}
