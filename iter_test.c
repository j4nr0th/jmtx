//
// Created by jan on 16.6.2022.
//
#include "sparse_row_compressed.h"
#include "jacobi_point_iteration.h"
#include "common.h"
#include "gauss_seidel_iteration.h"
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <time.h>


#define MATRIX_TEST_DIMS 10000
#define DX 0.1f
#define VIS 0.3f

int main()
{
    CrsMatrix matrix;
    matrix_crs_new(&matrix, MATRIX_TEST_DIMS, MATRIX_TEST_DIMS, 0);

    {
        const uint index[1] = {0};
        const scalar_t v[1] = {1};
        matrix_crs_set_row(&matrix, 0, 1, index, v);
    }

    const scalar_t values[3] = {-1/(2 * DX) + (1 - VIS) / (DX * DX), 2 * (VIS - 1) / (DX * DX) - 1, 1/(2*DX) + (1 - VIS)/(DX * DX)};
    uint indices[3] = {0, 1, 2};
    for (uint i = 1; i < MATRIX_TEST_DIMS - 1; ++i)
    {
        matrix_crs_set_row(&matrix, i, 3, indices, values);
        indices[0] += 1;
        indices[1] += 1;
        indices[2] += 1;
    }

    {
        const uint index[] = {MATRIX_TEST_DIMS - 1};
        const scalar_t v[] = {1};
        matrix_crs_set_row(&matrix, MATRIX_TEST_DIMS - 1, 1, index, v);
    }

    scalar_t x[MATRIX_TEST_DIMS];
    scalar_t y[MATRIX_TEST_DIMS] = {[0] = 1, [MATRIX_TEST_DIMS - 1] = 1};

    for (uint i = 0; i < MATRIX_TEST_DIMS; ++i)
    {
        y[i] = cosf((scalar_t)i / (MATRIX_TEST_DIMS - 1) * 2 * (scalar_t)M_PI);
    }

    uint n_iter;
    clock_t t0 = clock();
    gauss_seidel_crs_mt(&matrix, y, x, 1e-5f, 1000000, &n_iter, 8);
    clock_t t1 = clock();

    printf("Solution obtained after %u iterations (%g ms)\n", n_iter, ((double)(t1 - t0) / (double)CLOCKS_PER_SEC) * 1e3);
//    print_crs_matrix(&matrix);
    printf("\n[");
    for (uint i = 0; i < MATRIX_TEST_DIMS; ++i)
    {
        printf(" %g,", x[i]);
    }
    printf(" ]\n");
    matrix_crs_destroy(&matrix);
    return 0;
}
