//
// Created by jan on 16.6.2022.
//
#include "source/matrices/sparse_row_compressed.h"
#include "source/solvers/jacobi_point_iteration.h"
#include "tests/test_common.h"
#include "source/solvers/gauss_seidel_iteration.h"
#include "source/solvers/bicgstab_iteration.h"
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <stdlib.h>

#define MATRIX_TEST_DIMS 10000
#define DX 0.1f
#define VIS 0.3f

int main(int argc, char* argv[])
{
    for (uint32_t I = 1; I < argc; ++I)
    {
        char* end;
        const uint32_t dims = strtoul(argv[I], &end, 10);
        if (end == argv[I])
        {
            printf("Could not translate \"%s\" into an integer\n", argv[I]);
            continue;
        }
//            const uint32_t dims = 4096;
        jmtx_matrix_crs matrix = { 0};
        printf("Creating a %u by %u matrix\n", dims, dims);
        matrix_crs_new(&matrix, dims, dims, 3 * dims, NULL);
        printf("Setting the top BC\n");
        {
            const uint32_t index[1] = {0};
            const float v[1] = { 1};
            matrix_crs_build_row(&matrix, 0, 1, index, v);
        }
        printf("Setting the middle rows\n");
        const float values[3] = { -1 / (2 * DX) + (1 - VIS) / (DX * DX), 2 * (VIS - 1) / (DX * DX) - 1, 1 / (2 * DX) + (1 - VIS) / (DX * DX)};
        uint32_t indices[3] = {0, 1, 2};
//            clock_t t0, t1, t2, t3;

        for (uint32_t i = 1; i < dims - 1; ++i)
        {
            matrix_crs_build_row(&matrix, i, 3, indices, values);
            indices[0] += 1;
            indices[1] += 1;
            indices[2] += 1;
        }

        printf("Setting the bottom BC\n");
        {
            const uint32_t index[] = {dims - 1};
            const float v[] = { 1};
            matrix_crs_build_row(&matrix, dims - 1, 1, index, v);
        }


        printf("Allocating memory for x\n");
        float* const x = calloc(dims * 2, sizeof*x);
        if (!x)
        {
            perror("Could not allocate memory for vector x");
            exit(EXIT_FAILURE);
        }
        float* const x2 = x + dims;

        printf("Allocating memory for y\n");
        float* const y = calloc(dims, sizeof*y);
        if (!y)
        {
            perror("Could not allocate memory for vector y");
            exit(EXIT_FAILURE);
        }
        y[0] = 1;
        y[dims - 1] = 1;

        printf("Computing y\n");
        for (uint32_t i = 0; i < dims; ++i)
        {
            y[i] = cosf((float)i / (float)(dims - 1) * 2 * (float)M_PI);
        }

        uint32_t n_iter;
        float error;
        printf("Solving for a problem of size %u\n", dims);
        struct timespec t0, t1, t2, t3;

        clock_gettime(CLOCK_MONOTONIC, &t0);
        bicgstab_crs(&matrix, y, x, 1e-5f, 1000, &n_iter, &error, NULL, NULL);
        clock_gettime(CLOCK_MONOTONIC, &t1);
        matrix_crs_vector_multiply(&matrix, x, x2);
        float residual_magnitude = 0;
        for (uint32_t i = 0; i < dims; ++i)
        {
            float residual = y[i] - x2[i];
            residual_magnitude += residual * residual;
        }
        residual_magnitude = sqrtf(residual_magnitude);

        printf("Solution obtained with error %f after %u iterations (%g s)\n", error, n_iter, ((double)(t1.tv_sec - t0.tv_sec) + (double)(t1.tv_nsec - t0.tv_nsec) * 1e-9) );
        printf("Total residual magnitude %g\nAverage residual magnitude %g\n", (residual_magnitude), (residual_magnitude) / (float)dims);

        clock_gettime(CLOCK_MONOTONIC, &t2);
        bicgstab_crs_mt(&matrix, y, x2, 1e-5f, 1000, &n_iter, &error, NULL, NULL, 8);
        clock_gettime(CLOCK_MONOTONIC, &t3);
        matrix_crs_vector_multiply(&matrix, x2, x);
        residual_magnitude = 0;
        for (uint32_t i = 0; i < dims; ++i)
        {
            float residual = y[i] - x[i];
            residual_magnitude += residual * residual;
        }
        residual_magnitude = sqrtf(residual_magnitude);

        printf("Solution obtained with error %f after %u iterations (%g s)\n", error, n_iter, ((double)(t3.tv_sec - t2.tv_sec) + (double)(t3.tv_nsec - t2.tv_nsec) * 1e-9));
        printf("Total residual magnitude %g\nAverage residual magnitude %g\n", (residual_magnitude), (residual_magnitude) / (float)dims);



        printf("Freeing the matrix\n");
        matrix_crs_destroy(&matrix);
        printf("Freeing x\n");
        free(x);
        printf("Freeing y\n");
        free(y);
    }
//    THREAD_END;
    return 0;
}
