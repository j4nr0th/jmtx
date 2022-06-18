//
// Created by jan on 16.6.2022.
//
#include "sparse_row_compressed.h"
#include "jacobi_point_iteration.h"
#include "common.h"
#include "gauss_seidel_iteration.h"
#include "bicgstab_iteration.h"
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <time.h>


#define MATRIX_TEST_DIMS 10000
#define DX 0.1f
#define VIS 0.3f

int main(int argc, char* argv[])
    {
        for (uint I = 1; I < argc; ++I)
        {
            char* end;
            const uint dims = strtoul(argv[I], &end, 10);
            if (end == argv[I])
            {
                printf("Could not translate \"%s\" into an integer\n", argv[I]);
                continue;
            }
//            const uint dims = 4096;
            CrsMatrix matrix;
            printf("Creating a %u by %u matrix\n", dims, dims);
            matrix_crs_new(&matrix, dims, dims, 3 * dims);
            printf("Setting the top BC\n");
            {
                const uint index[1] = {0};
                const scalar_t v[1] = {1};
                matrix_crs_build_row(&matrix, 0, 1, index, v);
            }
            printf("Setting the middle rows\n");
            const scalar_t values[3] = {-1/(2 * DX) + (1 - VIS) / (DX * DX), 2 * (VIS - 1) / (DX * DX) - 1, 1/(2*DX) + (1 - VIS)/(DX * DX)};
            uint indices[3] = {0, 1, 2};
            clock_t t0, t1, t2, t3;

            for (uint i = 1; i < dims - 1; ++i)
            {
                matrix_crs_build_row(&matrix, i, 3, indices, values);
                indices[0] += 1;
                indices[1] += 1;
                indices[2] += 1;
            }

            printf("Setting the bottom BC\n");
            {
                const uint index[] = {dims - 1};
                const scalar_t v[] = {1};
                matrix_crs_build_row(&matrix, dims - 1, 1, index, v);
            }


            printf("Allocating memory for x\n");
            scalar_t* const x = calloc(dims * 2, sizeof*x);
            if (!x)
            {
                perror("Could not allocate memory for vector x");
                exit(EXIT_FAILURE);
            }
            scalar_t* const x2 = x + dims;

            printf("Allocating memory for y\n");
            scalar_t* const y = calloc(dims, sizeof*y);
            if (!y)
            {
                perror("Could not allocate memory for vector y");
                exit(EXIT_FAILURE);
            }
            y[0] = 1;
            y[dims - 1] = 1;

            printf("Computing y\n");
            for (uint i = 0; i < dims; ++i)
            {
                y[i] = cosf((scalar_t)i / (scalar_t)(dims - 1) * 2 * (scalar_t)M_PI);
            }

            uint n_iter;
            scalar_t error;
            printf("Solving for a problem of size %u\n", dims);
            t0 = clock();
            bicgstab_crs(&matrix, y, x, 1e-6f, 1000, &n_iter, &error);
            t1 = clock();

            printf("Solution obtained with error %f after %u iterations (%g ms)\n", error, n_iter, ((double)(t1 - t0) / (double)CLOCKS_PER_SEC) * 1e3);

            t0 = clock();
            bicgstab_crs_mt(&matrix, y, x2, 1e-6f, 1000, &n_iter, &error, 8);
            t1 = clock();

            printf("Solution obtained with error %f after %u iterations (%g ms)\n", error, n_iter, ((double)(t1 - t0) / (double)CLOCKS_PER_SEC) * 1e3);

            scalar_t err_s = 0, err = 0;
            for (uint i = 0; i < dims; ++i)
            {
                const scalar_t d = x[i] - x2[i];
                err += fabsf(d);
                err_s += d * d;
            }
            printf("Errors:\n\tSum of absolute: %g\n\tMagnitude: %g\n\tRMS: %g\n", err, sqrtf(err_s), sqrtf(err_s / (scalar_t)dims));

//            printf("\n[");
//            for (uint i = 0; i < dims; ++i)
//            {
//                printf(" %g,", x[i]);
//            }
//            printf(" ]\n");


            printf("Freeing the matrix\n");
            matrix_crs_destroy(&matrix);
            printf("Freeing x\n");
            free(x);
            printf("Freeing y\n");
            free(y);
        }

    return 0;
}
