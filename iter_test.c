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
            matrix_crs_new(&matrix, dims, dims, 0);

            {
                const uint index[1] = {0};
                const scalar_t v[1] = {1};
                matrix_crs_set_row(&matrix, 0, 1, index, v);
            }

            const scalar_t values[3] = {-1/(2 * DX) + (1 - VIS) / (DX * DX), 2 * (VIS - 1) / (DX * DX) - 1, 1/(2*DX) + (1 - VIS)/(DX * DX)};
            uint indices[3] = {0, 1, 2};
            for (uint i = 1; i < dims - 1; ++i)
            {
                matrix_crs_set_row(&matrix, i, 3, indices, values);
                indices[0] += 1;
                indices[1] += 1;
                indices[2] += 1;
            }

            {
                const uint index[] = {dims - 1};
                const scalar_t v[] = {1};
                matrix_crs_set_row(&matrix, dims - 1, 1, index, v);
            }


            scalar_t* const x = calloc(dims, sizeof*x);
            if (!x)
            {
                perror("Could not allocate memory for vector x");
                exit(EXIT_FAILURE);
            }
            scalar_t* const y = calloc(dims, sizeof*y);
            if (!y)
            {
                perror("Could not allocate memory for vector y");
                exit(EXIT_FAILURE);
            }
            y[0] = 1;
            y[dims - 1] = 1;

            for (uint i = 0; i < dims; ++i)
            {
                y[i] = cosf((scalar_t)i / (scalar_t)(dims - 1) * 2 * (scalar_t)M_PI);
            }

            uint n_iter;
            printf("Solving for a problem of size %u\n", dims);
            clock_t t0 = clock();
            gauss_seidel_crs_mt(&matrix, y, x, 1e-5f, 20000, &n_iter, 8);
            clock_t t1 = clock();

            printf("Solution obtained after %u iterations (%g ms)\n", n_iter, ((double)(t1 - t0) / (double)CLOCKS_PER_SEC) * 1e3);
            //    print_crs_matrix(&matrix);
//            printf("\n[");
//            for (uint i = 0; i < dims; ++i)
//            {
//                printf(" %g,", x[i]);
//            }
//            printf(" ]\n");
            matrix_crs_destroy(&matrix);
            free(x);
            free(y);
        }

    return 0;
}
