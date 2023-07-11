//
// Created by jan on 15.6.2022.
//
#include <stdio.h>
#include "source/matrices/sparse_column_compressed.h"
#include "source/matrices/sparse_row_compressed.h"
#include "tests/test_common.h"
#include <math.h>

#define TEST_MATRIX_MAJOR_DIM 255
#define TEST_MATRIX_MINOR_DIM 50


int main()
{
    fRNG rng = fRNG_create(0xBADBABEB00B1E5);
    scalar_t X[TEST_MATRIX_MINOR_DIM];
    scalar_t Y_crs[TEST_MATRIX_MAJOR_DIM], Y_ccs[TEST_MATRIX_MAJOR_DIM];
    CrsMatrix crs_matrix_1;
    CrsMatrix crs_matrix_2;
    matrix_crs_new(&crs_matrix_1, TEST_MATRIX_MINOR_DIM, TEST_MATRIX_MAJOR_DIM, 0);
    CcsMatrix ccs_matrix_1;
//    CcsMatrix ccs_matrix_2;
    matrix_ccs_new(&ccs_matrix_1, TEST_MATRIX_MAJOR_DIM, TEST_MATRIX_MINOR_DIM, 0);

    for (uint32_t i = 0; i < TEST_MATRIX_MAJOR_DIM; ++i)
    {
        X[i] = fRNG_float_range(&rng, -10, 10);
        scalar_t val[TEST_MATRIX_MINOR_DIM];
        uint32_t     ind[TEST_MATRIX_MINOR_DIM];
        const uint32_t n_elements = (uint32_t) TEST_MATRIX_MINOR_DIM;
        for (uint32_t j = 0; j < n_elements; ++j)
        {
            val[j] = fRNG_float_range(&rng, -10, 10);
            ind[j] = j + (TEST_MATRIX_MINOR_DIM - n_elements) / 2;
        }
        matrix_crs_set_row(&crs_matrix_1, i, n_elements, ind, val);
        matrix_ccs_set_col(&ccs_matrix_1, i, n_elements, ind, val);
    }

    print_crs_matrix(&crs_matrix_1);
    print_ccs_matrix(&ccs_matrix_1);

    matrix_crs_shrink(&crs_matrix_1);
    matrix_ccs_shrink(&ccs_matrix_1);

    matrix_crs_vector_multiply(&crs_matrix_1, X, Y_crs);
    matrix_ccs_vector_multiply(&ccs_matrix_1, X, Y_ccs);

    scalar_t residual = 0.0f;
    for (uint32_t i = 0; i < TEST_MATRIX_MAJOR_DIM; ++i)
    {
        residual += fabsf(Y_ccs[i] - Y_crs[i]);
    }
    printf("Residual of vector multiplication: %g\n", residual);


    matrix_crs_transpose(&crs_matrix_1, &crs_matrix_2);
    residual = 0.0f;
    for (uint32_t j = 0; j < TEST_MATRIX_MAJOR_DIM; ++j)
    {
        for (uint32_t i = 0; i < TEST_MATRIX_MINOR_DIM; ++i)
        {
            scalar_t x1, x2;
            matrix_crs_get_element(&crs_matrix_2, i, j, &x1);
            matrix_ccs_get_element(&ccs_matrix_1, i, j, &x2);
            residual += fabsf(x1 - x2);
        }
    }

    printf("Residual of crs transpose: %g\n", residual);


    matrix_crs_destroy(&crs_matrix_2);
    matrix_crs_destroy(&crs_matrix_1);
    matrix_ccs_destroy(&ccs_matrix_1);
    return 0;
}


