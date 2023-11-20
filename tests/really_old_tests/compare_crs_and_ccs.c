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
    float X[TEST_MATRIX_MAJOR_DIM];
    float Y_crs[TEST_MATRIX_MAJOR_DIM], Y_ccs[TEST_MATRIX_MAJOR_DIM];
    jmtx_matrix_crs* crs_matrix_1;
    jmtx_matrix_crs* crs_matrix_2;
    jmtx_matrix_crs_new(&crs_matrix_1, TEST_MATRIX_MINOR_DIM, TEST_MATRIX_MAJOR_DIM, 0, NULL);
    jmtx_matrix_ccs* ccs_matrix_1;
//    jmtx_matrix_ccs ccs_matrix_2;
    jmtx_matrix_ccs_new(&ccs_matrix_1, TEST_MATRIX_MAJOR_DIM, TEST_MATRIX_MINOR_DIM, 0, NULL);

    for (uint32_t i = 0; i < TEST_MATRIX_MAJOR_DIM; ++i)
    {
        X[i] = fRNG_float_range(&rng, -10, 10);
        float val[TEST_MATRIX_MINOR_DIM];
        uint32_t     ind[TEST_MATRIX_MINOR_DIM];
        const uint32_t n_elements = (uint32_t) TEST_MATRIX_MINOR_DIM;
        for (uint32_t j = 0; j < n_elements; ++j)
        {
            val[j] = fRNG_float_range(&rng, -10, 10);
            ind[j] = j + (TEST_MATRIX_MINOR_DIM - n_elements) / 2;
        }
        jmtx_matrix_crs_set_row(crs_matrix_1, i, n_elements, ind, val);
        jmtx_matrix_ccs_set_col(ccs_matrix_1, i, n_elements, ind, val);
    }

//    print_crs_matrix(crs_matrix_1);
//    print_ccs_matrix(&ccs_matrix_1);

    jmtx_matrix_crs_shrink(crs_matrix_1);
    jmtx_matrix_ccs_shrink(ccs_matrix_1);

    jmtx_matrix_crs_vector_multiply(crs_matrix_1, X, Y_crs);
    jmtx_matrix_ccs_vector_multiply(ccs_matrix_1, X, Y_ccs);

    float residual = 0.0f;
    for (uint32_t i = 0; i < TEST_MATRIX_MAJOR_DIM; ++i)
    {
        residual += fabsf(Y_ccs[i] - Y_crs[i]);
    }
    printf("Residual of vector multiplication: %g\n", residual);


    jmtx_matrix_crs_transpose(crs_matrix_1, &crs_matrix_2, NULL);
    residual = 0.0f;
    for (uint32_t j = 0; j < TEST_MATRIX_MAJOR_DIM; ++j)
    {
        for (uint32_t i = 0; i < TEST_MATRIX_MINOR_DIM; ++i)
        {
            float x1, x2;
            x1 = jmtx_matrix_crs_get_entry(crs_matrix_2, i, j);
            x2 = jmtx_matrix_ccs_get_entry(ccs_matrix_1, i, j);
            residual += fabsf(x1 - x2);
        }
    }

    printf("Residual of crs transpose: %g\n", residual);


    jmtx_matrix_crs_destroy(crs_matrix_2);
    jmtx_matrix_crs_destroy(crs_matrix_1);
    jmtx_matrix_ccs_destroy(ccs_matrix_1);
    return 0;
}


