//
// Created by jan on 11.7.2023.
//
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../test_common.h"
#include "../../../include/jmtx/float/matrices/sparse_column_compressed.h"
#include "../../../include/jmtx/float/matrices/sparse_column_compressed_safe.h"

int main()
{
    jmtx_matrix_ccs* mtx;
    const uint32_t n_rows = 16, n_cols = 16;
    jmtx_result mtx_res;
    MATRIX_TEST_CALL(jmtxs_matrix_ccs_new(&mtx, n_rows, n_cols, 4, NULL));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);

    //  Empty matrix should have only zeros everywhere
    for (uint32_t row = 0; row < n_rows; ++row)
    {
        for (uint32_t col = 0; col < n_cols; ++col)
        {
            float v;
            mtx_res = jmtxs_matrix_ccs_get_entry(mtx, row, col, &v);
            ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
            ASSERT(v == 0.0f);
//            printf("%7g ", (double)v);
        }
//        printf("\n");
    }

    float v;
    //  Attempt to access any value outside of  matrix should return and error code
    MATRIX_TEST_CALL(jmtxs_matrix_ccs_get_entry(mtx, n_rows, n_cols - 1, &v));
    ASSERT(mtx_res == JMTX_RESULT_INDEX_OUT_OF_BOUNDS);
    MATRIX_TEST_CALL(jmtxs_matrix_ccs_get_entry(mtx, n_rows - 1, n_cols, &v));
    ASSERT(mtx_res == JMTX_RESULT_INDEX_OUT_OF_BOUNDS);
    MATRIX_TEST_CALL(jmtxs_matrix_ccs_get_entry(mtx, n_rows, n_cols, &v));
    ASSERT(mtx_res == JMTX_RESULT_INDEX_OUT_OF_BOUNDS);
    MATRIX_TEST_CALL(jmtxs_matrix_ccs_get_entry(mtx, 1241284124, 1212412, &v));
    ASSERT(mtx_res == JMTX_RESULT_INDEX_OUT_OF_BOUNDS);
    MATRIX_TEST_CALL(jmtxs_matrix_ccs_get_entry(mtx, 0, -1, &v));
    ASSERT(mtx_res == JMTX_RESULT_INDEX_OUT_OF_BOUNDS);
    MATRIX_TEST_CALL(jmtxs_matrix_ccs_set_entry(mtx, n_rows, n_cols - 1, 6.9f));
    ASSERT(mtx_res == JMTX_RESULT_INDEX_OUT_OF_BOUNDS);
    MATRIX_TEST_CALL(jmtxs_matrix_ccs_set_entry(mtx, n_rows - 1, n_cols, 420.0f));
    ASSERT(mtx_res == JMTX_RESULT_INDEX_OUT_OF_BOUNDS);
    MATRIX_TEST_CALL(jmtxs_matrix_ccs_set_entry(mtx, n_rows, n_cols, 420.0f));
    ASSERT(mtx_res == JMTX_RESULT_INDEX_OUT_OF_BOUNDS);
    MATRIX_TEST_CALL(jmtxs_matrix_ccs_set_entry(mtx, 1241284124, 1212412, 420.0f));
    ASSERT(mtx_res == JMTX_RESULT_INDEX_OUT_OF_BOUNDS);
    MATRIX_TEST_CALL(jmtxs_matrix_ccs_set_entry(mtx, 0, -1, 420.0f));
    ASSERT(mtx_res == JMTX_RESULT_INDEX_OUT_OF_BOUNDS);

    //  Setting a whole row
    float some_values[15] =
            {
            1.23f, 123.0f, 4830.0f, 21.43f, 32414.92f,
            940.0f, 924.34f, 450.0f, 342.03f, 5.30f,
            -344.0f, 34.0f, 240.0f, 32.1f, 324.0f,
            };
    uint32_t indices[15] =
            {
            0, 1, 2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
            };
    MATRIX_TEST_CALL(jmtxs_matrix_ccs_set_col(mtx, 3, 15, indices, some_values));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);

    //  fourth row should have correct values now
    for (uint32_t row = 0; row < n_rows; ++row)
    {
        for (uint32_t col = 0; col < n_cols; ++col)
        {
            mtx_res = jmtxs_matrix_ccs_get_entry(mtx, row, col, &v);
            ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
            if (col != 3)
            {
                ASSERT(v == 0.0f);
            }
            else
            {
                uint32_t i;
                for (i = 0; i < 15; ++i)
                {
                    if (indices[i] == row) break;
                }
                if (i != 15)
                {
                    ASSERT(v == some_values[i]);
                }
                else
                {
                    ASSERT(v == 0.0f);
                }
            }
        }
    }


    //  Set the diagonal to negative zeros
    for (uint32_t i = 0; i < n_rows && i < n_cols; ++i)
    {
        mtx_res = jmtxs_matrix_ccs_set_entry(mtx, i, i, -0.0f);
        ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    }

    MATRIX_TEST_CALL(jmtxs_matrix_ccs_remove_zeros(mtx));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);

    //  Matrix should not contain any negative zeros
    for (uint32_t row = 0; row < n_rows; ++row)
    {
        for (uint32_t col = 0; col < n_cols; ++col)
        {
            mtx_res = jmtxs_matrix_ccs_get_entry(mtx, row, col, &v);
            ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
            //  No negative zeros!
            ASSERT(v != 0.0f || !signbit(v));
        }
    }

    const float other_values[16] =
            {
            124.756f, -123e5f, -1.35f, 31.092f,
            490.0231e-4f, 123.21f, -232e9f, +195.2f,
            -412.0f, 5556.0f, -513.04f, 4494.342f,
            95.3053f, 441.034f, 596.3f, 3059.03f,
            };
    const uint32_t row_indices[16] =
            {
             1,  5,  2, 12,
             3,  4,  6,  5,
             7, 15, 15,  2,
            11, 12,  3,  4,
            };
    const uint32_t col_indices[16] =
            {
             0,  0,  2,  5,
            14, 15, 13, 13,
             2,  4, 15,  8,
             7,  7,  2,  6,
            };

    for (uint32_t i = 0; i < 16; ++i)
    {
//        print_ccs_matrix(&mtx);
        mtx_res = jmtxs_matrix_ccs_set_entry(mtx, row_indices[i], col_indices[i], other_values[i]);
        ASSERT(mtx_res == JMTX_RESULT_SUCCESS);

//        print_ccs_matrix(&mtx);
        for (uint32_t j = 0; j <= i; ++j)
        {
            mtx_res = jmtxs_matrix_ccs_get_entry(mtx, row_indices[j], col_indices[j], &v);
            ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
            ASSERT(v == other_values[j]);
        }
    }
//    print_ccs_matrix(&mtx);

    for (uint32_t i = 0; i < 16; ++i)
    {
        mtx_res = jmtxs_matrix_ccs_get_entry(mtx, row_indices[i], col_indices[i], &v);
        ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
        ASSERT(v == other_values[i]);
    }

    //  Some transpose action
    for (uint32_t i = 0; i < 16; ++i)
    {
        mtx_res = jmtxs_matrix_ccs_set_entry(mtx, col_indices[i], row_indices[i], other_values[i]);
        ASSERT(mtx_res == JMTX_RESULT_SUCCESS);

        for (uint32_t j = 0; j <= i; ++j)
        {
            mtx_res = jmtxs_matrix_ccs_get_entry(mtx, col_indices[j], row_indices[j], &v);
            ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
            ASSERT(v == other_values[j]);
        }
    }
//    print_ccs_matrix(&mtx);

    for (uint32_t i = 0; i < 16; ++i)
    {
        mtx_res = jmtxs_matrix_ccs_get_entry(mtx, col_indices[i], row_indices[i], &v);
        ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
        ASSERT(v == other_values[i]);
    }


//    print_ccs_matrix(&mtx);


    MATRIX_TEST_CALL(jmtxs_matrix_ccs_destroy(mtx));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    return 0;
}