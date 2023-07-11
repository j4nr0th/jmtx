//
// Created by jan on 11.7.2023.
//
#include <stdio.h>
#include <stdlib.h>
#include "test_common.h"
#include "../source/matrices/sparse_row_compressed.h"

#ifndef NDEBUG
#   ifdef __GNUC__
#       define DBG_BREAK __builtin_trap()
#   endif
#else
#define DBG_BREAK (void)0
#endif

#define ASSERT(x) if (!(x)) {fprintf(stderr, "Failed assertion \"" #x "\" on line %u, file %s\n", __LINE__, __FILE__); DBG_BREAK; exit(EXIT_FAILURE);} (void)0
#define MATRIX_TEST_CALL(x) printf("Called:\t"#x" -> %s\n", jmtx_result_to_str((mtx_res = (x))))
int main()
{
    jmtx_matrix_crs mtx;
    const uint32_t n_rows = 16, n_cols = 16;
    jmtx_result mtx_res;
    MATRIX_TEST_CALL(matrix_crs_new(&mtx, n_rows, n_cols, 4, NULL));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);

    //  Empty matrix should have only zeros everywhere
    for (uint32_t row = 0; row < n_rows; ++row)
    {
        for (uint32_t col = 0; col < n_cols; ++col)
        {
            jmtx_scalar_t v;
            mtx_res = matrix_crs_get_element(&mtx, row, col, &v);
            ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
            ASSERT(v == 0.0f);
//            printf("%7g ", (double)v);
        }
//        printf("\n");
    }

    jmtx_scalar_t v;
    //  Attempt to access any value outside of  matrix should return and error code
    MATRIX_TEST_CALL(matrix_crs_get_element(&mtx, n_rows, n_cols - 1, &v));
    ASSERT(mtx_res == JMTX_RESULT_INDEX_OUT_OF_BOUNDS);
    MATRIX_TEST_CALL(matrix_crs_get_element(&mtx, n_rows - 1, n_cols, &v));
    ASSERT(mtx_res == JMTX_RESULT_INDEX_OUT_OF_BOUNDS);
    MATRIX_TEST_CALL(matrix_crs_get_element(&mtx, n_rows, n_cols, &v));
    ASSERT(mtx_res == JMTX_RESULT_INDEX_OUT_OF_BOUNDS);
    MATRIX_TEST_CALL(matrix_crs_get_element(&mtx, 1241284124, 1212412, &v));
    ASSERT(mtx_res == JMTX_RESULT_INDEX_OUT_OF_BOUNDS);
    MATRIX_TEST_CALL(matrix_crs_get_element(&mtx, 0, -1, &v));
    ASSERT(mtx_res == JMTX_RESULT_INDEX_OUT_OF_BOUNDS);
    MATRIX_TEST_CALL(matrix_crs_set_element(&mtx, n_rows, n_cols - 1, 6.9f));
    ASSERT(mtx_res == JMTX_RESULT_INDEX_OUT_OF_BOUNDS);
    MATRIX_TEST_CALL(matrix_crs_set_element(&mtx, n_rows - 1, n_cols, 420.0f));
    ASSERT(mtx_res == JMTX_RESULT_INDEX_OUT_OF_BOUNDS);
    MATRIX_TEST_CALL(matrix_crs_set_element(&mtx, n_rows, n_cols, 420.0f));
    ASSERT(mtx_res == JMTX_RESULT_INDEX_OUT_OF_BOUNDS);
    MATRIX_TEST_CALL(matrix_crs_set_element(&mtx, 1241284124, 1212412, 420.0f));
    ASSERT(mtx_res == JMTX_RESULT_INDEX_OUT_OF_BOUNDS);
    MATRIX_TEST_CALL(matrix_crs_set_element(&mtx, 0, -1, 420.0f));
    ASSERT(mtx_res == JMTX_RESULT_INDEX_OUT_OF_BOUNDS);

    //  Setting a whole row
    jmtx_scalar_t some_values[15] =
            {
            1.23f, 123.0f, 4830.0f, 21.43f, 32414.92f,
            940.0f, 924.34f, 450.0f, 342.03f, 5.30f,
            -344.0f, 34.0f, 240.0f, 32.1f, 324.0f,
            };
    uint32_t indices[15] =
            {
            0, 1, 2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
            };
    MATRIX_TEST_CALL(matrix_crs_set_row(&mtx, 3, 15, indices, some_values));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);

    //  fourth row should have correct values now
    for (uint32_t row = 0; row < n_rows; ++row)
    {
        for (uint32_t col = 0; col < n_cols; ++col)
        {
            mtx_res = matrix_crs_get_element(&mtx, row, col, &v);
            ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
            if (row != 3)
            {
                ASSERT(v == 0.0f);
            }
            else
            {
                uint32_t i;
                for (i = 0; i < 15; ++i)
                {
                    if (indices[i] == col) break;
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
        mtx_res = matrix_crs_set_element(&mtx, i, i, -0.0f);
        ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    }



    print_crs_matrix(&mtx);


    MATRIX_TEST_CALL(matrix_crs_destroy(&mtx));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    return 0;
}