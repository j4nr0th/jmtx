// Automatically generated from tests/float/matrix_ops_tests/crs_addition_test.c on Fri Dec  1 17:35:45 2023
//
// Created by jan on 21.7.2023.
//
#include <stdio.h>
#include "../../../include/jmtx/cfloat/matrices/sparse_row_compressed.h"
#include "../test_common.h"
#include "../../../include/jmtx/cfloat/matrices/sparse_row_compressed_safe.h"

int main()
{
    jmtxc_matrix_crs* mtx;
    jmtx_result mtx_res;
    MATRIX_TEST_CALL(jmtxcs_matrix_crs_new(&mtx, 4, 4, 4, NULL));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);

    const _Complex float m1[16] =
            {
            1, 1, 2, -3,
            0, 2, -4, 2,
            0, 0, 3, 0,
            0, 0, 0, 4,
            };

    const _Complex float m2[16] =
            {
            -3, 2, 0, 0,
            0, 1, 0, 0,
            0, 0, -3, 2,
            0, 0, 0, 1,
            };

    const _Complex float m3[16] =
            {
            1, 0, 0, 1,
            0, 2, 3, 0,
            0, 3, 2, 0,
            1, 0, 0, 1
            };
    int beef_status;

    //  Matrix is empty, so this is effectively just setting the entries
    for (uint32_t i = 0; i < 4; ++i)
    {
        for (uint32_t j = 0; j < 4; ++j)
        {
            _Complex float v = m1[j + i * 4];
            if (v == 0)
            {
                continue;
            }
            mtx_res = jmtxcs_matrix_crs_add_to_entry(mtx, i, j, v);
            ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
        }
    }
    printf("Matrix after adding m1:\n");
    print_crsc_matrix(mtx);
    printf("\n");

    for (uint32_t i = 0; i < 4; ++i)
    {
        for (uint32_t j = 0; j < 4; ++j)
        {
            _Complex float x = m1[j + i * 4];
            _Complex float v;
            mtx_res = jmtxcs_matrix_crs_get_entry(mtx, i, j, &v);
            ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
            ASSERT(v == x);
        }
    }

    //  Matrix is no longer empty, so this is adding m2 to m1
    for (uint32_t i = 0; i < 4; ++i)
    {
        for (uint32_t j = 0; j < 4; ++j)
        {
            _Complex float v = m2[j + i * 4];
            if (v == 0)
            {
                continue;
            }
            mtx_res = jmtxcs_matrix_crs_add_to_entry(mtx, i, j, v);
            ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
        }
    }
    printf("Matrix after adding m2:\n");
    print_crsc_matrix(mtx);
    printf("\n");

    for (uint32_t i = 0; i < 4; ++i)
    {
        for (uint32_t j = 0; j < 4; ++j)
        {
            _Complex float x = m1[j + i * 4] + m2[j + i * 4];
            _Complex float v;
            mtx_res = jmtxcs_matrix_crs_get_entry(mtx, i, j, &v);
            ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
            ASSERT(v == x);
        }
    }

    //  Matrix is no longer empty, so this is adding m1, m2, and m3 together
    for (uint32_t i = 0; i < 4; ++i)
    {
        for (uint32_t j = 0; j < 4; ++j)
        {
            _Complex float v = m3[j + i * 4];
            if (v == 0)
            {
                continue;
            }
            mtx_res = jmtxcs_matrix_crs_add_to_entry(mtx, i, j, v);
            ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
        }
    }
    printf("Matrix after adding m3:\n");
    print_crsc_matrix(mtx);
    printf("\n");

    for (uint32_t i = 0; i < 4; ++i)
    {
        for (uint32_t j = 0; j < 4; ++j)
        {
            _Complex float x = m1[j + i * 4] + m2[j + i * 4] + m3[j + i * 4];
            _Complex float v;
            mtx_res = jmtxcs_matrix_crs_get_entry(mtx, i, j, &v);
            ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
            ASSERT(v == x);
        }
    }

    //  Removing 2nd column
    MATRIX_TEST_CALL(jmtxcs_matrix_crs_remove_column(mtx, 1));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);

    printf("Matrix after removing 2nd column:\n");
    print_crsc_matrix(mtx);
    printf("\n");

    for (uint32_t i = 0; i < 4; ++i)
    {
        for (uint32_t j = 0; j < 3; ++j)
        {
            uint32_t c = j > 0 ? j + 1 : 0;
            _Complex float x = m1[c + i * 4] + m2[c + i * 4] + m3[c + i * 4];
            _Complex float v;
            mtx_res = jmtxcs_matrix_crs_get_entry(mtx, i, j, &v);
            ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
            ASSERT(v == x);
        }
    }

    //  Removing 1st row
    MATRIX_TEST_CALL(jmtxcs_matrix_crs_remove_row(mtx, 0));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);

    printf("Matrix after removing 1st row:\n");
    print_crsc_matrix(mtx);
    printf("\n");

    for (uint32_t i = 1; i < 4; ++i)
    {
        for (uint32_t j = 0; j < 3; ++j)
        {
            uint32_t c = j > 0 ? j + 1 : 0;
            _Complex float x = m1[c + i * 4] + m2[c + i * 4] + m3[c + i * 4];
            _Complex float v;
            mtx_res = jmtxcs_matrix_crs_get_entry(mtx, i - 1, j, &v);
            ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
            ASSERT(v == x);
        }
    }

    MATRIX_TEST_CALL(jmtxcs_matrix_crs_destroy(mtx));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    return 0;
}
