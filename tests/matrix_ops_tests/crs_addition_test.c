#include <stdio.h>
#include "../test_common.h"

int main()
{
    JMTX_NAME_TYPED(matrix_crs) * mtx;
    jmtx_result mtx_res;
    MATRIX_TEST_CALL(JMTX_NAME_TYPED(matrix_crs_new)(&mtx, 4, 4, 4, NULL));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);

    const JMTX_SCALAR_T m1[16] = {
        1, 1, 2, -3, 0, 2, -4, 2, 0, 0, 3, 0, 0, 0, 0, 4,
    };

    const JMTX_SCALAR_T m2[16] = {
        -3, 2, 0, 0, 0, 1, 0, 0, 0, 0, -3, 2, 0, 0, 0, 1,
    };

    const JMTX_SCALAR_T m3[16] = {1, 0, 0, 1, 0, 2, 3, 0, 0, 3, 2, 0, 1, 0, 0, 1};
    int beef_status;

    //  Matrix is empty, so this is effectively just setting the entries
    for (uint32_t i = 0; i < 4; ++i)
    {
        for (uint32_t j = 0; j < 4; ++j)
        {
            double v = m1[j + i * 4];
            if (v == 0)
            {
                continue;
            }
            mtx_res = JMTX_NAME_TYPED(matrix_crs_add_to_entry)(mtx, i, j, v);
            ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
        }
    }
    printf("Matrix after adding m1:\n");
    JMTX_NAME_TYPED(print_crs_matrix)(mtx);
    printf("\n");

    for (uint32_t i = 0; i < 4; ++i)
    {
        for (uint32_t j = 0; j < 4; ++j)
        {
            const JMTX_SCALAR_T x = m1[j + i * 4];
            const JMTX_SCALAR_T v = JMTX_NAME_TYPED(matrix_crs_get_entry)(mtx, i, j);
            ASSERT(v == x);
        }
    }

    //  Matrix is no longer empty, so this is adding m2 to m1
    for (uint32_t i = 0; i < 4; ++i)
    {
        for (uint32_t j = 0; j < 4; ++j)
        {
            const JMTX_SCALAR_T v = m2[j + i * 4];
            if (v == 0)
            {
                continue;
            }
            mtx_res = JMTX_NAME_TYPED(matrix_crs_add_to_entry)(mtx, i, j, v);
            ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
        }
    }
    printf("Matrix after adding m2:\n");
    JMTX_NAME_TYPED(print_crs_matrix)(mtx);
    printf("\n");

    for (uint32_t i = 0; i < 4; ++i)
    {
        for (uint32_t j = 0; j < 4; ++j)
        {
            const JMTX_SCALAR_T x = m1[j + i * 4] + m2[j + i * 4];
            const JMTX_SCALAR_T v = JMTX_NAME_TYPED(matrix_crs_get_entry)(mtx, i, j);
            ASSERT(v == x);
        }
    }

    //  Matrix is no longer empty, so this is adding m1, m2, and m3 together
    for (uint32_t i = 0; i < 4; ++i)
    {
        for (uint32_t j = 0; j < 4; ++j)
        {
            JMTX_SCALAR_T v = m3[j + i * 4];
            if (v == 0)
            {
                continue;
            }
            mtx_res = JMTX_NAME_TYPED(matrix_crs_add_to_entry)(mtx, i, j, v);
            ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
        }
    }
    printf("Matrix after adding m3:\n");
    JMTX_NAME_TYPED(print_crs_matrix)(mtx);
    printf("\n");

    for (uint32_t i = 0; i < 4; ++i)
    {
        for (uint32_t j = 0; j < 4; ++j)
        {
            const JMTX_SCALAR_T x = m1[j + i * 4] + m2[j + i * 4] + m3[j + i * 4];
            const JMTX_SCALAR_T v = JMTX_NAME_TYPED(matrix_crs_get_entry)(mtx, i, j);
            ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
            ASSERT(v == x);
        }
    }

    //  Removing 2nd column
    JMTX_NAME_TYPED(matrix_crs_remove_column)(mtx, 1);
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);

    printf("Matrix after removing 2nd column:\n");
    JMTX_NAME_TYPED(print_crs_matrix)(mtx);
    printf("\n");

    for (uint32_t i = 0; i < 4; ++i)
    {
        for (uint32_t j = 0; j < 3; ++j)
        {
            const uint32_t c = j > 0 ? j + 1 : 0;
            const JMTX_SCALAR_T x = m1[c + i * 4] + m2[c + i * 4] + m3[c + i * 4];
            const JMTX_SCALAR_T v = JMTX_NAME_TYPED(matrix_crs_get_entry)(mtx, i, j);
            ASSERT(v == x);
        }
    }

    //  Removing 1st row
    JMTX_NAME_TYPED(matrix_crs_remove_row)(mtx, 0);
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);

    printf("Matrix after removing 1st row:\n");
    JMTX_NAME_TYPED(print_crs_matrix)(mtx);
    printf("\n");

    for (uint32_t i = 1; i < 4; ++i)
    {
        for (uint32_t j = 0; j < 3; ++j)
        {
            const uint32_t c = j > 0 ? j + 1 : 0;
            const JMTX_SCALAR_T x = m1[c + i * 4] + m2[c + i * 4] + m3[c + i * 4];
            const JMTX_SCALAR_T v = JMTX_NAME_TYPED(matrix_crs_get_entry)(mtx, i - 1, j);
            ASSERT(v == x);
        }
    }

    JMTX_NAME_TYPED(matrix_crs_destroy)(mtx);
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    return 0;
}
