//
// Created by jan on 24.11.2023.
//

#include "../test_common.h"
#include "../../include/jmtx/matrices/band_row_major_safe.h"

static int are_arrays_the_same(unsigned len, const float a1[const static len], const float a2[const static len])
{
    for (unsigned i = 0; i < len; ++i)
    {
        if (a1[i] != a2[i])
        {
            return 0;
        }
    }
    return 1;
}

int main()
{
    jmtx_matrix_brm* test_matrix;
    jmtx_result mtx_res;
    MATRIX_TEST_CALL(jmtxs_matrix_brm_new(&test_matrix, 5, 5, 2, 1, NULL, NULL));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);

    //  Set the first row
    {
        MATRIX_TEST_CALL(jmtxs_matrix_brm_set_row(test_matrix, 0, (float[]){1.0f, 2.0f, 3.0f}));
        ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    }
    //  Set the second row
    {
        MATRIX_TEST_CALL(jmtxs_matrix_brm_set_row(test_matrix, 1, (float[]){-1.0f, -2.0f, -3.0f, -4.0f}));
        ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    }
    //  Set the third row
    {
        MATRIX_TEST_CALL(jmtxs_matrix_brm_set_row(test_matrix, 2, (float[]){1.0f, 2.0f, 3.0f, 4.0f}));
        ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    }
    //  Set the fourth row
    {
        MATRIX_TEST_CALL(jmtxs_matrix_brm_set_row(test_matrix, 3, (float[]){-5.0f, -6.0f, -7.0f}));
        ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    }
    //  Set the fifth row
    {
        MATRIX_TEST_CALL(jmtxs_matrix_brm_set_row(test_matrix, 4, (float[]){69.0f, 420.0f}));
        ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    }
    print_brm_matrix(test_matrix);

    float column_vals[4];
    uint_fast32_t cnt;

    ASSERT((cnt = jmtx_matrix_brm_get_col(test_matrix, 0, column_vals)) == 2);
    ASSERT(are_arrays_the_same(cnt, column_vals, (float[]){1.0f, -1.0f}));

    ASSERT((cnt = jmtx_matrix_brm_get_col(test_matrix, 1, column_vals)) == 3);
    ASSERT(are_arrays_the_same(cnt, column_vals, (float[]){2.0f, -2.0f, 1.0f}));

    ASSERT((cnt = jmtx_matrix_brm_get_col(test_matrix, 2, column_vals)) == 4);
    ASSERT(are_arrays_the_same(cnt, column_vals, (float[]){3.0f, -3.0f, 2.0f, -5.0f}));

    ASSERT((cnt = jmtx_matrix_brm_get_col(test_matrix, 3, column_vals)) == 4);
    ASSERT(are_arrays_the_same(cnt, column_vals, (float[]){-4.0f, 3.0f, -6.0f, 69.0f}));

    ASSERT((cnt = jmtx_matrix_brm_get_col(test_matrix, 4, column_vals)) == 3);
    ASSERT(are_arrays_the_same(cnt, column_vals, (float[]){4.0f, -7.0f, 420.0f}));

    const float test_vec[5] = {1.0f, -1.0f, 3.0f, 0.0f, -3.0f};
    float out_vec[5];

    jmtx_matrix_brm_vector_multiply(test_matrix, test_vec, out_vec);
    for (unsigned i = 0; i < 5; ++i)
    {
        printf("%g ", out_vec[i]);
    }
    printf("\n");

    jmtx_matrix_brm* transpose = NULL;
    MATRIX_TEST_CALL(jmtxs_matrix_brm_transpose(test_matrix, &transpose, NULL));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    print_brm_matrix(transpose);
    float* pv;
    ASSERT((cnt = jmtx_matrix_brm_get_row(transpose, 0, &pv)) == 2);
    ASSERT(are_arrays_the_same(cnt, pv, (float[]){1.0f, -1.0f}));

    ASSERT((cnt = jmtx_matrix_brm_get_row(transpose, 1, &pv)) == 3);
    ASSERT(are_arrays_the_same(cnt, pv, (float[]){2.0f, -2.0f, 1.0f}));

    ASSERT((cnt = jmtx_matrix_brm_get_row(transpose, 2, &pv)) == 4);
    ASSERT(are_arrays_the_same(cnt, pv, (float[]){3.0f, -3.0f, 2.0f, -5.0f}));

    ASSERT((cnt = jmtx_matrix_brm_get_row(transpose, 3, &pv)) == 4);
    ASSERT(are_arrays_the_same(cnt, pv, (float[]){-4.0f, 3.0f, -6.0f, 69.0f}));

    ASSERT((cnt = jmtx_matrix_brm_get_row(transpose, 4, &pv)) == 3);
    ASSERT(are_arrays_the_same(cnt, pv, (float[]){4.0f, -7.0f, 420.0f}));

    ASSERT((cnt = jmtx_matrix_brm_get_col(transpose, 0, column_vals)) == 3);
    ASSERT(are_arrays_the_same(cnt, column_vals, (float[]){1.0f, 2.0f, 3.0f}));

    ASSERT((cnt = jmtx_matrix_brm_get_col(transpose, 1, column_vals)) == 4);
    ASSERT(are_arrays_the_same(cnt, column_vals, (float[]){-1.0f, -2.0f, -3.0f, -4.0f}));

    ASSERT((cnt = jmtx_matrix_brm_get_col(transpose, 2, column_vals)) == 4);
    ASSERT(are_arrays_the_same(cnt, column_vals, (float[]){1.0f, 2.0f, 3.0f, 4.0f}));

    ASSERT((cnt = jmtx_matrix_brm_get_col(transpose, 3, column_vals)) == 3);
    ASSERT(are_arrays_the_same(cnt, column_vals, (float[]){-5.0f, -6.0f, -7.0f}));

    ASSERT((cnt = jmtx_matrix_brm_get_col(transpose, 4, column_vals)) == 2);
    ASSERT(are_arrays_the_same(cnt, column_vals, (float[]){69.0f, 420.0f}));

    for (unsigned i = 0; i < 5; ++i)
    {
        for (unsigned j = 0; j < 5; ++j)
        {
            float v1, v2;
            ASSERT(jmtxs_matrix_brm_get_entry(test_matrix, i, j, &v1) == JMTX_RESULT_SUCCESS);
            ASSERT(jmtxs_matrix_brm_get_entry(transpose, j, i, &v2) == JMTX_RESULT_SUCCESS);
            ASSERT(v1 == v2);
        }
    }

    MATRIX_TEST_CALL(jmtxs_matrix_brm_destroy(transpose));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    MATRIX_TEST_CALL(jmtxs_matrix_brm_destroy(test_matrix));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    return 0;
}
