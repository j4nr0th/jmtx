// Automatically generated from tests/float/matrix_ops_tests/base_cds_test.c on Fri Dec  1 17:35:45 2023
//
// Created by jan on 24.11.2023.
//

#include "../test_common.h"
#include "../../../include/jmtx/cfloat/matrices/sparse_diagonal_compressed.h"
#include "../../../include/jmtx/cfloat/matrices/sparse_diagonal_compressed_safe.h"
#include <complex.h>

static int are_value_arrays_the_same(unsigned len, const _Complex float a1[const static len], const _Complex float a2[const static len])
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

static int are_index_arrays_the_same(unsigned len, const uint32_t a1[const static len], const uint32_t a2[const static len])
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

int test_row(void)
{
    jmtxc_matrix_cds* test_matrix;
    jmtx_result mtx_res;
    MATRIX_TEST_CALL(jmtxcs_matrix_cds_new(&test_matrix, 5, 5, 3, (const int32_t[]){-1, 0, 1}, NULL));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    print_cds_matrix(test_matrix);

    //  Set the first row
    {
        MATRIX_TEST_CALL(jmtxcs_matrix_cds_set_row(test_matrix, 0, 3, (_Complex float[]){1.0f, 2.0f, 3.0f}, (uint32_t[]){0, 1, 2}));
        ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    }
    print_cds_matrix(test_matrix);
    //  Set the second row
    {
        MATRIX_TEST_CALL(jmtxcs_matrix_cds_set_row(test_matrix, 1, 4, (_Complex float[]){-1.0f, -2.0f, -3.0f, -4.0f}, (uint32_t[]){0, 1, 2, 3}));
        ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    }
    print_cds_matrix(test_matrix);
    //  Set the third row
    {
        MATRIX_TEST_CALL(jmtxcs_matrix_cds_set_row(test_matrix, 2, 4, (_Complex float[]){1.0f, 2.0f, 3.0f, 4.0f}, (uint32_t[]){1, 2, 3, 4}));
        ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    }
    print_cds_matrix(test_matrix);
    //  Set the fourth row
    {
        MATRIX_TEST_CALL(jmtxcs_matrix_cds_set_row(test_matrix, 3, 3, (_Complex float[]){-5.0f, -6.0f, -7.0f}, (uint32_t[]){2, 3, 4}));
        ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    }
    print_cds_matrix(test_matrix);
    //  Set the fifth row
    {
        MATRIX_TEST_CALL(jmtxcs_matrix_cds_set_row(test_matrix, 4, 2, (_Complex float[]){69.0f, 420.0f}, (uint32_t[]){3, 4}));
        ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    }
    print_cds_matrix(test_matrix);

    _Complex float column_vals[5];
    uint32_t indices[5];
    uint_fast32_t cnt;

    ASSERT((cnt = jmtxc_matrix_cds_get_col(test_matrix, 0, 5, column_vals, indices)) == 2);
    ASSERT(are_value_arrays_the_same(cnt, column_vals, (_Complex float[]){1.0f, -1.0f}));
    ASSERT(are_index_arrays_the_same(cnt, indices, (const uint32_t []) {0, 1}));
    ASSERT(cnt == jmtxc_matrix_cds_entries_in_col(test_matrix, 0));

    ASSERT((cnt = jmtxc_matrix_cds_get_col(test_matrix, 1, 5, column_vals, indices)) == 3);
    ASSERT(are_value_arrays_the_same(cnt, column_vals, (_Complex float[]){2.0f, -2.0f, 1.0f}));
    ASSERT(are_index_arrays_the_same(cnt, indices, (const uint32_t []) {0, 1, 2}));
    ASSERT(cnt == jmtxc_matrix_cds_entries_in_col(test_matrix, 1));

    ASSERT((cnt = jmtxc_matrix_cds_get_col(test_matrix, 2, 5, column_vals, indices)) == 4);
    ASSERT(are_value_arrays_the_same(cnt, column_vals, (_Complex float[]){3.0f, -3.0f, 2.0f, -5.0f}));
    ASSERT(are_index_arrays_the_same(cnt, indices, (const uint32_t []) {0, 1, 2, 3}));
    ASSERT(cnt == jmtxc_matrix_cds_entries_in_col(test_matrix, 2));

    ASSERT((cnt = jmtxc_matrix_cds_get_col(test_matrix, 3, 5, column_vals, indices)) == 4);
    ASSERT(are_value_arrays_the_same(cnt, column_vals, (_Complex float[]){-4.0f, 3.0f, -6.0f, 69.0f}));
    ASSERT(are_index_arrays_the_same(cnt, indices, (const uint32_t []) {1, 2, 3, 4}));
    ASSERT(cnt == jmtxc_matrix_cds_entries_in_col(test_matrix, 3));

    ASSERT((cnt = jmtxc_matrix_cds_get_col(test_matrix, 4, 5, column_vals, indices)) == 3);
    ASSERT(are_value_arrays_the_same(cnt, column_vals, (_Complex float[]){4.0f, -7.0f, 420.0f}));
    ASSERT(are_index_arrays_the_same(cnt, indices, (const uint32_t []) {2, 3, 4}));
    ASSERT(cnt == jmtxc_matrix_cds_entries_in_col(test_matrix, 4));

    const _Complex float test_vec[5] = {1.0f, -1.0f, 3.0f, 0.0f, -3.0f};
    _Complex float out_vec[5];

    jmtxc_matrix_cds_vector_multiply(test_matrix, test_vec, out_vec);
    for (unsigned i = 0; i < 5; ++i)
    {
        printf("%g%+g ", crealf(out_vec[i]), cimagf(out_vec[i]));
    }
    printf("\n");

    jmtxc_matrix_cds* transpose = NULL;
    MATRIX_TEST_CALL(jmtxcs_matrix_cds_transpose(test_matrix, &transpose, NULL));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    print_cds_matrix(transpose);

    MATRIX_TEST_CALL(jmtxcs_matrix_cds_destroy(transpose));
    MATRIX_TEST_CALL(jmtxcs_matrix_cds_destroy(test_matrix));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    return 0;
}

int test_col(void)
{
    jmtxc_matrix_cds* test_matrix;
    jmtx_result mtx_res;
    MATRIX_TEST_CALL(jmtxcs_matrix_cds_new(&test_matrix, 5, 5, 3, (const int32_t[]){-1, 0, 1}, NULL));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    print_cds_matrix(test_matrix);

    //  Set the first col
    {
        MATRIX_TEST_CALL(jmtxcs_matrix_cds_set_col(test_matrix, 0, 3, (_Complex float[]){1.0f, 2.0f, 3.0f}, (uint32_t[]){0, 1, 2}));
        ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    }
    print_cds_matrix(test_matrix);
    //  Set the second col
    {
        MATRIX_TEST_CALL(jmtxcs_matrix_cds_set_col(test_matrix, 1, 4, (_Complex float[]){-1.0f, -2.0f, -3.0f, -4.0f}, (uint32_t[]){0, 1, 2, 3}));
        ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    }
    print_cds_matrix(test_matrix);
    //  Set the third col
    {
        MATRIX_TEST_CALL(jmtxcs_matrix_cds_set_col(test_matrix, 2, 4, (_Complex float[]){1.0f, 2.0f, 3.0f, 4.0f}, (uint32_t[]){1, 2, 3, 4}));
        ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    }
    print_cds_matrix(test_matrix);
    //  Set the fourth col
    {
        MATRIX_TEST_CALL(jmtxcs_matrix_cds_set_col(test_matrix, 3, 3, (_Complex float[]){-5.0f, -6.0f, -7.0f}, (uint32_t[]){2, 3, 4}));
        ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    }
    print_cds_matrix(test_matrix);
    //  Set the fifth col
    {
        MATRIX_TEST_CALL(jmtxcs_matrix_cds_set_col(test_matrix, 4, 2, (_Complex float[]){69.0f, 420.0f}, (uint32_t[]){3, 4}));
        ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    }
    print_cds_matrix(test_matrix);

    _Complex float column_vals[5];
    uint32_t indices[5];
    uint_fast32_t cnt;

    ASSERT((cnt = jmtxc_matrix_cds_get_row(test_matrix, 0, 5, column_vals, indices)) == 2);
    ASSERT(are_value_arrays_the_same(cnt, column_vals, (_Complex float[]){1.0f, -1.0f}));
    ASSERT(are_index_arrays_the_same(cnt, indices, (const uint32_t []) {0, 1}));
    ASSERT(cnt == jmtxc_matrix_cds_entries_in_row(test_matrix, 0));

    ASSERT((cnt = jmtxc_matrix_cds_get_row(test_matrix, 1, 5, column_vals, indices)) == 3);
    ASSERT(are_value_arrays_the_same(cnt, column_vals, (_Complex float[]){2.0f, -2.0f, 1.0f}));
    ASSERT(are_index_arrays_the_same(cnt, indices, (const uint32_t []) {0, 1, 2}));
    ASSERT(cnt == jmtxc_matrix_cds_entries_in_row(test_matrix, 1));

    ASSERT((cnt = jmtxc_matrix_cds_get_row(test_matrix, 2, 5, column_vals, indices)) == 4);
    ASSERT(are_value_arrays_the_same(cnt, column_vals, (_Complex float[]){3.0f, -3.0f, 2.0f, -5.0f}));
    ASSERT(are_index_arrays_the_same(cnt, indices, (const uint32_t []) {0, 1, 2, 3}));
    ASSERT(cnt == jmtxc_matrix_cds_entries_in_row(test_matrix, 2));

    ASSERT((cnt = jmtxc_matrix_cds_get_row(test_matrix, 3, 5, column_vals, indices)) == 4);
    ASSERT(are_value_arrays_the_same(cnt, column_vals, (_Complex float[]){-4.0f, 3.0f, -6.0f, 69.0f}));
    ASSERT(are_index_arrays_the_same(cnt, indices, (const uint32_t []) {1, 2, 3, 4}));
    ASSERT(cnt == jmtxc_matrix_cds_entries_in_row(test_matrix, 3));

    ASSERT((cnt = jmtxc_matrix_cds_get_row(test_matrix, 4, 5, column_vals, indices)) == 3);
    ASSERT(are_value_arrays_the_same(cnt, column_vals, (_Complex float[]){4.0f, -7.0f, 420.0f}));
    ASSERT(are_index_arrays_the_same(cnt, indices, (const uint32_t []) {2, 3, 4}));
    ASSERT(cnt == jmtxc_matrix_cds_entries_in_row(test_matrix, 4));


    jmtxc_matrix_cds* transpose = NULL;
    MATRIX_TEST_CALL(jmtxcs_matrix_cds_transpose(test_matrix, &transpose, NULL));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    print_cds_matrix(transpose);

    MATRIX_TEST_CALL(jmtxcs_matrix_cds_destroy(transpose));
    MATRIX_TEST_CALL(jmtxcs_matrix_cds_destroy(test_matrix));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    return 0;
}

int main()
{
    ASSERT(test_col() == 0);
    ASSERT(test_row() == 0);
    return 0;
}
