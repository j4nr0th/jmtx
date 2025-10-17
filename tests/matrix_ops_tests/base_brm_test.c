#include INCLUDE_HEADER
#include "../test_common.h"

static int are_arrays_the_same(unsigned len, const JMTX_SCALAR_T a1[JMTX_ARRAY_ATTRIB(const static len)],
                               const JMTX_SCALAR_T a2[JMTX_ARRAY_ATTRIB(const static len)])
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
    JMTX_NAME_TYPED(matrix_brm) * test_matrix;
    jmtx_result mtx_res;
    MATRIX_TEST_CALL(JMTX_NAME_TYPED(matrix_brm_new)(&test_matrix, 5, 5, 2, 1, NULL, NULL));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);

    //  Set the first row
    JMTX_NAME_TYPED(matrix_brm_set_row)(test_matrix, 0, (JMTX_SCALAR_T[]){1.0f, 2.0f, 3.0f});
    //  Set the second row
    JMTX_NAME_TYPED(matrix_brm_set_row)(test_matrix, 1, (JMTX_SCALAR_T[]){-1.0f, -2.0f, -3.0f, -4.0f});
    //  Set the third row
    JMTX_NAME_TYPED(matrix_brm_set_row)(test_matrix, 2, (JMTX_SCALAR_T[]){1.0f, 2.0f, 3.0f, 4.0f});
    //  Set the fourth row
    JMTX_NAME_TYPED(matrix_brm_set_row)(test_matrix, 3, (JMTX_SCALAR_T[]){-5.0f, -6.0f, -7.0f});
    //  Set the fifth row
    JMTX_NAME_TYPED(matrix_brm_set_row)(test_matrix, 4, (JMTX_SCALAR_T[]){69.0f, 420.0f});
    JMTX_NAME_TYPED(print_brm_matrix)(test_matrix);

    JMTX_SCALAR_T column_vals[4];
    uint_fast32_t cnt;

    ASSERT((cnt = JMTX_NAME_TYPED(matrix_brm_get_col)(test_matrix, 0, column_vals)) == 2);
    ASSERT(are_arrays_the_same(cnt, column_vals, (JMTX_SCALAR_T[]){1.0f, -1.0f}));

    ASSERT((cnt = JMTX_NAME_TYPED(matrix_brm_get_col)(test_matrix, 1, column_vals)) == 3);
    ASSERT(are_arrays_the_same(cnt, column_vals, (JMTX_SCALAR_T[]){2.0f, -2.0f, 1.0f}));

    ASSERT((cnt = JMTX_NAME_TYPED(matrix_brm_get_col)(test_matrix, 2, column_vals)) == 4);
    ASSERT(are_arrays_the_same(cnt, column_vals, (JMTX_SCALAR_T[]){3.0f, -3.0f, 2.0f, -5.0f}));

    ASSERT((cnt = JMTX_NAME_TYPED(matrix_brm_get_col)(test_matrix, 3, column_vals)) == 4);
    ASSERT(are_arrays_the_same(cnt, column_vals, (JMTX_SCALAR_T[]){-4.0f, 3.0f, -6.0f, 69.0f}));

    ASSERT((cnt = JMTX_NAME_TYPED(matrix_brm_get_col)(test_matrix, 4, column_vals)) == 3);
    ASSERT(are_arrays_the_same(cnt, column_vals, (JMTX_SCALAR_T[]){4.0f, -7.0f, 420.0f}));

    const JMTX_SCALAR_T test_vec[5] = {1.0f, -1.0f, 3.0f, 0.0f, -3.0f};
    JMTX_SCALAR_T out_vec[5];

    JMTX_NAME_TYPED(matrix_brm_vector_multiply)(test_matrix, test_vec, out_vec);
    for (unsigned i = 0; i < 5; ++i)
    {
        printf("%g ", out_vec[i]);
    }
    printf("\n");

    JMTX_NAME_TYPED(matrix_brm) *transpose = NULL;
    MATRIX_TEST_CALL(JMTX_NAME_TYPED(matrix_brm_transpose)(test_matrix, &transpose, NULL));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    JMTX_NAME_TYPED(print_brm_matrix)(transpose);
    JMTX_SCALAR_T *pv;
    ASSERT((cnt = JMTX_NAME_TYPED(matrix_brm_get_row)(transpose, 0, &pv)) == 2);
    ASSERT(are_arrays_the_same(cnt, pv, (JMTX_SCALAR_T[]){1.0f, -1.0f}));

    ASSERT((cnt = JMTX_NAME_TYPED(matrix_brm_get_row)(transpose, 1, &pv)) == 3);
    ASSERT(are_arrays_the_same(cnt, pv, (JMTX_SCALAR_T[]){2.0f, -2.0f, 1.0f}));

    ASSERT((cnt = JMTX_NAME_TYPED(matrix_brm_get_row)(transpose, 2, &pv)) == 4);
    ASSERT(are_arrays_the_same(cnt, pv, (JMTX_SCALAR_T[]){3.0f, -3.0f, 2.0f, -5.0f}));

    ASSERT((cnt = JMTX_NAME_TYPED(matrix_brm_get_row)(transpose, 3, &pv)) == 4);
    ASSERT(are_arrays_the_same(cnt, pv, (JMTX_SCALAR_T[]){-4.0f, 3.0f, -6.0f, 69.0f}));

    ASSERT((cnt = JMTX_NAME_TYPED(matrix_brm_get_row)(transpose, 4, &pv)) == 3);
    ASSERT(are_arrays_the_same(cnt, pv, (JMTX_SCALAR_T[]){4.0f, -7.0f, 420.0f}));

    ASSERT((cnt = JMTX_NAME_TYPED(matrix_brm_get_col)(transpose, 0, column_vals)) == 3);
    ASSERT(are_arrays_the_same(cnt, column_vals, (JMTX_SCALAR_T[]){1.0f, 2.0f, 3.0f}));

    ASSERT((cnt = JMTX_NAME_TYPED(matrix_brm_get_col)(transpose, 1, column_vals)) == 4);
    ASSERT(are_arrays_the_same(cnt, column_vals, (JMTX_SCALAR_T[]){-1.0f, -2.0f, -3.0f, -4.0f}));

    ASSERT((cnt = JMTX_NAME_TYPED(matrix_brm_get_col)(transpose, 2, column_vals)) == 4);
    ASSERT(are_arrays_the_same(cnt, column_vals, (JMTX_SCALAR_T[]){1.0f, 2.0f, 3.0f, 4.0f}));

    ASSERT((cnt = JMTX_NAME_TYPED(matrix_brm_get_col)(transpose, 3, column_vals)) == 3);
    ASSERT(are_arrays_the_same(cnt, column_vals, (JMTX_SCALAR_T[]){-5.0f, -6.0f, -7.0f}));

    ASSERT((cnt = JMTX_NAME_TYPED(matrix_brm_get_col)(transpose, 4, column_vals)) == 2);
    ASSERT(are_arrays_the_same(cnt, column_vals, (JMTX_SCALAR_T[]){69.0f, 420.0f}));

    for (unsigned i = 0; i < 5; ++i)
    {
        for (unsigned j = 0; j < 5; ++j)
        {
            const JMTX_SCALAR_T v1 = JMTX_NAME_TYPED(matrix_brm_get_entry)(test_matrix, i, j);
            const JMTX_SCALAR_T v2 = JMTX_NAME_TYPED(matrix_brm_get_entry)(transpose, j, i);
            ASSERT(v1 == v2);
        }
    }

    JMTX_NAME_TYPED(matrix_brm_destroy)(transpose);
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    JMTX_NAME_TYPED(matrix_brm_destroy)(test_matrix);
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    return 0;
}
