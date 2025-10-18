#include <inttypes.h>
#include <stdio.h>
#include "../test_common.h"

enum
{
    TEST_COLS1 = 3,
    TEST_ROWS1 = 4,
};

int main()
{
    jmtx_result mtx_res;
    JMTX_NAME_TYPED(matrix_drm) *mtx1 = NULL;
    MATRIX_TEST_CALL(JMTX_NAME_TYPED(matrix_drm_new)(&mtx1, TEST_ROWS1, TEST_COLS1, NULL, NULL));
    JMTX_NAME_TYPED(print_drm_matrix)(mtx1);
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);

    JMTX_NAME_TYPED(matrix_drm_set_row)(mtx1, 0, (JMTX_SCALAR_T[]){1, 2, 4});
    JMTX_NAME_TYPED(matrix_drm_set_row)(mtx1, 1, (JMTX_SCALAR_T[]){-2, -3, 1});
    JMTX_NAME_TYPED(matrix_drm_set_row)(mtx1, 2, (JMTX_SCALAR_T[]){2, 3, 0});
    JMTX_NAME_TYPED(matrix_drm_set_row)(mtx1, 3, (JMTX_SCALAR_T[]){0, 0, -4});
    JMTX_NAME_TYPED(print_drm_matrix)(mtx1);

    uint32_t lbw, ubw;
    JMTX_NAME_TYPED(matrix_drm_get_bandwidths)(mtx1, &ubw, &lbw);
    printf("Upper and lower bandwidths of the matrix are %" PRIu32 " and %" PRIu32 "\n", ubw, lbw);

    JMTX_NAME_TYPED(matrix_drm_set_permutation)(mtx1, (uint32_t[]){1, 3, 0, 2});
    JMTX_NAME_TYPED(print_drm_matrix)(mtx1);

    JMTX_SCALAR_T test_vec1_in[TEST_COLS1] = {2, -1, 3};
    JMTX_SCALAR_T test_vec1_out[TEST_ROWS1] = {-69, -68, -67};

    JMTX_NAME_TYPED(matrix_drm_vector_multiply)(mtx1, test_vec1_in, test_vec1_out);
    JMTX_NAME_TYPED(print_vec)(TEST_COLS1, test_vec1_in);
    JMTX_NAME_TYPED(print_vec)(TEST_ROWS1, test_vec1_out);

    JMTX_NAME_TYPED(matrix_drm_set_permutation)(mtx1, (uint32_t[]){1, 2, 0, 3});
    JMTX_NAME_TYPED(matrix_drm_swap_rows)(mtx1, 1, 3);
    JMTX_NAME_TYPED(print_drm_matrix)(mtx1);
    //    float extra_row[TEST_COLS1];
    //    jmtxd_matrix_drm_commit_permutations2(mtx1, extra_row);
    JMTX_NAME_TYPED(matrix_drm_commit_permutations)(mtx1);
    JMTX_NAME_TYPED(print_drm_matrix)(mtx1);

    JMTX_NAME_TYPED(matrix_drm_destroy)(mtx1);
    return 0;
}
