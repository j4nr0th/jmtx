//
// Created by jan on 3.1.2024.
//
#include <inttypes.h>
#include <stdio.h>
#include "../test_common.h"

#include "../../../include/jmtx/float/matrices/dense_row_major.h"

enum
{
    TEST_COLS1 = 3, TEST_ROWS1 = 4,
};

int main()
{
    jmtx_result mtx_res;
    jmtx_matrix_drm* mtx1 = NULL;
    MATRIX_TEST_CALL(jmtx_matrix_drm_new(&mtx1, TEST_ROWS1, TEST_COLS1, NULL, NULL));
    print_drm_matrix(mtx1);
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);

    jmtx_matrix_drm_set_row(mtx1, 0, (float[]){ 1,  2, 4});
    jmtx_matrix_drm_set_row(mtx1, 1, (float[]){-2, -3, 1});
    jmtx_matrix_drm_set_row(mtx1, 2, (float[]){ 2,  3, 0});
    jmtx_matrix_drm_set_row(mtx1, 3, (float[]){ 0,  0,-4});
    print_drm_matrix(mtx1);

    uint32_t lbw, ubw;
    jmtx_matrix_drm_get_bandwidths(mtx1, &ubw, &lbw);
    printf("Upper and lower bandwidths of the matrix are %"PRIu32" and %"PRIu32"\n", ubw, lbw);

    jmtx_matrix_drm_set_permutation(mtx1, (uint32_t[]){1, 3, 0, 2});
    print_drm_matrix(mtx1);

    float test_vec1_in[TEST_COLS1] = {2, -1, 3};
    float test_vec1_out[TEST_ROWS1] = {-69, -68, -67};

    jmtx_matrix_drm_vector_multiply(mtx1, test_vec1_in, test_vec1_out);
    print_vec(TEST_COLS1, test_vec1_in);
    print_vec(TEST_ROWS1, test_vec1_out);

    jmtx_matrix_drm_set_permutation(mtx1, (uint32_t[]){1, 2, 0, 3});
    jmtx_matrix_drm_swap_rows(mtx1, 1, 3);
    print_drm_matrix(mtx1);
//    float extra_row[TEST_COLS1];
//    jmtx_matrix_drm_commit_permutations2(mtx1, extra_row);
    jmtx_matrix_drm_commit_permutations(mtx1);
    print_drm_matrix(mtx1);

    jmtx_matrix_drm_destroy(mtx1);
    return 0;
}
