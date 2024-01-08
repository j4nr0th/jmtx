//
// Created by jan on 3.1.2024.
//
#include <inttypes.h>
#include <stdio.h>
#include "../test_common.h"

#include "../../../include/jmtx/cfloat/matrices/dense_row_major.h"

enum
{
    TEST_COLS1 = 3, TEST_ROWS1 = 4,
};

int main()
{
    jmtx_result mtx_res;
    jmtxc_matrix_drm* mtx1 = NULL;
    MATRIX_TEST_CALL(jmtxc_matrix_drm_new(&mtx1, TEST_COLS1, TEST_ROWS1, NULL, NULL));
    print_drmc_matrix(mtx1);
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);

    jmtxc_matrix_drm_set_row(mtx1, 0, (_Complex float[]){ 1,  2, 4});
    jmtxc_matrix_drm_set_row(mtx1, 1, (_Complex float[]){-2, -3, 1});
    jmtxc_matrix_drm_set_row(mtx1, 2, (_Complex float[]){ 2,  3, 0});
    jmtxc_matrix_drm_set_row(mtx1, 3, (_Complex float[]){ 0,  0,-4});
    print_drmc_matrix(mtx1);

    uint32_t lbw, ubw;
    jmtxc_matrix_drm_get_bandwidths(mtx1, &ubw, &lbw);
    printf("Upper and lower bandwidths of the matrix are %"PRIu32" and %"PRIu32"\n", ubw, lbw);

    jmtxc_matrix_drm_set_permutation(mtx1, (uint32_t[]){1, 3, 0, 2});
    print_drmc_matrix(mtx1);

    _Complex float test_vec1_in[TEST_COLS1] = {2, -1, 3};
    _Complex float test_vec1_out[TEST_ROWS1] = {-69, -68, -67};

    jmtxc_matrix_drm_vector_multiply(mtx1, test_vec1_in, test_vec1_out);
    print_vecc(TEST_COLS1, test_vec1_in);
    print_vecc(TEST_ROWS1, test_vec1_out);

    jmtxc_matrix_drm_set_permutation(mtx1, (uint32_t[]){1, 2, 0, 3});
    jmtxc_matrix_drm_swap_rows(mtx1, 1, 3);
    print_drmc_matrix(mtx1);
//    _Complex float extra_row[TEST_COLS1];
//    jmtxc_matrix_drm_commit_permutations2(mtx1, extra_row);
    jmtxc_matrix_drm_commit_permutations(mtx1);
    print_drmc_matrix(mtx1);

    jmtxc_matrix_drm_destroy(mtx1);
    return 0;
}
