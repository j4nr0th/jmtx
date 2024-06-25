//
// Created by jan on 3.1.2024.
//
#include <inttypes.h>
#include <stdio.h>
#include "../test_common.h"

#include "../../../include/jmtx/double/matrices/dense_row_major.h"

enum
{
    TEST_COLS1 = 3, TEST_ROWS1 = 4,
};

int main()
{
    jmtx_result mtx_res;
    jmtxd_matrix_drm* mtx1 = NULL;
    MATRIX_TEST_CALL(jmtxd_matrix_drm_new(&mtx1, TEST_ROWS1, TEST_COLS1, NULL, NULL));
    print_drmd_matrix(mtx1);
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);

    jmtxd_matrix_drm_set_row(mtx1, 0, (double[]){ 1,  2, 4});
    jmtxd_matrix_drm_set_row(mtx1, 1, (double[]){-2, -3, 1});
    jmtxd_matrix_drm_set_row(mtx1, 2, (double[]){ 2,  3, 0});
    jmtxd_matrix_drm_set_row(mtx1, 3, (double[]){ 0,  0,-4});
    print_drmd_matrix(mtx1);

    uint32_t lbw, ubw;
    jmtxd_matrix_drm_get_bandwidths(mtx1, &ubw, &lbw);
    printf("Upper and lower bandwidths of the matrix are %"PRIu32" and %"PRIu32"\n", ubw, lbw);

    jmtxd_matrix_drm_set_permutation(mtx1, (uint32_t[]){1, 3, 0, 2});
    print_drmd_matrix(mtx1);

    double test_vec1_in[TEST_COLS1] = {2, -1, 3};
    double test_vec1_out[TEST_ROWS1] = {-69, -68, -67};

    jmtxd_matrix_drm_vector_multiply(mtx1, test_vec1_in, test_vec1_out);
    print_vecd(TEST_COLS1, test_vec1_in);
    print_vecd(TEST_ROWS1, test_vec1_out);

    jmtxd_matrix_drm_set_permutation(mtx1, (uint32_t[]){1, 2, 0, 3});
    jmtxd_matrix_drm_swap_rows(mtx1, 1, 3);
    print_drmd_matrix(mtx1);
//    float extra_row[TEST_COLS1];
//    jmtxd_matrix_drm_commit_permutations2(mtx1, extra_row);
    jmtxd_matrix_drm_commit_permutations(mtx1);
    print_drmd_matrix(mtx1);

    jmtxd_matrix_drm_destroy(mtx1);
    return 0;
}
