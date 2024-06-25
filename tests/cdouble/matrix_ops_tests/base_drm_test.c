//
// Created by jan on 3.1.2024.
//
#include <inttypes.h>
#include <stdio.h>
#include "../test_common.h"

#include "../../../include/jmtx/cdouble/matrices/dense_row_major.h"

enum
{
    TEST_COLS1 = 3, TEST_ROWS1 = 4,
};

int main()
{
    jmtx_result mtx_res;
    jmtxz_matrix_drm* mtx1 = NULL;
    MATRIX_TEST_CALL(jmtxz_matrix_drm_new(&mtx1, TEST_ROWS1, TEST_COLS1, NULL, NULL));
    print_drmz_matrix(mtx1);
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);

    jmtxz_matrix_drm_set_row(mtx1, 0, (_Complex double[]){ 1,  2, 4});
    jmtxz_matrix_drm_set_row(mtx1, 1, (_Complex double[]){-2, -3, 1});
    jmtxz_matrix_drm_set_row(mtx1, 2, (_Complex double[]){ 2,  3, 0});
    jmtxz_matrix_drm_set_row(mtx1, 3, (_Complex double[]){ 0,  0,-4});
    print_drmz_matrix(mtx1);

    uint32_t lbw, ubw;
    jmtxz_matrix_drm_get_bandwidths(mtx1, &ubw, &lbw);
    printf("Upper and lower bandwidths of the matrix are %"PRIu32" and %"PRIu32"\n", ubw, lbw);

    jmtxz_matrix_drm_set_permutation(mtx1, (uint32_t[]){1, 3, 0, 2});
    print_drmz_matrix(mtx1);

    _Complex double test_vec1_in[TEST_COLS1] = {2, -1, 3};
    _Complex double test_vec1_out[TEST_ROWS1] = {-69, -68, -67};

    jmtxz_matrix_drm_vector_multiply(mtx1, test_vec1_in, test_vec1_out);
    print_vecz(TEST_COLS1, test_vec1_in);
    print_vecz(TEST_ROWS1, test_vec1_out);

    jmtxz_matrix_drm_set_permutation(mtx1, (uint32_t[]){1, 2, 0, 3});
    jmtxz_matrix_drm_swap_rows(mtx1, 1, 3);
    print_drmz_matrix(mtx1);
//    _Complex double extra_row[TEST_COLS1];
//    jmtxz_matrix_drm_commit_permutations2(mtx1, extra_row);
    jmtxz_matrix_drm_commit_permutations(mtx1);
    print_drmz_matrix(mtx1);

    jmtxz_matrix_drm_destroy(mtx1);
    return 0;
}
