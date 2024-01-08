#include "../../../include/jmtx/cfloat/decompositions/dense_lu_decomposition.h"
#include "../../../include/jmtx/cfloat/solvers/lu_solving.h"
#include "../test_common.h"
#include <time.h>
#include <inttypes.h>
#include <math.h>
#include <complex.h>

enum {PROBLEM_DIMS = 4};

int main()
{

    struct timespec ts0, ts1;

    jmtx_result mtx_res;
    jmtxc_matrix_drm* mtx, *decomp, *decompnp;
    MATRIX_TEST_CALL(jmtxc_matrix_drm_new(&mtx, PROBLEM_DIMS, PROBLEM_DIMS, NULL, NULL));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    MATRIX_TEST_CALL(jmtxc_matrix_drm_new(&decomp, PROBLEM_DIMS, PROBLEM_DIMS, NULL, NULL));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    MATRIX_TEST_CALL(jmtxc_matrix_drm_new(&decompnp, PROBLEM_DIMS, PROBLEM_DIMS, NULL, NULL));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);

    jmtxc_matrix_drm_set_row(mtx, 0, (_Complex float[PROBLEM_DIMS]){ 1, 3, 0, 0});
    jmtxc_matrix_drm_set_row(mtx, 1, (_Complex float[PROBLEM_DIMS]){-.001f, 2, 40, 0});
    jmtxc_matrix_drm_set_row(mtx, 2, (_Complex float[PROBLEM_DIMS]){ 0, 1, 1, 0});
    jmtxc_matrix_drm_set_row(mtx, 3, (_Complex float[PROBLEM_DIMS]){ 0, 200,-400, 1});
    print_drmc_matrix(mtx);
    uint32_t ubw, lbw;
    jmtxc_matrix_drm_get_bandwidths(mtx, &ubw, &lbw);
    printf("UBW: %"PRIu32", LBW: %"PRIu32"\n", ubw, lbw);

    _Complex float test_vector_in1[PROBLEM_DIMS] = {1.f, -.3f, .2f, -3.f};
    _Complex float test_vector_out1[PROBLEM_DIMS];
    _Complex float test_vector_reverse1[PROBLEM_DIMS];
    _Complex float test_vector_reverse2[PROBLEM_DIMS];
    jmtxc_matrix_drm_vector_multiply(mtx, test_vector_in1, test_vector_out1);


    clock_gettime(CLOCK_MONOTONIC, &ts0);
    jmtxc_decompose_lu_pivot_drm(mtx, decomp);
    clock_gettime(CLOCK_MONOTONIC, &ts1);
    printf("Time taken with pivoting DRM %g seconds\n", (double)(ts1.tv_sec - ts0.tv_sec) + (double)(ts1.tv_nsec - ts0.tv_nsec) * 1e-9);
    print_drmc_matrix(decomp);

    clock_gettime(CLOCK_MONOTONIC, &ts0);
    jmtxc_decompose_lu_drm(mtx, decompnp);
    clock_gettime(CLOCK_MONOTONIC, &ts1);
    printf("Time taken without pivoting DRM %g seconds\n", (double)(ts1.tv_sec - ts0.tv_sec) + (double)(ts1.tv_nsec - ts0.tv_nsec) * 1e-9);
    print_drmc_matrix(decompnp);

    jmtxc_solve_direct_lu_drm(decomp, test_vector_out1, test_vector_reverse1, (_Complex float[PROBLEM_DIMS]){});
    jmtxc_solve_direct_lu_drm(decompnp, test_vector_out1, test_vector_reverse2, (_Complex float[PROBLEM_DIMS]){});

    printf("\nOriginal vector:\n");
    print_vecc(PROBLEM_DIMS, test_vector_in1);
    printf("\nMultiplied vector:\n");
    print_vecc(PROBLEM_DIMS, test_vector_out1);

    float err_with_pivot = 0, err_without_pivot = 0;
    for (uint32_t i = 0; i < PROBLEM_DIMS; ++i)
    {
        const _Complex float real = (_Complex float)test_vector_in1[i];
        const _Complex float ewith = real - (_Complex float)test_vector_reverse1[i];
        const _Complex float ewithout = real - (_Complex float)test_vector_reverse2[i];
        err_with_pivot += conjf(ewith) * ewith;
        err_without_pivot += conjf(ewithout) * ewithout;
    }
    err_without_pivot = sqrt(err_without_pivot);
    err_with_pivot = sqrt(err_with_pivot);
    printf("\nSolved vector with pivots (L2 error norm: %e):\n", err_with_pivot);
    print_vecc(PROBLEM_DIMS, test_vector_reverse1);
    printf("\nSolved vector without pivots (L2 error norm: %e):\n", err_without_pivot);
    print_vecc(PROBLEM_DIMS, test_vector_reverse2);


    jmtxc_matrix_drm_destroy(decompnp);
    jmtxc_matrix_drm_destroy(decomp);
    jmtxc_matrix_drm_destroy(mtx);

    return 0;
}
