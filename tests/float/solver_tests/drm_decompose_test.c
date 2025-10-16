#include "../../../include/jmtx/float/decompositions/dense_lu_decomposition.h"
#include "../../../include/jmtx/float/solvers/lu_solving.h"
#include "../test_common.h"
#include <time.h>
#include <inttypes.h>
#include <math.h>

enum
{
    PROBLEM_DIMS = 4
};

int main()
{

    struct timespec ts0, ts1;

    jmtx_result mtx_res;
    jmtxf_matrix_drm *mtx, *decomp, *decompnp;
    MATRIX_TEST_CALL(jmtxf_matrix_drm_new(&mtx, PROBLEM_DIMS, PROBLEM_DIMS, NULL, NULL));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    MATRIX_TEST_CALL(jmtxf_matrix_drm_new(&decomp, PROBLEM_DIMS, PROBLEM_DIMS, NULL, NULL));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    MATRIX_TEST_CALL(jmtxf_matrix_drm_new(&decompnp, PROBLEM_DIMS, PROBLEM_DIMS, NULL, NULL));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);

    jmtxf_matrix_drm_set_row(mtx, 0, (float[PROBLEM_DIMS]){1, 3, 0, 0});
    jmtxf_matrix_drm_set_row(mtx, 1, (float[PROBLEM_DIMS]){-.001f, 2, 40, 0});
    jmtxf_matrix_drm_set_row(mtx, 2, (float[PROBLEM_DIMS]){0, 1, 1, 0});
    jmtxf_matrix_drm_set_row(mtx, 3, (float[PROBLEM_DIMS]){0, 200, -400, 1});
    print_drm_matrix(mtx);
    uint32_t ubw, lbw;
    jmtxf_matrix_drm_get_bandwidths(mtx, &ubw, &lbw);
    printf("UBW: %" PRIu32 ", LBW: %" PRIu32 "\n", ubw, lbw);

    float test_vector_in1[PROBLEM_DIMS] = {1.f, -.3f, .2f, -3.f};
    float test_vector_out1[PROBLEM_DIMS];
    float test_vector_reverse1[PROBLEM_DIMS];
    float test_vector_reverse2[PROBLEM_DIMS];
    jmtxf_matrix_drm_vector_multiply(mtx, test_vector_in1, test_vector_out1);

    clock_gettime(CLOCK_MONOTONIC, &ts0);
    jmtxf_decompose_lu_pivot_drm(mtx, decomp);
    clock_gettime(CLOCK_MONOTONIC, &ts1);
    printf("Time taken with pivoting DRM %g seconds\n",
           (double)(ts1.tv_sec - ts0.tv_sec) + (double)(ts1.tv_nsec - ts0.tv_nsec) * 1e-9);
    print_drm_matrix(decomp);

    clock_gettime(CLOCK_MONOTONIC, &ts0);
    jmtxf_decompose_lu_drm(mtx, decompnp);
    clock_gettime(CLOCK_MONOTONIC, &ts1);
    printf("Time taken without pivoting DRM %g seconds\n",
           (double)(ts1.tv_sec - ts0.tv_sec) + (double)(ts1.tv_nsec - ts0.tv_nsec) * 1e-9);
    print_drm_matrix(decompnp);

    jmtx_solve_direct_lu_drm(decomp, test_vector_out1, test_vector_reverse1, (float[PROBLEM_DIMS]){});
    jmtx_solve_direct_lu_drm(decompnp, test_vector_out1, test_vector_reverse2, (float[PROBLEM_DIMS]){});

    printf("\nOriginal vector:\n");
    print_vec(PROBLEM_DIMS, test_vector_in1);
    printf("\nMultiplied vector:\n");
    print_vec(PROBLEM_DIMS, test_vector_out1);

    double err_with_pivot = 0, err_without_pivot = 0;
    for (uint32_t i = 0; i < PROBLEM_DIMS; ++i)
    {
        const double real = (double)test_vector_in1[i];
        const double ewith = real - (double)test_vector_reverse1[i];
        const double ewithout = real - (double)test_vector_reverse2[i];
        err_with_pivot += ewith * ewith;
        err_without_pivot += ewithout * ewithout;
    }
    err_without_pivot = sqrt(err_without_pivot);
    err_with_pivot = sqrt(err_with_pivot);
    printf("\nSolved vector with pivots (L2 error norm: %e):\n", err_with_pivot);
    print_vec(PROBLEM_DIMS, test_vector_reverse1);
    printf("\nSolved vector without pivots (L2 error norm: %e):\n", err_without_pivot);
    print_vec(PROBLEM_DIMS, test_vector_reverse2);

    jmtxf_matrix_drm_destroy(decompnp);
    jmtxf_matrix_drm_destroy(decomp);
    jmtxf_matrix_drm_destroy(mtx);

    return 0;
}
