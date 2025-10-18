#include "../test_common.h"
#include "decompositions/dense_lu_decomposition.h"
#include "solvers/lu_solving.h"
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
    JMTX_NAME_TYPED(matrix_drm) * mtx, *decomp, *decompnp;
    MATRIX_TEST_CALL(JMTX_NAME_TYPED(matrix_drm_new)(&mtx, PROBLEM_DIMS, PROBLEM_DIMS, NULL, NULL));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    MATRIX_TEST_CALL(JMTX_NAME_TYPED(matrix_drm_new)(&decomp, PROBLEM_DIMS, PROBLEM_DIMS, NULL, NULL));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    MATRIX_TEST_CALL(JMTX_NAME_TYPED(matrix_drm_new)(&decompnp, PROBLEM_DIMS, PROBLEM_DIMS, NULL, NULL));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);

    JMTX_NAME_TYPED(matrix_drm_set_row)(mtx, 0, (JMTX_SCALAR_T[PROBLEM_DIMS]){1, 3, 0, 0});
    JMTX_NAME_TYPED(matrix_drm_set_row)(mtx, 1, (JMTX_SCALAR_T[PROBLEM_DIMS]){-.001f, 2, 40, 0});
    JMTX_NAME_TYPED(matrix_drm_set_row)(mtx, 2, (JMTX_SCALAR_T[PROBLEM_DIMS]){0, 1, 1, 0});
    JMTX_NAME_TYPED(matrix_drm_set_row)(mtx, 3, (JMTX_SCALAR_T[PROBLEM_DIMS]){0, 200, -400, 1});
    JMTX_NAME_TYPED(print_drm_matrix)(mtx);
    uint32_t ubw, lbw;
    JMTX_NAME_TYPED(matrix_drm_get_bandwidths)(mtx, &ubw, &lbw);
    printf("UBW: %" PRIu32 ", LBW: %" PRIu32 "\n", ubw, lbw);

    const JMTX_SCALAR_T test_vector_in1[PROBLEM_DIMS] = {1.f, -.3f, .2f, -3.f};
    JMTX_SCALAR_T test_vector_out1[PROBLEM_DIMS];
    JMTX_SCALAR_T test_vector_reverse1[PROBLEM_DIMS];
    JMTX_SCALAR_T test_vector_reverse2[PROBLEM_DIMS];
    JMTX_NAME_TYPED(matrix_drm_vector_multiply)(mtx, test_vector_in1, test_vector_out1);

    clock_gettime(CLOCK_MONOTONIC, &ts0);
    JMTX_NAME_TYPED(decompose_lu_pivot_drm)(mtx, decomp);
    clock_gettime(CLOCK_MONOTONIC, &ts1);
    printf("Time taken with pivoting DRM %g seconds\n",
           (double)(ts1.tv_sec - ts0.tv_sec) + (double)(ts1.tv_nsec - ts0.tv_nsec) * 1e-9);
    JMTX_NAME_TYPED(print_drm_matrix)(decomp);

    clock_gettime(CLOCK_MONOTONIC, &ts0);
    JMTX_NAME_TYPED(decompose_lu_drm)(mtx, decompnp);
    clock_gettime(CLOCK_MONOTONIC, &ts1);
    printf("Time taken without pivoting DRM %g seconds\n",
           (double)(ts1.tv_sec - ts0.tv_sec) + (double)(ts1.tv_nsec - ts0.tv_nsec) * 1e-9);
    JMTX_NAME_TYPED(print_drm_matrix)(decompnp);

    JMTX_NAME_TYPED(solve_direct_lu_drm)(decomp, test_vector_out1, test_vector_reverse1,
                                         (JMTX_SCALAR_T[PROBLEM_DIMS]){});
    JMTX_NAME_TYPED(solve_direct_lu_drm)(decompnp, test_vector_out1, test_vector_reverse2,
                                         (JMTX_SCALAR_T[PROBLEM_DIMS]){});

    printf("\nOriginal vector:\n");
    JMTX_NAME_TYPED(print_vec)(PROBLEM_DIMS, test_vector_in1);
    printf("\nMultiplied vector:\n");
    JMTX_NAME_TYPED(print_vec)(PROBLEM_DIMS, test_vector_out1);

    double err_with_pivot = 0, err_without_pivot = 0;
    for (uint32_t i = 0; i < PROBLEM_DIMS; ++i)
    {
        const JMTX_SCALAR_T real = test_vector_in1[i];
        const JMTX_SCALAR_T ewith = real - test_vector_reverse1[i];
        const JMTX_SCALAR_T ewithout = real - test_vector_reverse2[i];
        err_with_pivot += JMTX_DOT(ewith, ewith);
        err_without_pivot += JMTX_DOT(ewithout, ewithout);
    }
    err_without_pivot = sqrt(err_without_pivot);
    err_with_pivot = sqrt(err_with_pivot);
    printf("\nSolved vector with pivots (L2 error norm: %e):\n", err_with_pivot);
    JMTX_NAME_TYPED(print_vec)(PROBLEM_DIMS, test_vector_reverse1);
    printf("\nSolved vector without pivots (L2 error norm: %e):\n", err_without_pivot);
    JMTX_NAME_TYPED(print_vec)(PROBLEM_DIMS, test_vector_reverse2);

    JMTX_NAME_TYPED(matrix_drm_destroy)(decompnp);
    JMTX_NAME_TYPED(matrix_drm_destroy)(decomp);
    JMTX_NAME_TYPED(matrix_drm_destroy)(mtx);

    return 0;
}
