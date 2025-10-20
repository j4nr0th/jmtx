#include "../test_common.h"
#include "matrices/sparse_conversion.h"
#include "matrices/sparse_multiplication.h"
#include "solvers/lu_solving.h"
#include "decompositions/band_lu_decomposition.h"
#include <float.h>
#include <math.h>
#include <assert.h>

static const double default_r_tol = FLT_EPSILON * 400;
static const double default_a_tol = FLT_EPSILON * 100;

enum
{
    PROBLEM_SIZE = 5
};

int main()
{
    JMTX_NAME_TYPED(matrix_crs) * lower, *upper, *combined, *multiplied;
    JMTX_NAME_TYPED(matrix_brm) * lower_brm, *upper_brm, *combined_brm, *multiplied_brm;
    JMTX_NAME_TYPED(matrix_ccs) * cu;
    jmtx_result mtx_res;

    MATRIX_TEST_CALL(JMTX_NAME_TYPED(matrix_crs_new)(&lower, PROBLEM_SIZE, PROBLEM_SIZE, 0, NULL));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    MATRIX_TEST_CALL(JMTX_NAME_TYPED(matrix_crs_new)(&upper, PROBLEM_SIZE, PROBLEM_SIZE, 0, NULL));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    MATRIX_TEST_CALL(JMTX_NAME_TYPED(matrix_crs_new)(&combined, PROBLEM_SIZE, PROBLEM_SIZE, 0, NULL));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);

    MATRIX_TEST_CALL(JMTX_NAME_TYPED(matrix_brm_new)(&lower_brm, PROBLEM_SIZE, PROBLEM_SIZE, 0, 4, 0, NULL));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    MATRIX_TEST_CALL(JMTX_NAME_TYPED(matrix_brm_new)(&upper_brm, PROBLEM_SIZE, PROBLEM_SIZE, 4, 0, 0, NULL));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    MATRIX_TEST_CALL(JMTX_NAME_TYPED(matrix_brm_new)(&combined_brm, PROBLEM_SIZE, PROBLEM_SIZE, 4, 4, 0, NULL));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);

    //  Make L based on predefined values
    {
        JMTX_SCALAR_T v1[1] = {1.0f};
        uint32_t i1[1] = {0};

        JMTX_SCALAR_T v2[2] = {3.0f, 1.0f};
        uint32_t i2[2] = {0, 1};

        JMTX_SCALAR_T v3[1] = {1.0f};
        uint32_t i3[1] = {2};

        JMTX_SCALAR_T v4[4] = {2.0f, 1.0f, -2.0f, 1.0f};
        uint32_t i4[4] = {0, 1, 2, 3};

        JMTX_SCALAR_T v5[3] = {1.0f, -3.0f, 1.0f};
        uint32_t i5[3] = {0, 2, 4};

        uint32_t counts[PROBLEM_SIZE] = {1, 2, 1, 4, 3};
        uint32_t *indices[PROBLEM_SIZE] = {i1, i2, i3, i4, i5};
        JMTX_SCALAR_T *values[PROBLEM_SIZE] = {v1, v2, v3, v4, v5};

        for (uint32_t i = 0; i < PROBLEM_SIZE; ++i)
        {
            mtx_res = JMTX_NAME_TYPED(matrix_crs_build_row)(lower, i, counts[i], indices[i], values[i]);
            ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
        }
    }
    {
        JMTX_NAME_TYPED(matrix_brm_set_row)(lower_brm, 0, (JMTX_SCALAR_T[]){1.0f});
        JMTX_NAME_TYPED(matrix_brm_set_row)(lower_brm, 1, (JMTX_SCALAR_T[]){3.0f, 1.0f});
        JMTX_NAME_TYPED(matrix_brm_set_row)(lower_brm, 2, (JMTX_SCALAR_T[]){0.0f, 0.0f, 1.0f});
        JMTX_NAME_TYPED(matrix_brm_set_row)(lower_brm, 3, (JMTX_SCALAR_T[]){2.0f, 1.0f, -2.0f, 1.0f});
        JMTX_NAME_TYPED(matrix_brm_set_row)(lower_brm, 4, (JMTX_SCALAR_T[]){1.0f, 0.0f, -3.0f, 0.0f, 1.0f});
    }
    JMTX_NAME_TYPED(print_crs_matrix)(lower);
    JMTX_NAME_TYPED(print_brm_matrix)(lower_brm);

    //  Make U based on predefined values
    {
        JMTX_SCALAR_T v1[2] = {3.0f, 2.0f};
        uint32_t i1[2] = {0, 4};

        JMTX_SCALAR_T v2[2] = {2.0f, 1.0f};
        uint32_t i2[2] = {1, 2};

        JMTX_SCALAR_T v3[2] = {1.0f, 3.0f};
        uint32_t i3[2] = {2, 4};

        JMTX_SCALAR_T v4[2] = {4.0f, 1.0f};
        uint32_t i4[2] = {3, 4};

        JMTX_SCALAR_T v5[1] = {-1.0f};
        uint32_t i5[1] = {4};

        uint32_t counts[PROBLEM_SIZE] = {2, 2, 2, 2, 1};
        uint32_t *indices[PROBLEM_SIZE] = {i1, i2, i3, i4, i5};
        JMTX_SCALAR_T *values[PROBLEM_SIZE] = {v1, v2, v3, v4, v5};

        for (uint32_t i = 0; i < PROBLEM_SIZE; ++i)
        {
            mtx_res = JMTX_NAME_TYPED(matrix_crs_build_row)(upper, i, counts[i], indices[i], values[i]);
            ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
        }
    }
    {
        JMTX_NAME_TYPED(matrix_brm_set_row)(upper_brm, 0, (JMTX_SCALAR_T[]){3.0f, 0.0f, 0.0f, 0.0f, 2.0f});
        JMTX_NAME_TYPED(matrix_brm_set_row)(upper_brm, 1, (JMTX_SCALAR_T[]){2.0f, 1.0f, 0.0f, 0.0f});
        JMTX_NAME_TYPED(matrix_brm_set_row)(upper_brm, 2, (JMTX_SCALAR_T[]){1.0f, 0.0f, 3.0f});
        JMTX_NAME_TYPED(matrix_brm_set_row)(upper_brm, 3, (JMTX_SCALAR_T[]){4.0f, 1.0f});
        JMTX_NAME_TYPED(matrix_brm_set_row)(upper_brm, 4, (JMTX_SCALAR_T[]){-1.0f});
    }
    JMTX_NAME_TYPED(print_crs_matrix)(upper);
    JMTX_NAME_TYPED(print_brm_matrix)(upper_brm);

    MATRIX_TEST_CALL(JMTX_NAME_TYPED(convert_crs_to_ccs)(upper, &cu, NULL));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);

    MATRIX_TEST_CALL(JMTX_NAME_TYPED(multiply_matrix_brm)(lower_brm, upper_brm, &multiplied_brm, NULL));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    MATRIX_TEST_CALL(JMTX_NAME_TYPED(multiply_matrix_crs)(lower, cu, &multiplied, NULL));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    JMTX_NAME_TYPED(matrix_ccs_destroy)(cu);

    const JMTX_SCALAR_T exact_multiplied[PROBLEM_SIZE][PROBLEM_SIZE] = {
        {3, 0, 0, 0, 2}, {9, 2, 1, 0, 6}, {0, 0, 1, 0, 3}, {6, 2, -1, 4, -1}, {3, 0, -3, 0, -8}};
    for (unsigned i = 0; i < PROBLEM_SIZE; ++i)
    {
        JMTX_NAME_TYPED(matrix_brm_set_row)(combined_brm, i, exact_multiplied[i]);
        for (unsigned j = 0; j < PROBLEM_SIZE; ++j)
        {
            ASSERT(are_close(JMTX_NAME_TYPED(matrix_crs_get_entry)(multiplied, i, j), exact_multiplied[i][j],
                             default_r_tol, default_r_tol));
            ASSERT(are_close(JMTX_NAME_TYPED(matrix_brm_get_entry)(multiplied_brm, i, j), exact_multiplied[i][j],
                             default_r_tol, default_r_tol));
        }
    }
    JMTX_NAME_TYPED(print_brm_matrix)(multiplied_brm);
    JMTX_NAME_TYPED(print_crs_matrix)(multiplied);
    JMTX_NAME_TYPED(print_brm_matrix)(combined_brm);

    JMTX_NAME_TYPED(matrix_brm) * du, *dl;
    MATRIX_TEST_CALL(JMTX_NAME_TYPED(decompose_lu_brm)(combined_brm, &dl, &du, NULL));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    JMTX_NAME_TYPED(print_crs_matrix)(upper);
    JMTX_NAME_TYPED(print_brm_matrix)(du);
    JMTX_NAME_TYPED(print_crs_matrix)(lower);
    JMTX_NAME_TYPED(print_brm_matrix)(dl);

    for (unsigned i = 0; i < PROBLEM_SIZE; ++i)
    {
        for (unsigned j = 0; j < PROBLEM_SIZE; ++j)
        {
            ASSERT(are_close(JMTX_NAME_TYPED(matrix_brm_get_entry)(du, i, j),
                             JMTX_NAME_TYPED(matrix_brm_get_entry)(upper_brm, i, j), default_r_tol, default_r_tol));
            ASSERT(are_close(JMTX_NAME_TYPED(matrix_brm_get_entry)(dl, i, j),
                             JMTX_NAME_TYPED(matrix_brm_get_entry)(lower_brm, i, j), default_r_tol, default_r_tol));
        }
    }

    JMTX_NAME_TYPED(matrix_brm_destroy)(dl);
    JMTX_NAME_TYPED(matrix_brm_destroy)(du);

    const JMTX_SCALAR_T x_exact[PROBLEM_SIZE] = {1.0f, -2.0f, 3.0f, -4.0f, 5.0f};
    const JMTX_SCALAR_T y_exact[PROBLEM_SIZE] = {13, 38, 18, -22, -46};
    JMTX_SCALAR_T y[PROBLEM_SIZE];
    JMTX_SCALAR_T x[PROBLEM_SIZE];
    const JMTX_SCALAR_T xb_exact[PROBLEM_SIZE] = {1.0f, -2.0f, 3.0f, -4.0f, 5.0f};
    const JMTX_SCALAR_T yb_exact[PROBLEM_SIZE] = {13, 38, 18, -22, -46};
    JMTX_SCALAR_T yb[PROBLEM_SIZE];
    JMTX_SCALAR_T xb[PROBLEM_SIZE];

    JMTX_NAME_TYPED(matrix_crs_vector_multiply)(multiplied, x_exact, y);
    JMTX_NAME_TYPED(matrix_brm_vector_multiply)(multiplied_brm, xb_exact, yb);
    for (unsigned i = 0; i < PROBLEM_SIZE; ++i)
    {
        printf("|yex_%u - y_%u|: %e\n", i, i, (double)JMTX_ABS(y_exact[i] - y[i]));
        ASSERT(are_close(y_exact[i], y[i], default_r_tol, default_a_tol));
        ASSERT(are_close(yb_exact[i], yb[i], default_r_tol, default_a_tol));
    }

    JMTX_NAME_TYPED(solve_direct_lu_crs)(lower, upper, y, x);
    JMTX_NAME_TYPED(solve_direct_lu_brm)(lower_brm, upper_brm, yb, xb);

    for (unsigned i = 0; i < PROBLEM_SIZE; ++i)
    {
        printf("|xex_%u - x_%u|: %e\n", i, i, (double)JMTX_ABS(x_exact[i] - x[i]));
        ASSERT(are_close(x_exact[i], x[i], default_r_tol, default_a_tol));
        ASSERT(are_close(xb_exact[i], xb[i], default_r_tol, default_a_tol));
    }

    JMTX_NAME_TYPED(matrix_crs_destroy)(multiplied);
    JMTX_NAME_TYPED(matrix_crs_destroy)(combined);
    JMTX_NAME_TYPED(matrix_crs_destroy)(upper);
    JMTX_NAME_TYPED(matrix_crs_destroy)(lower);

    JMTX_NAME_TYPED(matrix_brm_destroy)(multiplied_brm);
    JMTX_NAME_TYPED(matrix_brm_destroy)(combined_brm);
    JMTX_NAME_TYPED(matrix_brm_destroy)(upper_brm);
    JMTX_NAME_TYPED(matrix_brm_destroy)(lower_brm);

    return 0;
}
