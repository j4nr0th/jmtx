// Automatically generated from tests/float/matrix_ops_tests/lu_solve_test.c on Fri Dec  1 17:35:45 2023
//
// Created by jan on 23.11.2023.
//
#include "../test_common.h"
#include "../../../include/jmtx/cfloat/matrices/sparse_row_compressed_safe.h"
#include "../../../include/jmtx/cfloat/matrices/sparse_column_compressed_safe.h"
#include "../../../include/jmtx/cfloat/solvers/lu_solving.h"
#include "../../../include/jmtx/cfloat/matrices/sparse_conversion.h"
#include "../../../include/jmtx/cfloat/matrices/sparse_multiplication.h"
#include <float.h>
#include <math.h>
#include <assert.h>
#include <complex.h>
int are_close(_Complex float v1, _Complex float v2, float relative_tol, float abs_tol)
{
    if (v1 == v2)
    {
        return 1;
    }
    assert((abs_tol) >= 0.0f);
    assert((relative_tol) >= 0.0f);
    const _Complex float dif = (v1 - v2);
    if (cabsf(dif) < abs_tol)
    {
        return 1;
    }
    if (cabsf(v1) > cabsf(v2))
    {
        if (cabsf(v2 + dif) > cabsf(v1))
        {
            return 1;
        }
    }
    else
    {
        if (cabsf(v1 + dif) > cabsf(v2))
        {
            return 1;
        }
    }
    return 0;
}

static const _Complex float default_r_tol = FLT_EPSILON * 400;
static const _Complex float default_a_tol = FLT_EPSILON * 100;

enum {PROBLEM_SIZE = 5};

int main()
{
    jmtxc_matrix_crs* lower, *upper, *combined, *multiplied;
    jmtxc_matrix_ccs* cu;
    jmtx_result mtx_res;

    MATRIX_TEST_CALL(jmtxcs_matrix_crs_new(&lower, PROBLEM_SIZE, PROBLEM_SIZE, 0, NULL));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    MATRIX_TEST_CALL(jmtxcs_matrix_crs_new(&upper, PROBLEM_SIZE, PROBLEM_SIZE, 0, NULL));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    MATRIX_TEST_CALL(jmtxcs_matrix_crs_new(&combined, PROBLEM_SIZE, PROBLEM_SIZE, 0, NULL));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);

    //  Make L based on predefined values
    {
        _Complex float v1[1] = {1.0f};
        uint32_t i1[1] = {0};


        _Complex float v2[2] = {3.0f, 1.0f};
        uint32_t i2[2] = {0, 1};

        _Complex float v3[1] = {1.0f};
        uint32_t i3[1] = {2};

        _Complex float v4[4] = {2.0f, 1.0f, -2.0f, 1.0f};
        uint32_t i4[4] = {0, 1, 2, 3};

        _Complex float v5[3] = {1.0f, -3.0f, 1.0f};
        uint32_t i5[3] = {0, 2, 4};

        uint32_t counts[PROBLEM_SIZE] =
                {
                1, 2, 1, 4, 3
                };
        uint32_t* indices[PROBLEM_SIZE] =
                {
                i1, i2, i3, i4, i5
                };
        _Complex float* values[PROBLEM_SIZE] =
                {
                v1, v2, v3, v4, v5
                };

        for (uint32_t i = 0; i < PROBLEM_SIZE; ++i)
        {
            mtx_res = jmtxc_matrix_crs_build_row(lower, i, counts[i], indices[i], values[i]);
            ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
        }
    }
    print_crsc_matrix(lower);

    //  Make U based on predefined values
    {
        _Complex float v1[2] = {3.0f, 2.0f};
        uint32_t i1[2] = {0, 4};

        _Complex float v2[2] = {2.0f, 1.0f};
        uint32_t i2[2] = {1, 2};

        _Complex float v3[2] = {1.0f, 3.0f};
        uint32_t i3[2] = {2, 4};

        _Complex float v4[2] = {4.0f, 1.0f};
        uint32_t i4[2] = {3, 4};

        _Complex float v5[1] = {-1.0f};
        uint32_t i5[1] = {4};

        uint32_t counts[PROBLEM_SIZE] =
                {
                        2, 2, 2, 2, 1
                };
        uint32_t* indices[PROBLEM_SIZE] =
                {
                        i1, i2, i3, i4, i5
                };
        _Complex float* values[PROBLEM_SIZE] =
                {
                        v1, v2, v3, v4, v5
                };

        for (uint32_t i = 0; i < PROBLEM_SIZE; ++i)
        {
            mtx_res = jmtxc_matrix_crs_build_row(upper, i, counts[i], indices[i], values[i]);
            ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
        }
    }
    print_crsc_matrix(upper);

    MATRIX_TEST_CALL(jmtxc_convert_crs_to_ccs(upper, &cu, NULL));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);

    MATRIX_TEST_CALL(jmtxc_multiply_matrix_crs(lower, cu, &multiplied, NULL));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    MATRIX_TEST_CALL(jmtxcs_matrix_ccs_destroy(cu));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);

    _Complex float exact_multiplied[PROBLEM_SIZE][PROBLEM_SIZE] =
            {
                    {3, 0, 0, 0, 2},
                    {9, 2, 1, 0, 6},
                    {0, 0, 1, 0, 3},
                    {6, 2,-1, 4,-1},
                    {3, 0,-3, 0,-8}
            };
    for (unsigned i = 0; i < PROBLEM_SIZE; ++i)
    {
        for (unsigned j = 0; j < PROBLEM_SIZE; ++j)
        {
            ASSERT(are_close(jmtxc_matrix_crs_get_entry(multiplied, i, j), exact_multiplied[i][j], default_r_tol, default_r_tol));
        }
    }

    print_crsc_matrix(multiplied);

    const _Complex float x_exact[PROBLEM_SIZE] = {1.0f, -2.0f, 3.0f, -4.0f, 5.0f};
    const _Complex float y_exact[PROBLEM_SIZE] = {13, 38, 18,-22,-46};
    _Complex float y[PROBLEM_SIZE];
    _Complex float x[PROBLEM_SIZE];

    jmtxc_matrix_crs_vector_multiply(multiplied, x_exact, y);
    for (unsigned i = 0; i < PROBLEM_SIZE; ++i)
    {
        printf("yex_%u: %g%+g\ty_%u: %g%+g\n", i, crealf(y_exact[i]), cimagf(y_exact[i]), i, crealf(y[i]), cimagf(y[i]));
        ASSERT(are_close(y_exact[i], y[i], default_r_tol, default_a_tol));
    }

    jmtxc_solve_direct_lu_crs(lower, upper, y, x);

    for (unsigned i = 0; i < PROBLEM_SIZE; ++i)
    {
        printf("xex_%u: %g%+g\tx_%u: %g%+g\n", i, crealf(x_exact[i]), cimagf(x_exact[i]), i, crealf(x[i]), cimagf(x[i]));
        ASSERT(are_close(x_exact[i], x[i], default_r_tol, default_a_tol));
    }



    MATRIX_TEST_CALL(jmtxcs_matrix_crs_destroy(multiplied));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    MATRIX_TEST_CALL(jmtxcs_matrix_crs_destroy(combined));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    MATRIX_TEST_CALL(jmtxcs_matrix_crs_destroy(upper));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    MATRIX_TEST_CALL(jmtxcs_matrix_crs_destroy(lower));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);

    return 0;
}