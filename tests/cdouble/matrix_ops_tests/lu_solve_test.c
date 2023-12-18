// Automatically generated from tests/cfloat/matrix_ops_tests/lu_solve_test.c on Fri Dec  1 18:48:10 2023
// Automatically generated from tests/cdouble/matrix_ops_tests/lu_solve_test.c on Fri Dec  1 17:35:45 2023
//
// Created by jan on 23.11.2023.
//
#include "../test_common.h"
#include "../../../include/jmtx/cdouble/matrices/sparse_row_compressed_safe.h"
#include "../../../include/jmtx/cdouble/matrices/sparse_column_compressed_safe.h"
#include "../../../include/jmtx/cdouble/matrices/sparse_conversion.h"
#include "../../../include/jmtx/cdouble/matrices/sparse_multiplication.h"
#include "../../../include/jmtx/cdouble/solvers/lu_solving.h"
#include <float.h>
#include <math.h>
#include <assert.h>
#include <complex.h>
int are_close(_Complex double v1, _Complex double v2, double relative_tol, double abs_tol)
{
    if (v1 == v2)
    {
        return 1;
    }
    assert((abs_tol) >= 0.0f);
    assert((relative_tol) >= 0.0f);
    const _Complex double dif = (v1 - v2);
    if (cabs(dif) < abs_tol)
    {
        return 1;
    }
    if (cabs(v1) > cabs(v2))
    {
        if (cabs(v2 + dif) > cabs(v1))
        {
            return 1;
        }
    }
    else
    {
        if (cabs(v1 + dif) > cabs(v2))
        {
            return 1;
        }
    }
    return 0;
}

static const _Complex double default_r_tol = FLT_EPSILON * 400;
static const _Complex double default_a_tol = FLT_EPSILON * 100;

enum {PROBLEM_SIZE = 5};

int main()
{
    jmtxz_matrix_crs* lower, *upper, *combined, *multiplied;
    jmtxz_matrix_ccs* cu;
    jmtx_result mtx_res;

    MATRIX_TEST_CALL(jmtxzs_matrix_crs_new(&lower, PROBLEM_SIZE, PROBLEM_SIZE, 0, NULL));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    MATRIX_TEST_CALL(jmtxzs_matrix_crs_new(&upper, PROBLEM_SIZE, PROBLEM_SIZE, 0, NULL));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    MATRIX_TEST_CALL(jmtxzs_matrix_crs_new(&combined, PROBLEM_SIZE, PROBLEM_SIZE, 0, NULL));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);

    //  Make L based on predefined values
    {
        _Complex double v1[1] = {1.0f};
        uint32_t i1[1] = {0};


        _Complex double v2[2] = {3.0f, 1.0f};
        uint32_t i2[2] = {0, 1};

        _Complex double v3[1] = {1.0f};
        uint32_t i3[1] = {2};

        _Complex double v4[4] = {2.0f, 1.0f, -2.0f, 1.0f};
        uint32_t i4[4] = {0, 1, 2, 3};

        _Complex double v5[3] = {1.0f, -3.0f, 1.0f};
        uint32_t i5[3] = {0, 2, 4};

        uint32_t counts[PROBLEM_SIZE] =
                {
                1, 2, 1, 4, 3
                };
        uint32_t* indices[PROBLEM_SIZE] =
                {
                i1, i2, i3, i4, i5
                };
        _Complex double* values[PROBLEM_SIZE] =
                {
                v1, v2, v3, v4, v5
                };

        for (uint32_t i = 0; i < PROBLEM_SIZE; ++i)
        {
            mtx_res = jmtxz_matrix_crs_build_row(lower, i, counts[i], indices[i], values[i]);
            ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
        }
    }
    print_crs_matrix(lower);

    //  Make U based on predefined values
    {
        _Complex double v1[2] = {3.0f, 2.0f};
        uint32_t i1[2] = {0, 4};

        _Complex double v2[2] = {2.0f, 1.0f};
        uint32_t i2[2] = {1, 2};

        _Complex double v3[2] = {1.0f, 3.0f};
        uint32_t i3[2] = {2, 4};

        _Complex double v4[2] = {4.0f, 1.0f};
        uint32_t i4[2] = {3, 4};

        _Complex double v5[1] = {-1.0f};
        uint32_t i5[1] = {4};

        uint32_t counts[PROBLEM_SIZE] =
                {
                        2, 2, 2, 2, 1
                };
        uint32_t* indices[PROBLEM_SIZE] =
                {
                        i1, i2, i3, i4, i5
                };
        _Complex double* values[PROBLEM_SIZE] =
                {
                        v1, v2, v3, v4, v5
                };

        for (uint32_t i = 0; i < PROBLEM_SIZE; ++i)
        {
            mtx_res = jmtxz_matrix_crs_build_row(upper, i, counts[i], indices[i], values[i]);
            ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
        }
    }
    print_crs_matrix(upper);

    MATRIX_TEST_CALL(jmtxz_convert_crs_to_ccs(upper, &cu, NULL));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);

    MATRIX_TEST_CALL(jmtxz_multiply_matrix_crs(lower, cu, &multiplied, NULL));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    MATRIX_TEST_CALL(jmtxzs_matrix_ccs_destroy(cu));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);

    _Complex double exact_multiplied[PROBLEM_SIZE][PROBLEM_SIZE] =
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
            ASSERT(are_close(jmtxz_matrix_crs_get_entry(multiplied, i, j), exact_multiplied[i][j], default_r_tol, default_r_tol));
        }
    }

    print_crs_matrix(multiplied);

    const _Complex double x_exact[PROBLEM_SIZE] = {1.0f, -2.0f, 3.0f, -4.0f, 5.0f};
    const _Complex double y_exact[PROBLEM_SIZE] = {13, 38, 18,-22,-46};
    _Complex double y[PROBLEM_SIZE];
    _Complex double x[PROBLEM_SIZE];

    jmtxz_matrix_crs_vector_multiply(multiplied, x_exact, y);
    for (unsigned i = 0; i < PROBLEM_SIZE; ++i)
    {
        printf("yex_%u: %g%+g\ty_%u: %g%+g\n", i, creal(y_exact[i]), cimag(y_exact[i]), i, creal(y[i]), cimag(y[i]));
        ASSERT(are_close(y_exact[i], y[i], default_r_tol, default_a_tol));
    }

    jmtxz_solve_direct_lu_crs(lower, upper, y, x);

    for (unsigned i = 0; i < PROBLEM_SIZE; ++i)
    {
        printf("xex_%u: %g%+g\tx_%u: %g%+g\n", i, creal(x_exact[i]), cimag(x_exact[i]), i, creal(x[i]), cimag(x[i]));
        ASSERT(are_close(x_exact[i], x[i], default_r_tol, default_a_tol));
    }



    MATRIX_TEST_CALL(jmtxzs_matrix_crs_destroy(multiplied));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    MATRIX_TEST_CALL(jmtxzs_matrix_crs_destroy(combined));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    MATRIX_TEST_CALL(jmtxzs_matrix_crs_destroy(upper));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    MATRIX_TEST_CALL(jmtxzs_matrix_crs_destroy(lower));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);

    return 0;
}