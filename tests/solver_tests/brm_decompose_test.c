#include "../test_common.h"
#include "decompositions/band_lu_decomposition.h"
#include "solvers/lu_solving.h"
#include <time.h>
#include <math.h>
#undef MATRIX_TEST_CALL
#define MATRIX_TEST_CALL(x) ASSERT(x == JMTX_RESULT_SUCCESS)

static JMTX_SCALAR_T get_rms_error_for_element_count_brm(const unsigned n)
{
    const JMTX_SCALAR_T dx = 1.0f / (JMTX_SCALAR_T)n;
    const JMTX_SCALAR_T r_dx2 = 1.0f / (dx * dx);

    JMTX_NAME_TYPED(matrix_brm) * system_matrix;
    MATRIX_TEST_CALL(JMTX_NAME_TYPED(matrix_brm_new)(&system_matrix, n, n, 1, 1, NULL, NULL));

    //  Build the matrix
    {
        //  First row
        enum
        {
            ROW_SIZE = 2
        };
        JMTX_SCALAR_T val[ROW_SIZE] = {-2.0f * r_dx2, r_dx2};
        JMTX_NAME_TYPED(matrix_brm_set_row)(system_matrix, 0, val);
    }
    {
        //  Interior rows
        enum
        {
            ROW_SIZE = 3
        };
        JMTX_SCALAR_T val[ROW_SIZE] = {r_dx2, -2.0f * r_dx2, r_dx2};
        for (uint32_t row = 1; row < n - 1; ++row)
        {
            JMTX_NAME_TYPED(matrix_brm_set_row)(system_matrix, row, val);
        }
    }
    {
        //  Last row
        enum
        {
            ROW_SIZE = 2
        };
        JMTX_SCALAR_T val[ROW_SIZE] = {2.0f * r_dx2, -2.0f * r_dx2};
        JMTX_NAME_TYPED(matrix_brm_set_row)(system_matrix, n - 1, val);
    }

    JMTX_SCALAR_T *const f = calloc(n + 1, sizeof(*f));
    ASSERT(f);
    JMTX_SCALAR_T *const sol = calloc(n + 1, sizeof(*sol));
    ASSERT(sol);
    sol[0] = 0;

    for (uint32_t i = 0; i < n + 1; ++i)
    {
        const double x = (double)i / (double)n;
        f[i] = sin(x * 1.5f * M_PI);
    }

    JMTX_NAME_TYPED(matrix_brm) * l, *u;

    MATRIX_TEST_CALL(JMTX_NAME_TYPED(decompose_lu_brm)(system_matrix, &l, &u, NULL));
    //    print_brmd_matrix(l);
    //    print_brmd_matrix(u);
    JMTX_NAME_TYPED(solve_direct_lu_brm)(l, u, f + 1, sol + 1);
    JMTX_NAME_TYPED(matrix_brm_destroy)(u);
    JMTX_NAME_TYPED(matrix_brm_destroy)(l);
    free(f);
    sol[0] = 0.0f;

    enum
    {
        EXACT_POINTS = 257
    };
    JMTX_SCALAR_T *const err = calloc(n + 1, sizeof(*err));
    ASSERT(err);
    JMTX_NAME_TYPED(matrix_brm_destroy)(system_matrix);

    JMTX_SCALAR_T rms = 0;
    for (uint32_t i = 0; i < n + 1; ++i)
    {
        const double x = (double)i / (double)n;
        err[i] = -sol[i] - 4.0f * sin(1.5f * M_PI * x) / (9.0f * (M_PI * M_PI));
        rms += err[i] * err[i];
    }
    rms = sqrt(rms / (JMTX_SCALAR_T)(n + 1));
    free(err);
    free(sol);
    return rms;
}

int main()
{
    //  Copied after work of Marc Gerritsma

    /*
     * Sample problem:
     *      d^2 u / dx^2 = f(x)
     *      u(0) = 0
     *      d u / dx (1) = 0
     * uniform mesh spacing on domain (0, 1)
     */

    enum
    {
        POWER_COUNT = 16
    };
    double errors[POWER_COUNT];
    // double size[POWER_COUNT];
    unsigned p = 1;

    struct timespec ts0, ts1;
    //    clock_gettime(CLOCK_MONOTONIC, &ts0);
    //    for (unsigned i = 0; i < POWER_COUNT; ++i)
    //    {
    //        p <<= 1;
    //        size[i] = log10f((JMTX_SCALAR_T)p);
    //        errors[i] = log10f(get_rms_error_for_element_count(p));
    //        printf("RMS error for %u elements was %g\n", p, errors[i]);
    //    }
    //    clock_gettime(CLOCK_MONOTONIC, &ts1);
    //    printf("Time taken using CRS/CCS %g seconds\n", (JMTX_SCALAR_T)(ts1.tv_sec - ts0.tv_sec) +
    //    (JMTX_SCALAR_T)(ts1.tv_nsec - ts0.tv_nsec) * 1e-9);

    clock_gettime(CLOCK_MONOTONIC, &ts0);
    for (unsigned i = 0; i < POWER_COUNT; ++i)
    {
        p <<= 1;
        // size[i] = log10(p);
        errors[i] = log10(get_rms_error_for_element_count_brm(p));
        printf("RMS error for %u elements was %g\n", p, errors[i]);
    }
    clock_gettime(CLOCK_MONOTONIC, &ts1);
    printf("Time taken using BRM %g seconds\n",
           (double)(ts1.tv_sec - ts0.tv_sec) + (double)(ts1.tv_nsec - ts0.tv_nsec) * 1e-9);

    return 0;
}
