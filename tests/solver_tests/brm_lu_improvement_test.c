#include "../test_common.h"
#include <time.h>
#include <math.h>
#include "decompositions/band_lu_decomposition.h"
#include "solvers/lu_solving.h"
#undef MATRIX_TEST_CALL
#define MATRIX_TEST_CALL(x) ASSERT(x == JMTX_RESULT_SUCCESS)

static JMTX_REAL_T get_rms_error_for_element_count_brm(const unsigned n)
{
    const JMTX_REAL_T dx = 1.0f / (JMTX_REAL_T)n;
    const JMTX_REAL_T r_dx2 = 1.0f / (dx * dx);

    JMTX_NAME_TYPED(matrix_brm) * system_matrix;
    jmtx_result mtx_res;
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

    JMTX_SCALAR_T *const x = calloc(n + 1, sizeof(*x));
    ASSERT(x);
    JMTX_SCALAR_T *const f = calloc(n + 1, sizeof(*f));
    ASSERT(f);
    JMTX_SCALAR_T *const sol = calloc(n + 1, sizeof(*sol));
    ASSERT(sol);
    sol[0] = 0;

    for (uint32_t i = 0; i < n + 1; ++i)
    {
        x[i] = (double)i / (double)n;
        f[i] = sin(x[i] * 1.5f * (double)M_PI);
    }

    JMTX_NAME_TYPED(matrix_brm) * l, *u;

    MATRIX_TEST_CALL(JMTX_NAME_TYPED(decompose_lu_brm)(system_matrix, &l, &u, NULL));
    //    print_brmd_matrix(l);
    //    print_brmd_matrix(u);

    JMTX_NAME_TYPED(solve_direct_lu_brm)(l, u, f + 1, sol + 1);
    JMTX_NAME_TYPED(matrix_brm_destroy)(u);
    (JMTX_NAME_TYPED(matrix_brm_destroy)(l));
    free(f);
    sol[0] = 0.0f;

    enum
    {
        EXACT_POINTS = 257
    };

    JMTX_NAME_TYPED(matrix_brm_destroy)(system_matrix);

    double rms = 0;
    for (uint32_t i = 0; i < n + 1; ++i)
    {
        const JMTX_SCALAR_T err = -sol[i] + (-4.0f * sin(1.5f * (double)M_PI * x[i]) / (9.0f * (double)(M_PI * M_PI)));
        rms += JMTX_DOT(err, err);
    }
    rms = sqrt(rms / (double)(n + 1));

    free(sol);
    free(x);

    return rms;
}

static JMTX_REAL_T get_rms_error_for_element_count_brm_improved(const unsigned n)
{
    const JMTX_REAL_T dx = 1.0f / (JMTX_REAL_T)n;
    const JMTX_REAL_T r_dx2 = 1.0f / (dx * dx);

    JMTX_NAME_TYPED(matrix_brm) * system_matrix;
    jmtx_result mtx_res;
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

    JMTX_SCALAR_T *const x = calloc(n + 1, sizeof(*x));
    ASSERT(x);
    JMTX_SCALAR_T *const f = calloc(n + 1, sizeof(*f));
    ASSERT(f);
    JMTX_SCALAR_T *const sol = calloc(n + 1, sizeof(*sol));
    ASSERT(sol);
    JMTX_SCALAR_T *const aux = calloc(n + 1, sizeof(*aux));
    ASSERT(aux);
    sol[0] = 0;

    for (uint32_t i = 0; i < n + 1; ++i)
    {
        x[i] = (double)i / (double)n;
        f[i] = sin(x[i] * 1.5f * (double)M_PI);
    }

    JMTX_NAME_TYPED(matrix_brm) * l, *u;

    MATRIX_TEST_CALL(JMTX_NAME_TYPED(decompose_lu_brm)(system_matrix, &l, &u, NULL));
    //    print_brmd_matrix(l);
    //    print_brmd_matrix(u);

    JMTX_NAME_TYPED(solver_arguments)
    args = {
        .in_max_iterations = 64, // n,
        .in_convergence_criterion = 1e-8f,
    };
    JMTX_NAME_TYPED(solve_iterative_lu_brm_refine)(system_matrix, l, u, f + 1, sol + 1, aux, &args);
    JMTX_NAME_TYPED(matrix_brm_destroy)(u);
    JMTX_NAME_TYPED(matrix_brm_destroy)(l);
    free(f);
    sol[0] = 0.0f;

    enum
    {
        EXACT_POINTS = 257
    };
    JMTX_NAME_TYPED(matrix_brm_destroy)(system_matrix);

    double rms = 0;
    for (uint32_t i = 0; i < n + 1; ++i)
    {
        const JMTX_REAL_T err = -sol[i] + (-4.0f * sin(1.5f * (double)M_PI * x[i]) / (9.0f * (double)(M_PI * M_PI)));
        rms += JMTX_DOT(err, err);
    }
    rms = sqrt(rms / (double)(n + 1));
    free(aux);
    free(sol);
    free(x);

    return rms;
}

int main(int argc, const char *argv[])
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
        POWER_COUNT = 10
    };
    double errors1[POWER_COUNT];
    double errors2[POWER_COUNT];
    double size[POWER_COUNT];
    unsigned p = 1;

    struct timespec ts0, ts1;
    //    clock_gettime(CLOCK_MONOTONIC, &ts0);
    //    for (unsigned i = 0; i < POWER_COUNT; ++i)
    //    {
    //        p <<= 1;
    //        size[i] = log10f((double)p);
    //        errors[i] = log10f(get_rms_error_for_element_count(p));
    //        printf("RMS error for %u elements was %g\n", p, errors[i]);
    //    }
    //    clock_gettime(CLOCK_MONOTONIC, &ts1);
    //    printf("Time taken using CRS/CCS %g seconds\n", (double)(ts1.tv_sec - ts0.tv_sec) + (double)(ts1.tv_nsec -
    //    ts0.tv_nsec) * 1e-9);

    clock_gettime(CLOCK_MONOTONIC, &ts0);
    p = 1;
    for (unsigned i = 0; i < POWER_COUNT; ++i)
    {
        p <<= 1;
        size[i] = log10f((double)p);
        errors1[i] = log10f(get_rms_error_for_element_count_brm(p));
        errors2[i] = log10f(get_rms_error_for_element_count_brm_improved(p));
        printf("RMS error for %u elements was %g for base and %g for improved\n", p, errors1[i], errors2[i]);
    }
    clock_gettime(CLOCK_MONOTONIC, &ts1);
    printf("Time taken using BRM %g seconds\n",
           (double)(ts1.tv_sec - ts0.tv_sec) + (double)(ts1.tv_nsec - ts0.tv_nsec) * 1e-9);

    return 0;
}
