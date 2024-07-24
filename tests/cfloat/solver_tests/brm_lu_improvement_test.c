// Automatically generated from tests/float/solver_tests/brm_lu_improvement_test.c on Fri Dec  1 17:35:45 2023
#include "../../../include/jmtx/cfloat/matrices/sparse_row_compressed_safe.h"
#include "../../../include/jmtx/cfloat/matrices/band_row_major_safe.h"
#include "../../../include/jmtx/cfloat/decompositions/incomplete_lu_decomposition.h"
#include "../../../include/jmtx/cfloat/decompositions/band_lu_decomposition.h"
#include "../../../include/jmtx/cfloat/solvers/lu_solving.h"
#include "../test_common.h"
#include <time.h>
#include <math.h>
#include <complex.h>
#undef MATRIX_TEST_CALL
#define MATRIX_TEST_CALL(x) ASSERT(x == JMTX_RESULT_SUCCESS)

static float get_rms_error_for_element_count_brm(const unsigned n)
{
    const float dx = 1.0f / (_Complex float)n;
    const float r_dx2 = 1.0f / (dx * dx);

    jmtxc_matrix_brm* system_matrix;
    MATRIX_TEST_CALL(jmtxcs_matrix_brm_new(&system_matrix, n, n, 1, 1, NULL, NULL));

    //  Build the matrix
    {
        //  First row
        enum {ROW_SIZE = 2};
        _Complex float val[ROW_SIZE] = {-2.0f * r_dx2, r_dx2};
        jmtxc_matrix_brm_set_row(system_matrix, 0, val);
    }
    {
        //  Interior rows
        enum {ROW_SIZE = 3};
        _Complex float val[ROW_SIZE] = {r_dx2, -2.0f * r_dx2, r_dx2};
        for (uint32_t row = 1; row < n - 1; ++row)
        {
            jmtxc_matrix_brm_set_row(system_matrix, row, val);
        }
    }
    {
        //  Last row
        enum {ROW_SIZE = 2};
        _Complex float val[ROW_SIZE] = {2.0f * r_dx2, -2.0f * r_dx2};
        jmtxc_matrix_brm_set_row(system_matrix, n - 1, val);
    }

    _Complex float* const x = calloc(n + 1, sizeof(*x));
    ASSERT(x);
    _Complex float* const f = calloc(n + 1, sizeof(*f));
    ASSERT(f);
    _Complex float* const sol = calloc(n + 1, sizeof(*sol));
    ASSERT(sol);
    sol[0] = 0;

    for (uint32_t i = 0; i < n + 1; ++i)
    {
        x[i] = (_Complex float)i / (_Complex float)n;
        f[i] = sinf(x[i] * 1.5f * (_Complex float)M_PI);
    }

    jmtxc_matrix_brm* l, *u;

    MATRIX_TEST_CALL(jmtxc_decompose_lu_brm(system_matrix, &l, &u, NULL));
//    print_brmc_matrix(l);
//    print_brmc_matrix(u);

    jmtxc_solve_direct_lu_brm(l, u, f + 1, sol + 1);
    MATRIX_TEST_CALL(jmtxcs_matrix_brm_destroy(u));
    MATRIX_TEST_CALL(jmtxcs_matrix_brm_destroy(l));
    free(f);
    sol[0] = 0.0f;

    enum {EXACT_POINTS = 257};
    _Complex float* const err = calloc(n + 1, sizeof(*err));
    ASSERT(err);
    MATRIX_TEST_CALL(jmtxcs_matrix_brm_destroy(system_matrix));

    float rms = 0;
    for (uint32_t i = 0; i < n + 1; ++i)
    {
        err[i] = -sol[i] + (-4.0f * sinf(1.5f * (_Complex float)M_PI * x[i]) / (9.0f * (_Complex float)(M_PI * M_PI)));
        rms += conjf(err[i]) * err[i];
    }
    rms = sqrtf(rms / (float)(n + 1));
    free(err);
    free(sol);
    free(x);

    return rms;
}

static float get_rms_error_for_element_count_brm_improved(const unsigned n)
{
    const _Complex float dx = 1.0f / (_Complex float)n;
    const _Complex float r_dx2 = 1.0f / (dx * dx);

    jmtxc_matrix_brm* system_matrix;
    MATRIX_TEST_CALL(jmtxcs_matrix_brm_new(&system_matrix, n, n, 1, 1, NULL, NULL));

    //  Build the matrix
    {
        //  First row
        enum {ROW_SIZE = 2};
        _Complex float val[ROW_SIZE] = {-2.0f * r_dx2, r_dx2};
        jmtxc_matrix_brm_set_row(system_matrix, 0, val);
    }
    {
        //  Interior rows
        enum {ROW_SIZE = 3};
        _Complex float val[ROW_SIZE] = {r_dx2, -2.0f * r_dx2, r_dx2};
        for (uint32_t row = 1; row < n - 1; ++row)
        {
            jmtxc_matrix_brm_set_row(system_matrix, row, val);
        }
    }
    {
        //  Last row
        enum {ROW_SIZE = 2};
        _Complex float val[ROW_SIZE] = {2.0f * r_dx2, -2.0f * r_dx2};
        jmtxc_matrix_brm_set_row(system_matrix, n - 1, val);
    }

    _Complex float* const x = calloc(n + 1, sizeof(*x));
    ASSERT(x);
    _Complex float* const f = calloc(n + 1, sizeof(*f));
    ASSERT(f);
    _Complex float* const sol = calloc(n + 1, sizeof(*sol));
    ASSERT(sol);
    _Complex float* const aux = calloc(n + 1, sizeof(*aux));
    ASSERT(aux);
    sol[0] = 0;

    for (uint32_t i = 0; i < n + 1; ++i)
    {
        x[i] = (_Complex float)i / (_Complex float)n;
        f[i] = sinf(x[i] * 1.5f * (_Complex float)M_PI);
    }

    jmtxc_matrix_brm* l, *u;

    MATRIX_TEST_CALL(jmtxc_decompose_lu_brm(system_matrix, &l, &u, NULL));
//    print_brmc_matrix(l);
//    print_brmc_matrix(u);

    jmtx_solver_arguments args =
            {
            .in_max_iterations = 64,//n,
            .in_convergence_criterion = 1e-8f,
            };
    jmtxc_solve_iterative_lu_brm_refine(system_matrix, l, u, f + 1, sol + 1, aux, &args);
    MATRIX_TEST_CALL(jmtxcs_matrix_brm_destroy(u));
    MATRIX_TEST_CALL(jmtxcs_matrix_brm_destroy(l));
    free(f);
    sol[0] = 0.0f;

    enum {EXACT_POINTS = 257};
    _Complex float* const err = calloc(n + 1, sizeof(*err));
    ASSERT(err);
    MATRIX_TEST_CALL(jmtxcs_matrix_brm_destroy(system_matrix));

    float rms = 0;
    for (uint32_t i = 0; i < n + 1; ++i)
    {
        err[i] = -sol[i] + (-4.0f * sinf(1.5f * (_Complex float)M_PI * x[i]) / (9.0f * (_Complex float)(M_PI * M_PI)));
        rms += conjf(err[i]) * err[i];
    }
    rms = sqrtf(rms / (float)(n + 1));
    free(aux);
    free(err);
    free(sol);
    free(x);

    return rms;
}

int main(int argc, const char* argv[])
{
    (void)argv; (void)argc;
    //  Copied after work of Marc Gerritsma

    /*
     * Sample problem:
     *      d^2 u / dx^2 = f(x)
     *      u(0) = 0
     *      d u / dx (1) = 0
     * uniform mesh spacing on domain (0, 1)
     */

    enum {POWER_COUNT = 10};
    float errors1[POWER_COUNT];
    float errors2[POWER_COUNT];
    float size[POWER_COUNT];
    unsigned p = 1;

    struct timespec ts0, ts1;
//    clock_gettime(CLOCK_MONOTONIC, &ts0);
//    for (unsigned i = 0; i < POWER_COUNT; ++i)
//    {
//        p <<= 1;
//        size[i] = log10f((_Complex float)p);
//        errors[i] = log10f(get_rms_error_for_element_count(p));
//        printf("RMS error for %u elements was %g\n", p, errors[i]);
//    }
//    clock_gettime(CLOCK_MONOTONIC, &ts1);
//    printf("Time taken using CRS/CCS %g seconds\n", (double)(ts1.tv_sec - ts0.tv_sec) + (double)(ts1.tv_nsec - ts0.tv_nsec) * 1e-9);
//    p = 1;

    clock_gettime(CLOCK_MONOTONIC, &ts0);
    for (unsigned i = 0; i < POWER_COUNT; ++i)
    {
        p <<= 1;
        size[i] = log10f((float)p);
        errors1[i] = log10f(get_rms_error_for_element_count_brm(p));
        errors2[i] = log10f(get_rms_error_for_element_count_brm_improved(p));
        printf("RMS error for %u elements was %g for base and %g for improved\n", p, errors1[i], errors2[i]);
    }
    clock_gettime(CLOCK_MONOTONIC, &ts1);
    printf("Time taken using BRM %g seconds\n", (double)(ts1.tv_sec - ts0.tv_sec) + (double)(ts1.tv_nsec - ts0.tv_nsec) * 1e-9);


    return 0;
}
