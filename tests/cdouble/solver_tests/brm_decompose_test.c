// Automatically generated from tests/cfloat/solver_tests/brm_decompose_test.c on Fri Dec  1 18:48:10 2023
// Automatically generated from tests/cdouble/solver_tests/brm_decompose_test.c on Fri Dec  1 17:35:45 2023
#include "../../../include/jmtx/cdouble/matrices/sparse_row_compressed_safe.h"
#include "../../../include/jmtx/cdouble/matrices/sparse_column_compressed_safe.h"
#include "../../../include/jmtx/cdouble/matrices/band_row_major_safe.h"
#include "../../../include/jmtx/cdouble/solvers/incomplete_lu_decomposition.h"
#include "../../../include/jmtx/cdouble/solvers/band_lu_decomposition.h"
#include "../../../include/jmtx/cdouble/matrices/sparse_conversion.h"
#include "../../../include/jmtx/cdouble/solvers/lu_solving.h"
#include "../test_common.h"
#include <time.h>
#include <math.h>
#include <complex.h>
#undef MATRIX_TEST_CALL
#define MATRIX_TEST_CALL(x) ASSERT(x == JMTX_RESULT_SUCCESS)

static _Complex double get_rms_error_for_element_count(const unsigned n)
{
    const double dx = 1.0f / (_Complex double)n;
    const double r_dx2 = 1.0f / (dx * dx);

    jmtxz_matrix_crs* system_matrix;
    MATRIX_TEST_CALL(jmtxzs_matrix_crs_new(&system_matrix, n, n, 3 * n < n * n ? 3 * n : 0, NULL));

    //  Build the matrix
    {
        //  First row
        enum {ROW_SIZE = 2};
        const uint32_t idx[ROW_SIZE] = {0, 1};
        const _Complex double val[ROW_SIZE] = {-2.0f * r_dx2, r_dx2};
        MATRIX_TEST_CALL(jmtxz_matrix_crs_build_row(system_matrix, 0, ROW_SIZE, idx, val));
    }
    {
        //  Interior rows
        enum {ROW_SIZE = 3};
        const _Complex double val[ROW_SIZE] = {r_dx2, -2.0f * r_dx2, r_dx2};
        uint32_t idx[ROW_SIZE];
        for (uint32_t row = 1; row < n - 1; ++row)
        {
            idx[0] = row - 1;
            idx[1] = row;
            idx[2] = row + 1;
            MATRIX_TEST_CALL(jmtxz_matrix_crs_build_row(system_matrix, row, ROW_SIZE, idx, val));
        }
    }
    {
        //  Last row
        enum {ROW_SIZE = 2};
        const uint32_t idx[ROW_SIZE] = {n - 2, n - 1};
        const _Complex double val[ROW_SIZE] = {2.0f * r_dx2, -2.0f * r_dx2};
        MATRIX_TEST_CALL(jmtxz_matrix_crs_build_row(system_matrix, n - 1, ROW_SIZE, idx, val));
    }

    _Complex double* const x = calloc(n + 1, sizeof(*x));
    ASSERT(x);
    _Complex double* const f = calloc(n + 1, sizeof(*f));
    ASSERT(f);
    _Complex double* const sol = calloc(n + 1, sizeof(*sol));
    ASSERT(sol);
    sol[0] = 0;

    for (uint32_t i = 0; i < n + 1; ++i)
    {
        x[i] = (_Complex double)i / (_Complex double)n;
        f[i] = sinf(x[i] * 1.5f * (_Complex double)M_PI);
    }

    jmtxz_matrix_crs* l, *u;
    jmtxz_matrix_ccs* cu;

    MATRIX_TEST_CALL(jmtxz_incomplete_lu_crs(system_matrix, &l, &cu, NULL));

    MATRIX_TEST_CALL(jmtxz_convert_ccs_to_crs(cu, &u, NULL));
    MATRIX_TEST_CALL(jmtxzs_matrix_ccs_destroy(cu));
//    print_crs_matrix(l);
//    print_crs_matrix(u);

    jmtxz_lu_solve_crs(l, u, f + 1, sol + 1);
    MATRIX_TEST_CALL(jmtxzs_matrix_crs_destroy(u));
    MATRIX_TEST_CALL(jmtxzs_matrix_crs_destroy(l));
    free(f);
    sol[0] = 0.0f;

    enum {EXACT_POINTS = 257};
    _Complex double* const err = calloc(n + 1, sizeof(*err));
    ASSERT(err);
    MATRIX_TEST_CALL(jmtxzs_matrix_crs_destroy(system_matrix));

    double rms = 0;
    for (uint32_t i = 0; i < n + 1; ++i)
    {
        err[i] = -sol[i] + (-4.0f * sinf(1.5f * (_Complex double)M_PI * x[i]) / (9.0f * (_Complex double)(M_PI * M_PI)));
        rms += conjf(err[i]) * err[i];
    }
    rms = sqrt(rms / (double)(n + 1));
    free(err);
    free(sol);
    free(x);

    return rms;
}


static _Complex double get_rms_error_for_element_count_brm(const unsigned n)
{
    const _Complex double dx = 1.0f / (_Complex double)n;
    const _Complex double r_dx2 = 1.0f / (dx * dx);

    jmtxz_matrix_brm* system_matrix;
    MATRIX_TEST_CALL(jmtxzs_matrix_brm_new(&system_matrix, n, n, 1, 1, NULL, NULL));

    //  Build the matrix
    {
        //  First row
        enum {ROW_SIZE = 2};
        _Complex double val[ROW_SIZE] = {-2.0f * r_dx2, r_dx2};
        jmtxz_matrix_brm_set_row(system_matrix, 0, val);
    }
    {
        //  Interior rows
        enum {ROW_SIZE = 3};
        _Complex double val[ROW_SIZE] = {r_dx2, -2.0f * r_dx2, r_dx2};
        for (uint32_t row = 1; row < n - 1; ++row)
        {
            jmtxz_matrix_brm_set_row(system_matrix, row, val);
        }
    }
    {
        //  Last row
        enum {ROW_SIZE = 2};
        _Complex double val[ROW_SIZE] = {2.0f * r_dx2, -2.0f * r_dx2};
        jmtxz_matrix_brm_set_row(system_matrix, n - 1, val);
    }

    _Complex double* const x = calloc(n + 1, sizeof(*x));
    ASSERT(x);
    _Complex double* const f = calloc(n + 1, sizeof(*f));
    ASSERT(f);
    _Complex double* const sol = calloc(n + 1, sizeof(*sol));
    ASSERT(sol);
    sol[0] = 0;

    for (uint32_t i = 0; i < n + 1; ++i)
    {
        x[i] = (_Complex double)i / (_Complex double)n;
        f[i] = sinf(x[i] * 1.5f * (_Complex double)M_PI);
    }

    jmtxz_matrix_brm* l, *u;

    MATRIX_TEST_CALL(jmtxz_band_lu_decomposition_brm(system_matrix, &l, &u, NULL));
//    print_brm_matrix(l);
//    print_brm_matrix(u);

    jmtxz_lu_solve_brm(l, u, f + 1, sol + 1);
    MATRIX_TEST_CALL(jmtxzs_matrix_brm_destroy(u));
    MATRIX_TEST_CALL(jmtxzs_matrix_brm_destroy(l));
    free(f);
    sol[0] = 0.0f;

    enum {EXACT_POINTS = 257};
    _Complex double* const err = calloc(n + 1, sizeof(*err));
    ASSERT(err);
    MATRIX_TEST_CALL(jmtxzs_matrix_brm_destroy(system_matrix));

    _Complex double rms = 0;
    for (uint32_t i = 0; i < n + 1; ++i)
    {
        err[i] = -sol[i] + (-4.0f * sinf(1.5f * (_Complex double)M_PI * x[i]) / (9.0f * (_Complex double)(M_PI * M_PI)));
        rms += err[i] * err[i];
    }
    rms = sqrt(rms / (_Complex double)(n + 1));
    free(err);
    free(sol);
    free(x);

    return rms;
}

int main(int argc, const char* argv[static argc])
{
    (void)argc; (void)argv;
    //  Copied after work of Marc Gerritsma

    /*
     * Sample problem:
     *      d^2 u / dx^2 = f(x)
     *      u(0) = 0
     *      d u / dx (1) = 0
     * uniform mesh spacing on domain (0, 1)
     */

    enum {POWER_COUNT = 16};
    double errors[POWER_COUNT];
    double size[POWER_COUNT];
    unsigned p = 1;

    struct timespec ts0, ts1;
//    clock_gettime(CLOCK_MONOTONIC, &ts0);
//    for (unsigned i = 0; i < POWER_COUNT; ++i)
//    {
//        p <<= 1;
//        size[i] = log10f((_Complex double)p);
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
        size[i] = log10f((double)p);
        errors[i] = log10f(get_rms_error_for_element_count_brm(p));
        printf("RMS error for %u elements was %g\n", p, errors[i]);
    }
    clock_gettime(CLOCK_MONOTONIC, &ts1);
    printf("Time taken using BRM %g seconds\n", (double)(ts1.tv_sec - ts0.tv_sec) + (double)(ts1.tv_nsec - ts0.tv_nsec) * 1e-9);


    return 0;
}
