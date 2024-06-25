#include "../../../include/jmtx/float/matrices/sparse_row_compressed_safe.h"
#include "../../../include/jmtx/float/matrices/sparse_column_compressed_safe.h"
#include "../../../include/jmtx/float/matrices/band_row_major_safe.h"
#include "../../../include/jmtx/float/decompositions/incomplete_lu_decomposition.h"
#include "../../../include/jmtx/float/decompositions/band_lu_decomposition.h"
#include "../../../include/jmtx/float/matrices/sparse_conversion.h"
#include "../../../include/jmtx/float/matrices/sparse_multiplication.h"
#include "../../../include/jmtx/float/solvers/lu_solving.h"
#include "../test_common.h"
#include <time.h>
#include <math.h>
#undef MATRIX_TEST_CALL
#define MATRIX_TEST_CALL(x) ASSERT(x == JMTX_RESULT_SUCCESS)

static float get_rms_error_for_element_count(const unsigned n)
{
    const float dx = 1.0f / (float)n;
    const float r_dx2 = 1.0f / (dx * dx);

    jmtx_matrix_crs* system_matrix;
    jmtx_result mtx_res;
    MATRIX_TEST_CALL(jmtxs_matrix_crs_new(&system_matrix, n, n, 3 * n < n * n ? 3 * n : 0, NULL));

    //  Build the matrix
    {
        //  First row
        enum {ROW_SIZE = 2};
        const uint32_t idx[ROW_SIZE] = {0, 1};
        const float val[ROW_SIZE] = {-2.0f * r_dx2, r_dx2};
        MATRIX_TEST_CALL(jmtx_matrix_crs_build_row(system_matrix, 0, ROW_SIZE, idx, val));
    }
    {
        //  Interior rows
        enum {ROW_SIZE = 3};
        const float val[ROW_SIZE] = {r_dx2, -2.0f * r_dx2, r_dx2};
        uint32_t idx[ROW_SIZE];
        for (uint32_t row = 1; row < n - 1; ++row)
        {
            idx[0] = row - 1;
            idx[1] = row;
            idx[2] = row + 1;
            MATRIX_TEST_CALL(jmtx_matrix_crs_build_row(system_matrix, row, ROW_SIZE, idx, val));
        }
    }
    {
        //  Last row
        enum {ROW_SIZE = 2};
        const uint32_t idx[ROW_SIZE] = {n - 2, n - 1};
        const float val[ROW_SIZE] = {2.0f * r_dx2, -2.0f * r_dx2};
        MATRIX_TEST_CALL(jmtx_matrix_crs_build_row(system_matrix, n - 1, ROW_SIZE, idx, val));
    }

    float* const x = calloc(n + 1, sizeof(*x));
    ASSERT(x);
    float* const f = calloc(n + 1, sizeof(*f));
    ASSERT(f);
    float* const sol = calloc(n + 1, sizeof(*sol));
    ASSERT(sol);
    sol[0] = 0;

    for (uint32_t i = 0; i < n + 1; ++i)
    {
        x[i] = (float)i / (float)n;
        f[i] = sinf(x[i] * 1.5f * (float)M_PI);
    }

    jmtx_matrix_crs* l, *u;
    jmtx_matrix_ccs* cu;

    MATRIX_TEST_CALL(jmtx_decompose_ilu_crs(system_matrix, &l, &cu, NULL));

    MATRIX_TEST_CALL(jmtx_convert_ccs_to_crs(cu, &u, NULL));
    MATRIX_TEST_CALL(jmtxs_matrix_ccs_destroy(cu));
//    print_crs_matrix(l);
//    print_crs_matrix(u);

    jmtx_solve_direct_lu_crs(l, u, f + 1, sol + 1);
    MATRIX_TEST_CALL(jmtxs_matrix_crs_destroy(u));
    MATRIX_TEST_CALL(jmtxs_matrix_crs_destroy(l));
    free(f);
    sol[0] = 0.0f;

    enum {EXACT_POINTS = 257};
    float* const err = calloc(n + 1, sizeof(*err));
    ASSERT(err);
    MATRIX_TEST_CALL(jmtxs_matrix_crs_destroy(system_matrix));

    float rms = 0;
    for (uint32_t i = 0; i < n + 1; ++i)
    {
        err[i] = -sol[i] + (-4.0f * sinf(1.5f * (float)M_PI * x[i]) / (9.0f * (float)(M_PI * M_PI)));
        rms += err[i] * err[i];
    }
    rms = sqrtf(rms / (float)(n + 1));
    free(err);
    free(sol);
    free(x);

    return rms;
}


static float get_rms_error_for_element_count_brm(const unsigned n)
{
    const float dx = 1.0f / (float)n;
    const float r_dx2 = 1.0f / (dx * dx);

    jmtx_matrix_brm* system_matrix;
    jmtx_result mtx_res;
    MATRIX_TEST_CALL(jmtxs_matrix_brm_new(&system_matrix, n, n, 1, 1, NULL, NULL));

    //  Build the matrix
    {
        //  First row
        enum {ROW_SIZE = 2};
        float val[ROW_SIZE] = {-2.0f * r_dx2, r_dx2};
        jmtx_matrix_brm_set_row(system_matrix, 0, val);
    }
    {
        //  Interior rows
        enum {ROW_SIZE = 3};
        float val[ROW_SIZE] = {r_dx2, -2.0f * r_dx2, r_dx2};
        for (uint32_t row = 1; row < n - 1; ++row)
        {
            jmtx_matrix_brm_set_row(system_matrix, row, val);
        }
    }
    {
        //  Last row
        enum {ROW_SIZE = 2};
        float val[ROW_SIZE] = {2.0f * r_dx2, -2.0f * r_dx2};
        jmtx_matrix_brm_set_row(system_matrix, n - 1, val);
    }

    float* const x = calloc(n + 1, sizeof(*x));
    ASSERT(x);
    float* const f = calloc(n + 1, sizeof(*f));
    ASSERT(f);
    float* const sol = calloc(n + 1, sizeof(*sol));
    ASSERT(sol);
    sol[0] = 0;

    for (uint32_t i = 0; i < n + 1; ++i)
    {
        x[i] = (float)i / (float)n;
        f[i] = sinf(x[i] * 1.5f * (float)M_PI);
    }

    jmtx_matrix_brm* l, *u;

    MATRIX_TEST_CALL(jmtx_decompose_lu_brm(system_matrix, &l, &u, NULL));
//    print_brm_matrix(l);
//    print_brm_matrix(u);

    jmtx_solve_direct_lu_brm(l, u, f + 1, sol + 1);
    MATRIX_TEST_CALL(jmtxs_matrix_brm_destroy(u));
    MATRIX_TEST_CALL(jmtxs_matrix_brm_destroy(l));
    free(f);
    sol[0] = 0.0f;

    enum {EXACT_POINTS = 257};
    float* const err = calloc(n + 1, sizeof(*err));
    ASSERT(err);
    MATRIX_TEST_CALL(jmtxs_matrix_brm_destroy(system_matrix));

    float rms = 0;
    for (uint32_t i = 0; i < n + 1; ++i)
    {
        err[i] = -sol[i] + (-4.0f * sinf(1.5f * (float)M_PI * x[i]) / (9.0f * (float)(M_PI * M_PI)));
        rms += err[i] * err[i];
    }
    rms = sqrtf(rms / (float)(n + 1));
    free(err);
    free(sol);
    free(x);

    return rms;
}

int main(int argc, const char* argv[static argc])
{
    //  Copied after work of Marc Gerritsma

    /*
     * Sample problem:
     *      d^2 u / dx^2 = f(x)
     *      u(0) = 0
     *      d u / dx (1) = 0
     * uniform mesh spacing on domain (0, 1)
     */

    enum {POWER_COUNT = 16};
    float errors[POWER_COUNT];
    float size[POWER_COUNT];
    unsigned p = 1;

    struct timespec ts0, ts1;
//    clock_gettime(CLOCK_MONOTONIC, &ts0);
//    for (unsigned i = 0; i < POWER_COUNT; ++i)
//    {
//        p <<= 1;
//        size[i] = log10f((float)p);
//        errors[i] = log10f(get_rms_error_for_element_count(p));
//        printf("RMS error for %u elements was %g\n", p, errors[i]);
//    }
//    clock_gettime(CLOCK_MONOTONIC, &ts1);
//    printf("Time taken using CRS/CCS %g seconds\n", (double)(ts1.tv_sec - ts0.tv_sec) + (double)(ts1.tv_nsec - ts0.tv_nsec) * 1e-9);

    clock_gettime(CLOCK_MONOTONIC, &ts0);
    p = 1;
    for (unsigned i = 0; i < POWER_COUNT; ++i)
    {
        p <<= 1;
        size[i] = log10f((float)p);
        errors[i] = log10f(get_rms_error_for_element_count_brm(p));
        printf("RMS error for %u elements was %g\n", p, errors[i]);
    }
    clock_gettime(CLOCK_MONOTONIC, &ts1);
    printf("Time taken using BRM %g seconds\n", (double)(ts1.tv_sec - ts0.tv_sec) + (double)(ts1.tv_nsec - ts0.tv_nsec) * 1e-9);


    return 0;
}
