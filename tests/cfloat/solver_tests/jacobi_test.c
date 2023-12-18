// Automatically generated from tests/float/solver_tests/jacobi_test.c on Fri Dec  1 17:35:45 2023
//
// Created by jan on 12.7.2023.
//

#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
//#include <cplot.h>
#include <math.h>
#include <time.h>
#include "../test_common.h"
#include "../../../include/jmtx/cfloat/matrices/sparse_row_compressed_safe.h"
#include "../../../include/jmtx/cfloat/solvers/jacobi_point_iteration.h"

static const _Complex float V = 10.0f;
static const _Complex float omega_0 = 5.0f;

#define steps 64
#define iterations 64

static double ts_difference(const struct timespec* t_begin, const struct timespec* t_end)
{
    return (double)(t_end->tv_sec - t_begin->tv_sec) + (double)(t_end->tv_nsec - t_begin->tv_nsec) * 1e-9;
}

int main()
{
    static _Complex float x_values[steps] = {0};
    static _Complex float x_aux1[steps] = {0};
    static _Complex float x_aux2[steps] = {0};
    static _Complex float y_exact[steps] = {0};
    static _Complex float y_approx[steps] = {0};
    static _Complex float y_relax[steps] = {0};
    static _Complex float y_approx10[steps] = {0};
    static _Complex float y_approx20[steps] = {0};
    static _Complex float y_approx30[steps] = {0};
    static _Complex float y_approx40[steps] = {0};
    static _Complex float f_exact[steps] = {0};


    const _Complex float dx = 1.0f / (_Complex float)(steps - 1);

    for (unsigned i = 0; i < steps; ++i)
    {
        _Complex float x = (_Complex float)i / (_Complex float)(steps - 1);
        x_values[i] = x;
        f_exact[i] = V * cosf(omega_0 * x);
        y_exact[i] = V / (1 + omega_0 * omega_0) * (omega_0 * sinf(omega_0 * x) - cosf(omega_0 * x) + expf(x));
    }

    jmtxc_matrix_crs* matrix;
    jmtx_result mtx_res;

    MATRIX_TEST_CALL(jmtxcs_matrix_crs_new(&matrix, steps, steps, 2 * steps, NULL));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    MATRIX_TEST_CALL(jmtxcs_matrix_crs_set_entry(matrix, 0, 0, 1.0f));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    for (uint32_t i = 1; i < steps; ++i)
    {
        uint32_t indices[2] = {i - 1, i};
        _Complex float values[2] = {- 1 / dx, 1 / dx - 1 };
        MATRIX_TEST_CALL(jmtxcs_matrix_crs_set_row(matrix, i, 2, indices, values));
        ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    }
//    print_crs_matrix(matrix);

    jmtx_solver_arguments solver_arguments =
            {
                    .in_max_iterations = 1 << 10,
                    .in_convergence_criterion = 1e-4f,
            };
    f_exact[0] = 0.0f;
    MATRIX_TEST_CALL(jmtxc_solve_iterative_jacobi_crs(matrix, f_exact, y_approx10, x_aux1, x_aux2, &solver_arguments));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS || mtx_res == JMTX_RESULT_NOT_CONVERGED);
    MATRIX_TEST_CALL(jmtxc_solve_iterative_jacobi_crs(matrix, f_exact, y_approx20, x_aux1, x_aux2, &solver_arguments));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS || mtx_res == JMTX_RESULT_NOT_CONVERGED);
    MATRIX_TEST_CALL(jmtxc_solve_iterative_jacobi_crs(matrix, f_exact, y_approx30, x_aux1, x_aux2, &solver_arguments));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS || mtx_res == JMTX_RESULT_NOT_CONVERGED);
    MATRIX_TEST_CALL(jmtxc_solve_iterative_jacobi_crs(matrix, f_exact, y_approx40, x_aux1, x_aux2, &solver_arguments));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS || mtx_res == JMTX_RESULT_NOT_CONVERGED);

    printf("Using normal Jacobi iteration\n");
    MATRIX_TEST_CALL(
            jmtxc_solve_iterative_jacobi_crs(matrix, f_exact, y_approx, x_aux1, x_aux2, &solver_arguments));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS || mtx_res == JMTX_RESULT_NOT_CONVERGED);
    printf("Afer %"PRIu32" iterations, the final error was %g\n", solver_arguments.out_last_iteration, solver_arguments.out_last_error);

    struct timespec ts_0, ts_1;
    memset(y_approx, 0, sizeof(y_approx));
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &ts_0);
    jmtxc_solve_iterative_jacobi_crs(matrix, f_exact, y_approx, x_aux1, x_aux2, &solver_arguments);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &ts_1);
    printf("Time needed for jmtxc_solve_iterative_jacobi_crs to solve the problem: %g ms\n", 1e3 * ts_difference(&ts_0, &ts_1));

    const float relax_factor = 1.0f;
    printf("Using a relaxation factor of %g\n", relax_factor);
    MATRIX_TEST_CALL(jmtxc_solve_iterative_jacobi_relaxed_crs(
            matrix, f_exact, y_relax, relax_factor, x_aux1, x_aux2, &solver_arguments));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS || mtx_res == JMTX_RESULT_NOT_CONVERGED);
    printf("Afer %"PRIu32" iterations, the final error was %g\n", solver_arguments.out_last_iteration, solver_arguments.out_last_error);

    memset(y_relax, 0, sizeof(y_relax));
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &ts_0);
    jmtxc_solve_iterative_jacobi_relaxed_crs(
            matrix, f_exact, y_relax, relax_factor, x_aux1, x_aux2, &solver_arguments);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &ts_1);
    printf("Time needed for jmtxc_solve_iterative_jacobi_relaxed_crs to solve the problem: %g ms\n", 1e3 * ts_difference(&ts_0, &ts_1));



    MATRIX_TEST_CALL(jmtxcs_matrix_crs_destroy(matrix));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);

    return 0;
}
