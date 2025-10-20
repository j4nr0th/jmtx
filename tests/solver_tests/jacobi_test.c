#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>

#include <math.h>
#include <time.h>
#include "../test_common.h"
#include "matrices/sparse_row_compressed.h"
#include "solvers/jacobi_point_iteration.h"

static const double V = 10.0f;
static const double omega_0 = 5.0f;

enum
{
    STEP_COUNT = 64,
    ITERATION_COUNT = 64,
};

static double ts_difference(const struct timespec *t_begin, const struct timespec *t_end)
{
    return (double)(t_end->tv_sec - t_begin->tv_sec) + (double)(t_end->tv_nsec - t_begin->tv_nsec) * 1e-9;
}

int main()
{
    static JMTX_SCALAR_T x_values[STEP_COUNT] = {0};
    static JMTX_SCALAR_T x_aux1[STEP_COUNT] = {0};
    static JMTX_SCALAR_T x_aux2[STEP_COUNT] = {0};
    static JMTX_SCALAR_T y_exact[STEP_COUNT] = {0};
    static JMTX_SCALAR_T y_approx[STEP_COUNT] = {0};
    static JMTX_SCALAR_T y_relax[STEP_COUNT] = {0};
    static JMTX_SCALAR_T y_approx10[STEP_COUNT] = {0};
    static JMTX_SCALAR_T y_approx20[STEP_COUNT] = {0};
    static JMTX_SCALAR_T y_approx30[STEP_COUNT] = {0};
    static JMTX_SCALAR_T y_approx40[STEP_COUNT] = {0};
    static JMTX_SCALAR_T f_exact[STEP_COUNT] = {0};

    const double dx = 1.0f / (double)(STEP_COUNT - 1);

    for (unsigned i = 0; i < STEP_COUNT; ++i)
    {
        double x = (double)i / (double)(STEP_COUNT - 1);
        x_values[i] = x;
        f_exact[i] = V * cos(omega_0 * x);
        y_exact[i] = V / (1 + omega_0 * omega_0) * (omega_0 * sin(omega_0 * x) - cos(omega_0 * x) + expf(x));
    }

    JMTX_NAME_TYPED(matrix_crs) * matrix;
    jmtx_result mtx_res;

    MATRIX_TEST_CALL(JMTX_NAME_TYPED(matrix_crs_new)(&matrix, STEP_COUNT, STEP_COUNT, 2 * STEP_COUNT, NULL));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    MATRIX_TEST_CALL(JMTX_NAME_TYPED(matrix_crs_set_entry)(matrix, 0, 0, 1.0f));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    for (uint32_t i = 1; i < STEP_COUNT; ++i)
    {
        const uint32_t indices[2] = {i - 1, i};
        const JMTX_SCALAR_T values[2] = {-1 / dx, 1 / dx - 1};
        MATRIX_TEST_CALL(JMTX_NAME_TYPED(matrix_crs_set_row)(matrix, i, 2, indices, values));
        ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    }
    //    print_crsd_matrix(matrix);

    JMTX_NAME_TYPED(solver_arguments)
    solver_arguments = {
        .in_max_iterations = 1 << 10,
        .in_convergence_criterion = 1e-4f,
    };
    f_exact[0] = 0.0f;
    MATRIX_TEST_CALL(
        JMTX_NAME_TYPED(solve_iterative_jacobi_crs)(matrix, f_exact, y_approx10, x_aux1, x_aux2, &solver_arguments));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS || mtx_res == JMTX_RESULT_NOT_CONVERGED);
    MATRIX_TEST_CALL(
        JMTX_NAME_TYPED(solve_iterative_jacobi_crs)(matrix, f_exact, y_approx20, x_aux1, x_aux2, &solver_arguments));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS || mtx_res == JMTX_RESULT_NOT_CONVERGED);
    MATRIX_TEST_CALL(
        JMTX_NAME_TYPED(solve_iterative_jacobi_crs)(matrix, f_exact, y_approx30, x_aux1, x_aux2, &solver_arguments));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS || mtx_res == JMTX_RESULT_NOT_CONVERGED);
    MATRIX_TEST_CALL(
        JMTX_NAME_TYPED(solve_iterative_jacobi_crs)(matrix, f_exact, y_approx40, x_aux1, x_aux2, &solver_arguments));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS || mtx_res == JMTX_RESULT_NOT_CONVERGED);

    printf("Using normal Jacobi iteration\n");
    MATRIX_TEST_CALL(
        JMTX_NAME_TYPED(solve_iterative_jacobi_crs)(matrix, f_exact, y_approx, x_aux1, x_aux2, &solver_arguments));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS || mtx_res == JMTX_RESULT_NOT_CONVERGED);
    printf("After %" PRIu32 " ITERATION_COUNT, the final error was %g\n", solver_arguments.out_last_iteration,
           solver_arguments.out_last_error);

    struct timespec ts_0, ts_1;
    memset(y_approx, 0, sizeof(y_approx));
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &ts_0);
    JMTX_NAME_TYPED(solve_iterative_jacobi_crs)(matrix, f_exact, y_approx, x_aux1, x_aux2, &solver_arguments);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &ts_1);
    printf("Time needed for JMTX_NAME_TYPED(solve_iterative_jacobi_crs) to solve the problem: %g ms\n",
           1e3 * ts_difference(&ts_0, &ts_1));

    const double relax_factor = 1.0f;
    printf("Using a relaxation factor of %g\n", relax_factor);
    MATRIX_TEST_CALL(JMTX_NAME_TYPED(solve_iterative_jacobi_relaxed_crs)(matrix, f_exact, y_relax, relax_factor, x_aux1,
                                                                         x_aux2, &solver_arguments));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS || mtx_res == JMTX_RESULT_NOT_CONVERGED);
    printf("After %" PRIu32 " ITERATION_COUNT, the final error was %g\n", solver_arguments.out_last_iteration,
           solver_arguments.out_last_error);

    memset(y_relax, 0, sizeof(y_relax));
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &ts_0);
    JMTX_NAME_TYPED(solve_iterative_jacobi_relaxed_crs)(matrix, f_exact, y_relax, relax_factor, x_aux1, x_aux2,
                                                        &solver_arguments);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &ts_1);
    printf("Time needed for JMTX_NAME_TYPED(solve_iterative_jacobi_relaxed_crs) to solve the problem: %g ms\n",
           1e3 * ts_difference(&ts_0, &ts_1));

    JMTX_NAME_TYPED(matrix_crs_destroy)(matrix);

    return 0;
}
