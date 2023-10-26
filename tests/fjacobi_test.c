//
// Created by jan on 12.7.2023.
//

#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
//#include <cplot.h>
#include <math.h>
#include <time.h>
#include "test_common.h"
#include "../source/solvers/jacobi_point_iteration.h"
#include "../source/matrices/sparse_row_compressed_internal.h"
#include "../source/f_exper/f_interface.h"

static const float V = 10.0f;
static const float omega_0 = 5.0f;

#define steps 2048
#define iterations 2048

static double ts_difference(const struct timespec* t_begin, const struct timespec* t_end)
{
    return (double)(t_end->tv_sec - t_begin->tv_sec) + (double)(t_end->tv_nsec - t_begin->tv_nsec) * 1e-9;
}

int main()
{
    static float x_values[steps] = {0};
    static float y_exact[steps] = {0};
    static float y_approx[steps] = {0};
    static float y_relax[steps] = {0};
    static float y_approx10[steps] = {0};
    static float y_approx20[steps] = {0};
    static float y_approx30[steps] = {0};
    static float y_approx40[steps] = {0};
    static float y_approx50[steps] = {0};
    static float y_approx60[steps] = {0};
    static float y_approx_f[steps] = {0};
    static float f_exact[steps] = {0};
    static float error_evol[iterations] = {0};
    static float error_relaxed[iterations] = {0};
    static float final_error;


    const float dx = 1.0f / (float)(steps - 1);

    for (unsigned i = 0; i < steps; ++i)
    {
        float x = (float)i / (float)(steps - 1);
        x_values[i] = x;
        f_exact[i] = V * cosf(omega_0 * x);
        y_exact[i] = V / (1 + omega_0 * omega_0) * (omega_0 * sinf(omega_0 * x) - cosf(omega_0 * x) + expf(x));
    }

    jmtx_matrix_crs* matrix;
    jmtx_result mtx_res;

    MATRIX_TEST_CALL(jmtx_matrix_crs_new(&matrix, steps, steps, 2 * steps, NULL));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    MATRIX_TEST_CALL(jmtx_matrix_crs_set_entry(matrix, 0, 0, 1.0f));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    for (uint32_t i = 1; i < steps; ++i)
    {
        uint32_t indices[2] = {i - 1, i};
        float values[2] = {- 1 / dx, 1 / dx - 1 };
        MATRIX_TEST_CALL(jmtx_matrix_crs_set_row(matrix, i, 2, indices, values));
        ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    }
//    print_crs_matrix(matrix);

    uint32_t iter_count, iter_relax;
    f_exact[0] = 0.0f;
    MATRIX_TEST_CALL(jmtx_jacobi_crs(matrix, f_exact, y_approx10, 1e-4f, (iterations / 5), NULL, NULL, NULL, NULL));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS || mtx_res == JMTX_RESULT_NOT_CONVERGED);
    MATRIX_TEST_CALL(jmtx_jacobi_crs(matrix, f_exact, y_approx20, 1e-4f, (iterations * 2 / 5), NULL, NULL, NULL, NULL));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS || mtx_res == JMTX_RESULT_NOT_CONVERGED);
    MATRIX_TEST_CALL(jmtx_jacobi_crs(matrix, f_exact, y_approx30, 1e-4f, (iterations * 3 / 5), NULL, NULL, NULL, NULL));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS || mtx_res == JMTX_RESULT_NOT_CONVERGED);
    MATRIX_TEST_CALL(jmtx_jacobi_crs(matrix, f_exact, y_approx40, 1e-4f, (iterations * 4 / 5), NULL, NULL, NULL, NULL));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS || mtx_res == JMTX_RESULT_NOT_CONVERGED);

    printf("Using normal Jacobi iteration\n");
    MATRIX_TEST_CALL(
            jmtx_jacobi_crs(matrix, f_exact, y_approx, 1e-4f, iterations, &iter_count, error_evol, &final_error, NULL));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS || mtx_res == JMTX_RESULT_NOT_CONVERGED);
    printf("Afer %"PRIu32" iterations, the final error was %g\n", iter_count, final_error);

    f_pt_jacobi((int32_t)matrix->n_entries, (int32_t)matrix->base.rows, (const int32_t*)matrix->end_of_row_offsets, (const int32_t*)matrix->indices, matrix->values, 1e-4f, iterations, y_approx_f, f_exact);

    struct timespec ts_0, ts_1;
    memset(y_approx, 0, sizeof(y_approx));
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &ts_0);
    jmtx_jacobi_crs(matrix, f_exact, y_approx, 1e-4f, iterations, &iter_count, error_evol, &final_error, NULL);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &ts_1);
    printf("Time needed for jmtx_jacobi_crs to solve the problem: %g ms\n", 1e3 * ts_difference(&ts_0, &ts_1));

    const float relax_factor = 1.0f;
    printf("Using a relaxation factor of %g\n", relax_factor);
    MATRIX_TEST_CALL(jmtx_jacobi_relaxed_crs(
            matrix, f_exact, y_relax, relax_factor, 1e-4f, iterations, &iter_relax, error_relaxed, &final_error, NULL));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS || mtx_res == JMTX_RESULT_NOT_CONVERGED);
    printf("Afer %"PRIu32" iterations, the final error was %g\n", iter_count, final_error);

    memset(y_relax, 0, sizeof(y_relax));
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &ts_0);
    jmtx_jacobi_relaxed_crs(
            matrix, f_exact, y_relax, relax_factor, 1e-4f, iterations, &iter_count, error_evol, &final_error, NULL);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &ts_1);
    printf("Time needed for jmtx_jacobi_relaxed_crs to solve the problem: %g ms\n", 1e3 * ts_difference(&ts_0, &ts_1));



    MATRIX_TEST_CALL(jmtx_matrix_crs_destroy(matrix));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
#ifdef PLOT_RESULTS
    cplot_context* cplot_ctx = cplot_context_create();
    ASSERT(cplot_ctx);
    cplot_figure* fig = cplot_figure_new(cplot_ctx, "Jacobi solution", 1600, 900);
    ASSERT(fig);
    
    cplot_subfigure* sfg_exact = cplot_subfigure_new(fig, 3, 2, 1, 3);
    ASSERT(sfg_exact);
    cplot_curve2d* exact_plot = cplot_plot(
            cplot_ctx, sfg_exact, steps, x_values, y_exact, (cplot_color) { .b = 0xFF, .a = 0xFF }, "Exact");
    ASSERT(exact_plot);
    (void)cplot_limit_x(sfg_exact, 0.0f, 1.0f);
    (void)cplot_limit_y(sfg_exact, -5.0f, 5.0f);
    const cplot_color black_color = (cplot_color){.a = 0xFF};
    (void)cplot_grid(sfg_exact, 1, &black_color);
    (void)cplot_title(sfg_exact, "Exact solution");
    
    cplot_subfigure* sfg_apporx = cplot_subfigure_new(fig, 3, 2, 2, 4);
    ASSERT(sfg_apporx);

    cplot_curve2d* approx10_plot = cplot_plot(
            cplot_ctx, sfg_apporx, steps, x_values, y_approx10, (cplot_color) { .b = 255, .r = 0, .a = 0xFF }, "20 %");
    ASSERT(approx10_plot);

    cplot_curve2d* approx20_plot = cplot_plot(
            cplot_ctx, sfg_apporx, steps, x_values, y_approx20, (cplot_color) { .b = 215, .r = 40, .a = 0xFF }, "40 %");
    ASSERT(approx20_plot);

    cplot_curve2d* approx30_plot = cplot_plot(
            cplot_ctx, sfg_apporx, steps, x_values, y_approx30, (cplot_color) { .b = 175, .r = 80, .a = 0xFF }, "60 %");
    ASSERT(approx30_plot);

    cplot_curve2d* approx40_plot = cplot_plot(
            cplot_ctx, sfg_apporx, steps, x_values, y_approx40, (cplot_color) { .b = 135, .r = 120,.a = 0xFF }, "80 %");
    ASSERT(approx40_plot);

    char name_buffer[64];
    snprintf(name_buffer, sizeof(name_buffer), "Converged: %d", iter_count);
    cplot_curve2d* approx_plot = cplot_plot(
            cplot_ctx, sfg_apporx, steps, x_values, y_approx, (cplot_color) { .r = 0xFF, .a = 0xFF }, name_buffer);
    ASSERT(approx_plot);

    (void)cplot_legend(sfg_apporx, cplot_legend_bottom|cplot_legend_left);



    (void)cplot_limit_x(sfg_apporx, 0.0f, 1.0f);
    (void)cplot_limit_y(sfg_apporx, -5.0f, 5.0f);
    (void)cplot_grid(sfg_apporx, 1, &black_color);
    (void)cplot_title(sfg_apporx, "Approximated solution");

    cplot_subfigure* sfg_error = cplot_subfigure_new(fig, 3, 2, 5, 6);
    ASSERT(sfg_error);
    float iter_count_array[iterations];
    for (uint32_t i = 0; i < iterations; ++i)
    {
        iter_count_array[i] = (float)(i + 1);
    }
    cplot_curve2d* err_plot = cplot_plot(
            cplot_ctx, sfg_error, iter_count, iter_count_array, error_evol, (cplot_color) { .g = 0xFF, .a = 0xFF },
            "Regular");
    ASSERT(err_plot);
    cplot_curve2d* rel_plot = cplot_plot(
            cplot_ctx, sfg_error, iter_relax, iter_count_array, error_relaxed,
            (cplot_color) { .g = 0xFF, .r = 0x80, .a = 0xFF }, "Relaxed");
    ASSERT(rel_plot);

    (void)cplot_limit_x(sfg_error, 1.0f, (float)iter_count);
//    (void)cplot_limit_y(sfg_error, 0.0f, 1.0f);
//    (void)cplot_grid(sfg_error, 1, &black_color);
    (void)cplot_title(sfg_error, "Evolution of error");
    (void)cplot_legend(sfg_error, cplot_legend_bottom|cplot_legend_left);

    int res = cplot_figure_show(cplot_ctx, fig);
    ASSERT(res == 0);
    cplot_context_destroy(cplot_ctx);
#endif
    return 0;
}
