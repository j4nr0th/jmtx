// Automatically generated from tests/float/solver_tests/omp_inclu_test.c on Fri Dec  1 06:43:09 2023
//
// Created by jan on 2.11.2023.
//
#include <omp.h>
#include <stdio.h>
#include "../test_common.h"

#include <math.h>
#include "inttypes.h"
#include "../../../include/jmtx/double/matrices/sparse_row_compressed_safe.h"
#include "../../../source/matrices/sparse_column_compressed_safe.h"
#include "../../../source/decompositions/incomplete_lu_decomposition.h"
#include "../../../source/matrices/sparse_multiplication.h"
#include "../../../source/solvers/lu_solving.h"
#include "../../../include/jmtx/double/matrices/sparse_conversion.h"

enum
{
    PROBLEM_SIZE_X = 1 << 2,
    PROBLEM_SIZE_Y = 1 << 2,
    INTERNAL_SIZE_X = PROBLEM_SIZE_X - 2,
    INTERNAL_SIZE_Y = PROBLEM_SIZE_Y - 2,
    PROBLEM_INTERNAL_PTS = INTERNAL_SIZE_X * INTERNAL_SIZE_Y,
    WORK_DIVISIONS = 4,
    MAXIMUM_ITERATIONS = (1 << 10),
};

static unsigned lexicographic_position(unsigned i, unsigned j)
{
    return INTERNAL_SIZE_X * i + j;
}

static void from_lexicographic(unsigned n, unsigned *pi, unsigned *pj)
{
    *pj = n % INTERNAL_SIZE_X;
    *pi = n / INTERNAL_SIZE_X;
}

int main()
{
    JMTX_NAME_TYPED(matrix_crs) * mtx;
    jmtx_result mtx_res;
    omp_set_dynamic(0);
    omp_set_num_threads(1);
    const int proc_count = omp_get_num_procs();
    const int max_threads = omp_get_max_threads();
    printf("OpenMP found %d processors, with a maximum of %d threads\n", proc_count, max_threads);

    const double dy = 1.0f / (PROBLEM_SIZE_Y - 1);
    const double dx = 1.0f / (PROBLEM_SIZE_X - 1);

    const double rdy2 = 1.0f / (dy * dy);
    const double rdx2 = 1.0f / (dx * dx);

    MATRIX_TEST_CALL(jmtxds_matrix_crs_new(&mtx, PROBLEM_INTERNAL_PTS, PROBLEM_INTERNAL_PTS,
                                           5 * PROBLEM_INTERNAL_PTS < PROBLEM_INTERNAL_PTS * PROBLEM_INTERNAL_PTS
                                               ? 5 * PROBLEM_INTERNAL_PTS
                                               : PROBLEM_INTERNAL_PTS * PROBLEM_INTERNAL_PTS,
                                           NULL));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    //  Serial construction
    const double t0_serial = omp_get_wtime();
    for (unsigned i = 0; i < INTERNAL_SIZE_Y; ++i)
    {
        for (unsigned j = 0; j < INTERNAL_SIZE_X; ++j)
        {
            //  Point is (DX * j, DY * i)
            unsigned k = 0;
            double values[5];
            uint32_t positions[5];
            if (i != 0)
            {
                // There's a bottom boundary
                values[k] = -rdy2;
                positions[k] = lexicographic_position(i - 1, j);
                k += 1;
            }

            if (j != 0)
            {
                // There's a left boundary
                values[k] = -rdx2;
                positions[k] = lexicographic_position(i, j - 1);
                k += 1;
            }

            values[k] = 2 * (rdx2 + rdy2);
            positions[k] = lexicographic_position(i, j);
            k += 1;

            if (j != INTERNAL_SIZE_X - 1)
            {
                // There's a right boundary
                values[k] = -rdx2;
                positions[k] = lexicographic_position(i, j + 1);
                k += 1;
            }

            if (i != INTERNAL_SIZE_Y - 1)
            {
                // There's a top boundary
                values[k] = -rdy2;
                positions[k] = lexicographic_position(i + 1, j);
                k += 1;
            }

            jmtxd_matrix_crs_build_row(mtx, lexicographic_position(i, j), k, positions, values);
        }
    }
    const double t1_serial = omp_get_wtime();

    printf("Serial construction of a %d by %d matrix took %g seconds\n", PROBLEM_INTERNAL_PTS, PROBLEM_INTERNAL_PTS,
           t1_serial - t0_serial);

    JMTX_NAME_TYPED(matrix_crs) *lower = NULL;
    JMTX_NAME_TYPED(matrix_ccs) *upper = NULL;
    const double t0_decomp = omp_get_wtime();
    MATRIX_TEST_CALL(jmtxd_decompose_ilu_crs(mtx, &lower, &upper, NULL));
    const double t1_decomp = omp_get_wtime();
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS || mtx_res == JMTX_RESULT_NOT_CONVERGED);

    printf("Decomposition took %g seconds and the result: %s\n", t1_decomp - t0_decomp, jmtx_result_to_str(mtx_res));

    double *const initial_vector = malloc(PROBLEM_INTERNAL_PTS * sizeof(*initial_vector));
    ASSERT(initial_vector != NULL);

    double *const forcing_vector = malloc(PROBLEM_INTERNAL_PTS * sizeof(*forcing_vector));
    ASSERT(forcing_vector != NULL);

    double *const approximate_vector = malloc(PROBLEM_INTERNAL_PTS * sizeof(*approximate_vector));
    ASSERT(approximate_vector != NULL);

    double *const auxiliary_vector = malloc(PROBLEM_INTERNAL_PTS * sizeof(*auxiliary_vector));
    ASSERT(auxiliary_vector != NULL);

    double mag_y = 0;
    for (unsigned i = 0; i < INTERNAL_SIZE_Y; ++i)
    {
        for (unsigned j = 0; j < INTERNAL_SIZE_X; ++j)
        {
            const double v = ((double)(i + 1) * dx) + ((double)(j + 1) * dy);
            initial_vector[i * INTERNAL_SIZE_X + j] = v;
            mag_y += v * v;
        }
    }
    mag_y = sqrt(mag_y);

    MATRIX_TEST_CALL(jmtxds_matrix_crs_vector_multiply(mtx, initial_vector, forcing_vector));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);

    JMTX_NAME_TYPED(matrix_crs) * upper_crs;
    MATRIX_TEST_CALL(jmtxd_convert_ccs_to_crs(upper, &upper_crs, NULL));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);

    JMTX_NAME_TYPED(solver_arguments)
    solve_args = {
        .in_max_iterations = MAXIMUM_ITERATIONS,
        .in_convergence_criterion = 1e-4f,
    };
    MATRIX_TEST_CALL(jmtxd_solve_iterative_ilu_crs_precomputed_parallel(
        mtx, lower, upper_crs, forcing_vector, approximate_vector, auxiliary_vector, &solve_args));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS || mtx_res == JMTX_RESULT_NOT_CONVERGED);
    printf("Solving using ILU took %" PRIu32 " iterations, with the final error of %g\n", solve_args.out_last_iteration,
           (double)solve_args.out_last_error);

    MATRIX_TEST_CALL(jmtxds_matrix_crs_destroy(upper_crs));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);

    printf("Comparison of approximate solution vs real solution:\n");
    double rms_err = 0;
    double residual = 0;
    for (unsigned i = 0; i < PROBLEM_INTERNAL_PTS; ++i)
    {
        forcing_vector[i] -= jmtxd_matrix_crs_vector_multiply_row(mtx, approximate_vector, i);
        const double err = (initial_vector[i] - approximate_vector[i]) / initial_vector[i];
        rms_err += err * err;
        //        printf("Element %u, real: %g, approx: %g, err: %g, residual: %g\n", i, initial_vector[i],
        //        approximate_vector[i],
        //               err, forcing_vector[i]);
        residual += forcing_vector[i] * forcing_vector[i];
    }
    rms_err = sqrt(rms_err / PROBLEM_INTERNAL_PTS);
    residual = sqrt(residual);
    printf("RMS Error: %g, residual/forcing: %g/%g = %g\n", rms_err, residual, mag_y, residual / mag_y);

    free(auxiliary_vector);
    free(approximate_vector);
    free(forcing_vector);
    free(initial_vector);

    MATRIX_TEST_CALL(jmtxds_matrix_crs_destroy(lower));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    MATRIX_TEST_CALL(jmtxds_matrix_ccs_destroy(upper));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);

    MATRIX_TEST_CALL(jmtxds_matrix_crs_destroy(mtx));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    return 0;
}
