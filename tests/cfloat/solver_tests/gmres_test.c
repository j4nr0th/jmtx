// Automatically generated from tests/float/solver_tests/omp_cg_test.c on Fri Dec  1 06:43:09 2023
//
// Created by jan on 22.10.2023.
//
#include <inttypes.h>
#include "../test_common.h"
#include "../../../include/jmtx/cfloat/matrices/sparse_row_compressed_safe.h"
#include "../../../include/jmtx/cfloat/solvers/generalized_minimum_residual_iteration.h"
#include <omp.h>

enum {PROBLEM_DIMS = (1 << 8), MAX_ITERATIONS = PROBLEM_DIMS, RESTART_INTERVAL = PROBLEM_DIMS, CG_ITERATION_ROUND = 1};

int main()
{
    //  Make the CRS matrix for the 1D Poisson equation
    jmtxc_matrix_crs* mtx = NULL;
    jmtx_result mtx_res;
    //  Problem to solve is d^2/dx^2 (u) = 1, with u(0) = 0 and u(1) = 0, on x in (0, 1)
    //  Exact solution is u(x) = x * (x - 1) / 2
    _Complex float* const exact_solution = calloc(sizeof(*exact_solution), PROBLEM_DIMS); // exact solution of u
    ASSERT(exact_solution);
    _Complex float* const forcing_vector = calloc(sizeof(*forcing_vector), PROBLEM_DIMS); // forcing vector for u (all values are 1)
    ASSERT(forcing_vector);
    _Complex float* const iterative_solution = calloc(sizeof(*iterative_solution), PROBLEM_DIMS);
    ASSERT(iterative_solution);
    _Complex float* const aux_v1 = calloc(sizeof(*aux_v1), RESTART_INTERVAL);
    ASSERT(aux_v1);
    _Complex float* const aux_v2 = calloc(sizeof(*aux_v2), RESTART_INTERVAL);
    ASSERT(aux_v2);
    _Complex float* const aux_v3 = calloc(sizeof(*aux_v3), RESTART_INTERVAL);
    ASSERT(aux_v3);
    _Complex float* const aux_v4 = calloc(sizeof(*aux_v4), RESTART_INTERVAL);
    ASSERT(aux_v4);
    _Complex float* const aux_v5 = calloc(sizeof(*aux_v5), RESTART_INTERVAL);
    ASSERT(aux_v5);
    _Complex float* const aux_vecs = calloc(sizeof(*aux_vecs), RESTART_INTERVAL * PROBLEM_DIMS);
    ASSERT(aux_vecs);

    for (unsigned i = 0; i < PROBLEM_DIMS; ++i)
    {
        const float x = (float)i / (float)(PROBLEM_DIMS - 1);
        exact_solution[i] = x * (x - 1) / 2;
    }
    const float dx = 1.0f / (float)(PROBLEM_DIMS - 1);

    for (unsigned i = 0; i < PROBLEM_DIMS; ++i)
    {
        forcing_vector[i] = 1.0f * dx * dx;
    }

    MATRIX_TEST_CALL(jmtxcs_matrix_crs_new(&mtx, PROBLEM_DIMS - 2, PROBLEM_DIMS - 2, 3 * PROBLEM_DIMS, NULL));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    //  Build the matrix
    {
        const uint32_t indices1[2] = {0, 1};
        const _Complex float values1[2] = {-2.0f, 1.0f};
        const uint32_t indices2[2] = {PROBLEM_DIMS - 4, PROBLEM_DIMS - 3};
        const _Complex float values2[2] = {1.0f, -2.0f};
        forcing_vector[0] += -0.0f;
        forcing_vector[PROBLEM_DIMS - 1] += -0.0f;
        ASSERT(mtx_res == (jmtxcs_matrix_crs_set_row(mtx, 0, 2, indices1, values1)));
        ASSERT(mtx_res == (jmtxcs_matrix_crs_set_row(mtx, PROBLEM_DIMS - 3, 2, indices2, values2)));
    }
    for (unsigned i = 1; i < PROBLEM_DIMS - 3; ++i)
    {
        const uint32_t indices[3] = {i - 1, i, i + 1};
        const _Complex float values[3] = { 1.0f, -2.0f, 1.0f };
        ASSERT(mtx_res == (jmtxcs_matrix_crs_set_row(mtx, i, 3, indices, values)));
    }
//    print_crsd_matrix(mtx);
    jmtxc_matrix_brm* r_mtx = NULL;
    MATRIX_TEST_CALL(jmtxc_matrix_brm_new(&r_mtx, RESTART_INTERVAL, RESTART_INTERVAL, RESTART_INTERVAL - 1, 0, NULL, NULL));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    uint32_t total_iterations = 0;
    double total_time = 0;
    jmtx_solver_arguments solver_arguments =
            {
            .in_convergence_criterion = 1e-4f,
            .in_max_iterations = MAX_ITERATIONS,
            };

    {
        const double t0 = omp_get_wtime();
        mtx_res = jmtxc_solve_iterative_gmresm_crs(
                mtx, forcing_vector, iterative_solution + 1, RESTART_INTERVAL, r_mtx, aux_v1, aux_v2, aux_v3, aux_v4, aux_v5, aux_vecs, &solver_arguments);
        const double t1 = omp_get_wtime();
        printf("Solution took %g seconds (%u iterations) for a problem of size %d (outcome: %s), error ratio: %g\n", t1 - t0, solver_arguments.out_last_iteration, PROBLEM_DIMS,
               jmtx_result_to_str(mtx_res), solver_arguments.out_last_error);
        ASSERT(mtx_res == JMTX_RESULT_SUCCESS || mtx_res == JMTX_RESULT_NOT_CONVERGED || mtx_res == JMTX_RESULT_STAGNATED);
        total_iterations += solver_arguments.out_last_iteration;
        total_time += t1 - t0;
    }
    jmtxc_matrix_brm_destroy(r_mtx);

    iterative_solution[0] = 0;
    iterative_solution[PROBLEM_DIMS - 1] = 0;
    printf("Iterative solution had final residual ratio of %g after %u iterations\n", solver_arguments.out_last_error, total_iterations);

    MATRIX_TEST_CALL(jmtxcs_matrix_crs_destroy(mtx));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);

    free(aux_vecs);
    free(aux_v5);
    free(aux_v4);
    free(aux_v3);
    free(aux_v2);
    free(aux_v1);
    free(iterative_solution);
    free(forcing_vector);
    free(exact_solution);
    return 0;
}