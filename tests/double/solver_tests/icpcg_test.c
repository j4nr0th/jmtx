// Automatically generated from tests/float/solver_tests/icpcg_test.c on Fri Dec  1 06:43:09 2023
//
// Created by jan on 22.10.2023.
//
#include <omp.h>
#include "../test_common.h"
#include "../../../include/jmtx/double/solvers/conjugate_gradient_iteration.h"
#include "../../../include/jmtx/double/matrices/sparse_row_compressed_safe.h"
#include "../../../include/jmtx/double/solvers/incomplete_cholesky_decomposition.h"

enum {PROBLEM_DIMS = (1 << 6), MAX_ITERATIONS = (PROBLEM_DIMS), CG_ITERATION_ROUND = 1};

int main()
{
    //  Make the CRS matrix for the 1D Poisson equation
    jmtxd_matrix_crs* mtx = NULL;
    jmtx_result mtx_res;
    //  Problem to solve is -d^2/dx^2 (u) = 1, with u(0) = 0 and u(1) = 0, on x in (0, 1)
    //  Exact solution is u(x) = x * (x - 1) / 2

    double* const exact_solution = calloc(PROBLEM_DIMS, sizeof(*exact_solution)); // exact solution of u
    ASSERT(exact_solution != NULL);
    double* const forcing_vector = calloc(PROBLEM_DIMS, sizeof(*forcing_vector)); // forcing vector for u (all values are 1)
    ASSERT(forcing_vector != NULL);
    double* const iterative_solution = calloc(PROBLEM_DIMS, sizeof(*iterative_solution));
    ASSERT(iterative_solution != NULL);
    double* const aux_v1 = calloc(PROBLEM_DIMS, sizeof(*aux_v1));
    ASSERT(aux_v1 != NULL);
    double* const aux_v2 = calloc(PROBLEM_DIMS, sizeof(*aux_v2));
    ASSERT(aux_v2 != NULL);
    double* const aux_v3 = calloc(PROBLEM_DIMS, sizeof(*aux_v3));
    ASSERT(aux_v3 != NULL);
    double* const aux_v4 = calloc(PROBLEM_DIMS, sizeof(*aux_v4));
    ASSERT(aux_v4 != NULL);

    omp_set_dynamic(1);
    const int proc_count = omp_get_num_procs();
    printf("OpenMP found %d processors\n", proc_count);
    
#pragma omp parallel for shared(exact_solution) default(none) schedule(static)
    for (unsigned i = 0; i < PROBLEM_DIMS; ++i)
    {
        const double x = (double)i / (double)(PROBLEM_DIMS - 1);
        exact_solution[i] = x * (x - 1) / 2;
    }
    const double dx = 1.0f / (double)(PROBLEM_DIMS - 1);
    const double dx2 = dx * dx;


    MATRIX_TEST_CALL(jmtxds_matrix_crs_new(&mtx, PROBLEM_DIMS - 2, PROBLEM_DIMS - 2, 3 * PROBLEM_DIMS, NULL));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    //  Build the matrix
    {
        const uint32_t indices1[2] = {0, 1};
        const double values1[2] = {2.0f / dx2, -1.0f / dx2};
        const uint32_t indices2[2] = {PROBLEM_DIMS - 4, PROBLEM_DIMS - 3};
        const double values2[2] = {-1.0f / dx2, 2.0f / dx2};
        forcing_vector[0] = 0.0f;
        forcing_vector[PROBLEM_DIMS - 1] = 0.0f;
        ASSERT(JMTX_RESULT_SUCCESS == (jmtxds_matrix_crs_set_row(mtx, 0, 2, indices1, values1)));
        ASSERT(JMTX_RESULT_SUCCESS == (jmtxds_matrix_crs_set_row(mtx, PROBLEM_DIMS - 3, 2, indices2, values2)));
    }
    for (unsigned i = 1; i < PROBLEM_DIMS - 3; ++i)
    {
        const uint32_t indices[3] = {i - 1, i, i + 1};
        const double values[3] = { -1.0f / dx2, 2.0f / dx2, -1.0f / dx2 };
        ASSERT(JMTX_RESULT_SUCCESS == (jmtxds_matrix_crs_set_row(mtx, i, 3, indices, values)));
    }

#pragma omp parallel for shared(forcing_vector, mtx, exact_solution) default(none) schedule(static)
    for (unsigned i = 0; i < PROBLEM_DIMS - 2; ++i)
    {
        forcing_vector[i + 1] = jmtxd_matrix_crs_vector_multiply_row(mtx, exact_solution + 1, i);
    }
    jmtxd_matrix_crs* cho = NULL;
    MATRIX_TEST_CALL(jmtxd_incomplete_cholensk_crs(mtx, &cho, NULL));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);

//    print_crs_matrix(cho);

    jmtxd_matrix_crs* cho_t;
    MATRIX_TEST_CALL(jmtxds_matrix_crs_transpose(cho, &cho_t, NULL));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);

    uint32_t total_iterations = 0;
    double total_time = 0;
    jmtxd_solver_arguments solver_arguments =
            {
            .in_convergence_criterion = 1e-6f,
            .in_max_iterations = MAX_ITERATIONS,
            };
    for (unsigned i = 0; i < CG_ITERATION_ROUND; ++i)
    {
        const double t0 = omp_get_wtime();
        mtx_res = jmtxd_incomplete_cholesky_preconditioned_conjugate_gradient_crs(
                mtx, cho, cho_t, forcing_vector + 1, iterative_solution + 1, 1e-6f, 1, aux_v1, aux_v2, aux_v3, aux_v4, &solver_arguments);
        const double t1 = omp_get_wtime();
        printf("Solution took %g seconds (%u iterations) for a problem of size %d (outcome: %s), error ratio: %g\n", t1 - t0, solver_arguments.out_last_iteration, PROBLEM_DIMS,
               jmtx_result_to_str(mtx_res), solver_arguments.out_last_error);
        ASSERT(mtx_res == JMTX_RESULT_SUCCESS || mtx_res == JMTX_RESULT_NOT_CONVERGED || mtx_res == JMTX_RESULT_STAGNATED);
        total_iterations += solver_arguments.out_last_iteration;
        total_time += t1 - t0;
    }

    iterative_solution[0] = 0;
    iterative_solution[PROBLEM_DIMS - 1] = 0;
    printf("Iterative solution had final residual ratio of %g after %u iterations\n", solver_arguments.out_last_error, total_iterations);

    for (unsigned i = 0; i < PROBLEM_DIMS; ++i)
    {
        const double x = (double)i / (double)(PROBLEM_DIMS - 1);
        printf("u_ex(%g) = %g, u_num(%g) = %g\n", x, exact_solution[i], x, iterative_solution[i]);
    }
    MATRIX_TEST_CALL(jmtxds_matrix_crs_destroy(cho_t));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    MATRIX_TEST_CALL(jmtxds_matrix_crs_destroy(cho));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    MATRIX_TEST_CALL(jmtxds_matrix_crs_destroy(mtx));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);

    free(exact_solution);
    free(forcing_vector);
    free(iterative_solution);
    free(aux_v1);
    free(aux_v2);
    free(aux_v3);
    free(aux_v4);

    return 0;
}