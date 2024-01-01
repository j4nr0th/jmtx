// Automatically generated from tests/float/solver_tests/bicgstab_cds_test.c on Sun Dec 17 16:07:06 2023
//
// Created by jan on 3.12.2023.
//
#include <inttypes.h>
#include <omp.h>
#include <math.h>
#include "../test_common.h"
#include "../../../include/jmtx/cfloat/matrices/sparse_diagonal_compressed_safe.h"
#include "../../../include/jmtx/cfloat/solvers/recursive_generalized_minimum_residual_iteration.h"

enum {PROBLEM_DIMS = (1 << 8), MAX_ITERATIONS = (PROBLEM_DIMS << 2), RESTART_INTERVAL = (1 << 4), TRUNCATION_INTERVAL = (1 << 4)};

int main()
{
    //  Make the CRS matrix for the 1D Poisson equation
    jmtxc_matrix_cds* mtx = NULL;
    jmtx_result mtx_res;
    //  Problem to solve is d^2/dx^2 (u) - d/dx (u) = -pi^2 sin(pi x) - pi cos(pi x), with u(0) = 0 and u(1) = 0,
    //  on x in (0, 1). Exact solution is u(x) = sin(pi x).
    _Complex float* const exact_solution = calloc(PROBLEM_DIMS, sizeof(*exact_solution)); // exact solution of u
    ASSERT(exact_solution);
    _Complex float* const forcing_vector = calloc(PROBLEM_DIMS, sizeof(*forcing_vector)); // forcing vector for u (all values are 1)
    ASSERT(forcing_vector);
    _Complex float* const iterative_solution = calloc(PROBLEM_DIMS, sizeof(*iterative_solution));
    ASSERT(iterative_solution);
    _Complex float* const aux_v1 = calloc(RESTART_INTERVAL, sizeof(*aux_v1));
    ASSERT(aux_v1);
    _Complex float* const aux_v2 = calloc(RESTART_INTERVAL, sizeof(*aux_v2));
    ASSERT(aux_v2);
    _Complex float* const aux_v3 = calloc(RESTART_INTERVAL, sizeof(*aux_v3));
    ASSERT(aux_v3);
    _Complex float* const aux_v4 = calloc(RESTART_INTERVAL, sizeof(*aux_v4));
    ASSERT(aux_v4);
    _Complex float* const aux_v5 = calloc(RESTART_INTERVAL, sizeof(*aux_v5));
    ASSERT(aux_v5);
    _Complex float* const aux_v6 = calloc(PROBLEM_DIMS, sizeof(*aux_v6));
    ASSERT(aux_v6);
    _Complex float* const aux_v7 = calloc(RESTART_INTERVAL, sizeof(*aux_v7));
    ASSERT(aux_v7);
    _Complex float* const aux_vecs = calloc(PROBLEM_DIMS * RESTART_INTERVAL, sizeof(*aux_vecs));
    ASSERT(aux_vecs);
    _Complex float* const aux_vecs2 = calloc(PROBLEM_DIMS * TRUNCATION_INTERVAL, sizeof(*aux_vecs2));
    ASSERT(aux_vecs2);
    _Complex float* const aux_vecs3 = calloc(PROBLEM_DIMS * TRUNCATION_INTERVAL, sizeof(*aux_vecs3));
    ASSERT(aux_vecs3);
    float* const err_evol = calloc(MAX_ITERATIONS, sizeof(*err_evol));
    ASSERT(err_evol);


    const float dx = 1.0f / (float)(PROBLEM_DIMS - 1);
    for (unsigned i = 0; i < PROBLEM_DIMS; ++i)
    {
        const float x = (float)i / (float)(PROBLEM_DIMS - 1);
        exact_solution[i] = sinf(M_PI * x);
        forcing_vector[i] = (float)(-M_PI * M_PI * sinf(M_PI * x) - M_PI * cosf(M_PI * x));
    }


    MATRIX_TEST_CALL(jmtxcs_matrix_cds_new(&mtx, PROBLEM_DIMS - 2, PROBLEM_DIMS - 2, 3, (int32_t[]){-1, 0, +1}, NULL));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    //  Build the matrix
    uint32_t len;
    _Complex float* dia = jmtxc_matrix_cds_allocate_diagonal(mtx, 0, &len);
    //  Main diagonal
    for (uint_fast32_t i = 0; i < len; ++i)
    {
        dia[i] = -2.0f / (dx * dx);
    }
    dia = jmtxc_matrix_cds_allocate_diagonal(mtx, +1, &len);
    //  Superdiagonal
    for (uint_fast32_t i = 0; i < len; ++i)
    {
        dia[i] = 1.0f / (dx * dx) + 1.0f / (2.0f * dx);
    }
    dia = jmtxc_matrix_cds_allocate_diagonal(mtx, -1, &len);
    //  Subdiagonal
    for (uint_fast32_t i = 0; i < len; ++i)
    {
        dia[i] = 1.0f / (dx * dx) - 1.0f / (2.0f * dx);
    }


    jmtxc_matrix_brm* r_mtx;
    MATRIX_TEST_CALL(jmtxc_matrix_brm_new(&r_mtx, RESTART_INTERVAL, RESTART_INTERVAL, RESTART_INTERVAL - 1, 0, NULL, NULL));

//    print_cds_matrix(mtx);
    uint32_t total_iterations = 0;
    double total_time = 0;
    jmtx_solver_arguments solver_arguments =
            {
            .in_convergence_criterion = 1e-6f,
            .in_max_iterations = MAX_ITERATIONS,
            .opt_error_evolution = err_evol,
            };

    {
        const double t0 = omp_get_wtime();
        mtx_res = jmtxc_solve_iterative_gmresr_cds(
                mtx, forcing_vector + 1, iterative_solution + 1, RESTART_INTERVAL, TRUNCATION_INTERVAL, r_mtx, aux_v1, aux_v2, aux_v3, aux_v4, aux_v5, aux_v6, aux_vecs, aux_vecs2, aux_vecs3, &solver_arguments);
        const double t1 = omp_get_wtime();
//        printf("Error evolution:\n");
//        for (uint_fast32_t j = 0; j < solver_arguments.out_last_iteration; ++j)
//        {
//            printf("%"PRIuFAST32": %.10e\n", j, err_evol[j]);
//        }
        ASSERT(mtx_res == JMTX_RESULT_SUCCESS || mtx_res == JMTX_RESULT_NOT_CONVERGED || mtx_res == JMTX_RESULT_STAGNATED);
        printf("Solution took %g seconds (%u iterations) for a problem of size %d (outcome: %s), error ratio: %g\n\n", t1 - t0, solver_arguments.out_last_iteration, PROBLEM_DIMS,
               jmtx_result_to_str(mtx_res), solver_arguments.out_last_error);

        total_iterations += solver_arguments.out_last_iteration;
        total_time += t1 - t0;
    }
    jmtxc_matrix_brm_destroy(r_mtx);

    iterative_solution[0] = 0;
    iterative_solution[PROBLEM_DIMS - 1] = 0;
    printf("Iterative solution had final residual ratio of %g after %u iterations (at most %u GMRES rounds)\n", solver_arguments.out_last_error, total_iterations, total_iterations * RESTART_INTERVAL);

    float rms_error = 0;
    for (unsigned i = 0; i < PROBLEM_DIMS; ++i)
    {
//        const float x = (float)i / (float)(PROBLEM_DIMS - 1);
//        printf("u_ex(%g) = %g, u_num(%g) = %g\n", x, exact_solution[i], x, iterative_solution[i]);
        const float local_error = exact_solution[i] - iterative_solution[i];
        rms_error += local_error * local_error;
    }
    rms_error = sqrtf(rms_error / PROBLEM_DIMS);
    printf("Iterative solution had final RMS error of %e after %u iterations\n", rms_error, total_iterations);

    MATRIX_TEST_CALL(jmtxcs_matrix_cds_destroy(mtx));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    free(err_evol);
    free(exact_solution);
    free(forcing_vector);
    free(iterative_solution);
    free(aux_v1);
    free(aux_v2);
    free(aux_v3);
    free(aux_v4);
    free(aux_v5);
    free(aux_v6);
    free(aux_v7);
    free(aux_vecs);
    free(aux_vecs2);
    free(aux_vecs3);
    return 0;
}