// Automatically generated from tests/float/solver_tests/bicgstab_cds_test.c on Sun Dec 17 16:07:06 2023
//
// Created by jan on 3.12.2023.
//
#include <inttypes.h>
#include <omp.h>
#include <math.h>
#include "../test_common.h"
#include "../../../include/jmtx/double/matrices/sparse_diagonal_compressed_safe.h"
#include "../../../include/jmtx/double/matrices/sparse_row_compressed_safe.h"
#include "../../../include/jmtx/double/matrices/sparse_column_compressed_safe.h"
#include "../../../include/jmtx/double/matrices/sparse_conversion.h"
#include "../../../include/jmtx/double/decompositions/incomplete_lu_decomposition.h"
#include "../../../include/jmtx/double/matrices/sparse_conversion.h"
#include "../../../include/jmtx/double/solvers/bicgstab_iteration.h"

enum {PROBLEM_DIMS = (1 << 20), MAX_ITERATIONS = PROBLEM_DIMS, CG_ITERATION_ROUND = 3};

int main()
{
    //  Make the CRS matrix for the 1D Poisson equation
    jmtxd_matrix_cds* mtx = NULL;
    jmtx_result mtx_res;
    //  Problem to solve is d^2/dx^2 (u) - d/dx (u) = -pi^2 sin(pi x) - pi cos(pi x), with u(0) = 0 and u(1) = 0,
    //  on x in (0, 1). Exact solution is u(x) = sin(pi x).
    double* const exact_solution = calloc(PROBLEM_DIMS, sizeof(*exact_solution)); // exact solution of u
    ASSERT(exact_solution);
    double* const forcing_vector = calloc(PROBLEM_DIMS, sizeof(*forcing_vector)); // forcing vector for u (all values are 1)
    ASSERT(forcing_vector);
    double* const iterative_solution = calloc(PROBLEM_DIMS, sizeof(*iterative_solution));
    ASSERT(iterative_solution);
    double* const aux_v1 = calloc(PROBLEM_DIMS, sizeof(*aux_v1));
    ASSERT(aux_v1);
    double* const aux_v2 = calloc(PROBLEM_DIMS, sizeof(*aux_v2));
    ASSERT(aux_v2);
    double* const aux_v3 = calloc(PROBLEM_DIMS, sizeof(*aux_v3));
    ASSERT(aux_v3);
    double* const aux_v4 = calloc(PROBLEM_DIMS, sizeof(*aux_v4));
    ASSERT(aux_v4);
    double* const aux_v5 = calloc(PROBLEM_DIMS, sizeof(*aux_v5));
    ASSERT(aux_v5);
    double* const aux_v6 = calloc(PROBLEM_DIMS, sizeof(*aux_v6));
    ASSERT(aux_v6);
    double* const aux_v7 = calloc(PROBLEM_DIMS, sizeof(*aux_v7));
    ASSERT(aux_v7);
    double* const aux_v8 = calloc(PROBLEM_DIMS, sizeof(*aux_v8));
    ASSERT(aux_v8);
    double* const err_evol = calloc(MAX_ITERATIONS, sizeof(*err_evol));
    ASSERT(err_evol);


    const double dx = 1.0f / (double)(PROBLEM_DIMS - 1);
    for (unsigned i = 0; i < PROBLEM_DIMS; ++i)
    {
        const double x = (double)i / (double)(PROBLEM_DIMS - 1);
        exact_solution[i] = sinf(M_PI * x);
        forcing_vector[i] = (double)(-M_PI * M_PI * sinf(M_PI * x) - M_PI * cosf(M_PI * x));
    }



    MATRIX_TEST_CALL(jmtxds_matrix_cds_new(&mtx, PROBLEM_DIMS, PROBLEM_DIMS, 3, (int32_t[]){-1, 0, +1}, NULL));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    //  Build the matrix
    uint32_t len;
    double* dia = jmtxd_matrix_cds_allocate_diagonal(mtx, 0, &len);
    //  Main diagonal
    dia[0] = 1;
    for (uint_fast32_t i = 1; i < len-1; ++i)
    {
        dia[i] = -2.0f / (dx * dx);
    }
    dia[len-1] = 1;
    dia = jmtxd_matrix_cds_allocate_diagonal(mtx, +1, &len);
    //  Superdiagonal
    dia[0] = 0;
    for (uint_fast32_t i = 1; i < len-1; ++i)
    {
        dia[i] = 1.0f / (dx * dx) + 1.0f / (2.0f * dx);
    }
    dia[len-1] = 0;
    dia = jmtxd_matrix_cds_allocate_diagonal(mtx, -1, &len);
    //  Subdiagonal
    dia[0] = 0;
    for (uint_fast32_t i = 1; i < len-1; ++i)
    {
        dia[i] = 1.0f / (dx * dx) - 1.0f / (2.0f * dx);
    }
    dia[len-1] = 0;

    // print_cdsd_matrix(mtx);
    jmtxd_matrix_cds_vector_multiply(mtx, exact_solution, forcing_vector);
    jmtxd_matrix_crs* newmtx;
    MATRIX_TEST_CALL(jmtxd_convert_cds_to_crs(mtx, &newmtx, NULL));
    jmtxd_matrix_crs* l, *u;
    double ti = omp_get_wtime();


    jmtxd_matrix_ccs *oldu;
    MATRIX_TEST_CALL(jmtxd_decompose_ilu_crs(newmtx, &l, &oldu, NULL));
    MATRIX_TEST_CALL(jmtxd_convert_ccs_to_crs(oldu, &u, NULL));
    jmtxd_matrix_ccs_destroy(oldu);

    double tf = omp_get_wtime();
    printf("ILU took %g\n", tf-ti);
    // exit(0);

//    print_cds_matrix(mtx);
    uint32_t total_iterations = 0;
    double total_time = 0;
    jmtxd_solver_arguments solver_arguments =
            {
            .in_convergence_criterion = 1e-20,
            .in_max_iterations = MAX_ITERATIONS,
            .opt_error_evolution = err_evol,
            };
    for (unsigned i = 0; i < CG_ITERATION_ROUND; ++i)
    {
        const double t0 = omp_get_wtime();
        mtx_res = jmtxd_solve_iterative_pilubicgstab_crs_parallel(
                newmtx, l, u, forcing_vector, iterative_solution, aux_v1, aux_v2, aux_v3, aux_v4, aux_v5, aux_v6, aux_v7, aux_v8, &solver_arguments);
        // mtx_res = jmtxd_solve_iterative_bicgstab_crs(
        //         newmtx, forcing_vector + 1, iterative_solution + 1, aux_v1, aux_v2, aux_v3, aux_v4, aux_v5, aux_v6, &solver_arguments);
        const double t1 = omp_get_wtime();
//        printf("Error evolution:\n");
//        for (uint_fast32_t j = 0; j < solver_arguments.out_last_iteration; ++j)
//        {
//            printf("%"PRIuFAST32": %.10e\n", j, err_evol[j]);
//        }
        ASSERT(mtx_res == JMTX_RESULT_SUCCESS || mtx_res == JMTX_RESULT_NOT_CONVERGED || mtx_res == JMTX_RESULT_STAGNATED);
        double rms_error = 0;
        for (uint_fast32_t j = 0; j < PROBLEM_DIMS; ++j)
        {
            const double local_error = exact_solution[j] - iterative_solution[j];
            rms_error += local_error * local_error;
        }
        rms_error = sqrt(rms_error / (double)PROBLEM_DIMS);
        printf("Solution took %g seconds (%u iterations) for a problem of size %d (outcome: %s), error ratio: %g\nReal RMS error was: %e\n\n", t1 - t0, solver_arguments.out_last_iteration, PROBLEM_DIMS,
               jmtx_result_to_str(mtx_res), solver_arguments.out_last_error, rms_error);

        total_iterations += solver_arguments.out_last_iteration;
        total_time += t1 - t0;
        if (solver_arguments.out_last_error < solver_arguments.in_convergence_criterion)
        {
            break;
        }
    }

    jmtxd_matrix_crs_destroy(l);
    jmtxd_matrix_crs_destroy(u);
    // iterative_solution[0] = 0;
    // iterative_solution[PROBLEM_DIMS - 1] = 0;
    printf("Iterative solution had final residual ratio of %g after %u iterations\n", solver_arguments.out_last_error, total_iterations);

    double rms_error = 0;
    for (unsigned i = 0; i < PROBLEM_DIMS; ++i)
    {
        // printf("%g, %g\n", exact_solution[i], iterative_solution[i]);
//        const double x = (double)i / (double)(PROBLEM_DIMS - 1);
//        printf("u_ex(%g) = %g, u_num(%g) = %g\n", x, exact_solution[i], x, iterative_solution[i]);
        const double local_error = exact_solution[i] - iterative_solution[i];
        rms_error += local_error * local_error;
    }
    rms_error = sqrt(rms_error / PROBLEM_DIMS);
    printf("Iterative solution had final RMS error of %e after %u iterations\n", rms_error, total_iterations);

    MATRIX_TEST_CALL(jmtxds_matrix_crs_destroy(newmtx));
    MATRIX_TEST_CALL(jmtxds_matrix_cds_destroy(mtx));
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
    free(aux_v8);
    return 0;
}