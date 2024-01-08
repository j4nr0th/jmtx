// Automatically generated from tests/cfloat/solver_tests/omp_jacobi_test.c on Fri Dec  1 18:48:10 2023
// Automatically generated from tests/cdouble/solver_tests/omp_jacobi_test.c on Fri Dec  1 17:35:45 2023
//
// Created by jan on 22.10.2023.
//
#include <omp.h>
#include "../test_common.h"
#include "../../../include/jmtx/cdouble/matrices/sparse_row_compressed_safe.h"
#include "../../../include/jmtx/cdouble/solvers/jacobi_point_iteration.h"

enum {PROBLEM_DIMS = (1 << 6), MAX_ITERATIONS = (1 << 6)};

int main()
{
    //  Make the CRS matrix for the 1D Poisson equation
    jmtxz_matrix_crs* mtx = NULL;
    jmtx_result mtx_res;
    //  Problem to solve is d^2/dx^2 (u) = 1, with u(0) = 0 and u(1) = 0, on x in (0, 1)
    //  Exact solution is u(x) = x * (x - 1) / 2
    _Complex double exact_solution[PROBLEM_DIMS]; // exact solution of u
    _Complex double forcing_vector[PROBLEM_DIMS]; // forcing vector for u (all values are 1)
    _Complex double iterative_solution[PROBLEM_DIMS] = {0};
    _Complex double aux_v1[PROBLEM_DIMS], aux_v2[PROBLEM_DIMS];

    omp_set_dynamic(1);
    const int proc_count = omp_get_num_procs();
    printf("OpenMP found %d processors\n", proc_count);
    
#pragma omp parallel for shared(exact_solution) default(none) schedule(static)
    for (unsigned i = 0; i < PROBLEM_DIMS; ++i)
    {
        const _Complex double x = (_Complex double)i / (_Complex double)(PROBLEM_DIMS - 1);
        exact_solution[i] = x * (x - 1) / 2;
    }

#pragma omp parallel for shared(forcing_vector) default(none) schedule(static)
    for (unsigned i = 0; i < PROBLEM_DIMS; ++i)
    {
        forcing_vector[i] = 1.0f;
    }

    MATRIX_TEST_CALL(jmtxzs_matrix_crs_new(&mtx, PROBLEM_DIMS - 2, PROBLEM_DIMS - 2, 3 * PROBLEM_DIMS, NULL));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    //  Build the matrix
    {
        const uint32_t indices1[2] = {0, 1};
        const _Complex double values1[2] = {-2.0f, 1.0f};
        const uint32_t indices2[2] = {PROBLEM_DIMS - 4, PROBLEM_DIMS - 3};
        const _Complex double values2[2] = {1.0f, -2.0f};
        forcing_vector[0] += -0.0f;
        forcing_vector[PROBLEM_DIMS - 1] += -0.0f;
        ASSERT(mtx_res == (jmtxzs_matrix_crs_set_row(mtx, 0, 2, indices1, values1)));
        ASSERT(mtx_res == (jmtxzs_matrix_crs_set_row(mtx, PROBLEM_DIMS - 3, 2, indices2, values2)));
    }
    for (unsigned i = 1; i < PROBLEM_DIMS - 3; ++i)
    {
        const uint32_t indices[3] = {i - 1, i, i + 1};
        const _Complex double values[3] = { 1.0f, -2.0f, 1.0f };
        ASSERT(mtx_res == (jmtxzs_matrix_crs_set_row(mtx, i, 3, indices, values)));
    }
//    print_crsz_matrix(mtx);
    jmtxd_solver_arguments solver_arguments =
            {
            .in_convergence_criterion = 1e-4f,
            .in_max_iterations = MAX_ITERATIONS,
            };
    const double t0 = omp_get_wtime();
    mtx_res = jmtxz_solve_iterative_jacobi_crs_parallel(mtx, forcing_vector, iterative_solution + 1, aux_v1, aux_v2, &solver_arguments);
    const double t1 = omp_get_wtime();
    printf("Solution took %g seconds for a problem of size %d\n", t1 - t0, PROBLEM_DIMS);
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS || mtx_res == JMTX_RESULT_NOT_CONVERGED);
    iterative_solution[0] = 0;
    iterative_solution[PROBLEM_DIMS - 1] = 0;
    printf("Iterative solution had final residual ratio of %g after %u iterations\n", solver_arguments.out_last_error, solver_arguments.out_last_iteration);
    const _Complex double dx = 1.0f / (_Complex double)(PROBLEM_DIMS - 1);
#pragma omp parallel for default(none) shared(iterative_solution) shared(dx)
    for (unsigned i = 0; i < PROBLEM_DIMS; ++i)
    {
        iterative_solution[i] *= dx * dx;
    }
    MATRIX_TEST_CALL(jmtxzs_matrix_crs_destroy(mtx));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    return 0;
}