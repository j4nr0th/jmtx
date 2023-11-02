//
// Created by jan on 22.10.2023.
//
#include <omp.h>
#include <inttypes.h>
#include "../test_common.h"
#include "../../source/solvers/gauss_seidel_iteration.h"

enum {PROBLEM_DIMS = (1 << 8), MAX_ITERATIONS = (1 << 14)};

int main()
{
    //  Make the CRS matrix for the 1D Poisson equation
    jmtx_matrix_crs* mtx = NULL;
    jmtx_result mtx_res;
    //  Problem to solve is d^2/dx^2 (u) = 1, with u(0) = 0 and u(1) = 0, on x in (0, 1)
    //  Exact solution is u(x) = x * (x - 1) / 2
    float exact_solution[PROBLEM_DIMS]; // exact solution of u
    float forcing_vector[PROBLEM_DIMS]; // forcing vector for u (all elements are 1)
    float iterative_solution[PROBLEM_DIMS] = {0};
    float aux_v1[PROBLEM_DIMS];

    omp_set_dynamic(1);
    const int proc_count = omp_get_num_procs();
    printf("OpenMP found %d processors\n", proc_count);
    
#pragma omp parallel for shared(exact_solution) default(none) schedule(static)
    for (unsigned i = 0; i < PROBLEM_DIMS; ++i)
    {
        const float x = (float)i / (float)(PROBLEM_DIMS - 1);
        exact_solution[i] = x * (x - 1) / 2;
    }

#pragma omp parallel for shared(forcing_vector) default(none) schedule(static)
    for (unsigned i = 0; i < PROBLEM_DIMS; ++i)
    {
        forcing_vector[i] = 1.0f;
    }

    MATRIX_TEST_CALL(jmtx_matrix_crs_new(&mtx, PROBLEM_DIMS - 2, PROBLEM_DIMS - 2, 3 * PROBLEM_DIMS, NULL));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    //  Build the matrix
    {
        const uint32_t indices1[2] = {0, 1};
        const float values1[2] = {-2.0f, 1.0f};
        const uint32_t indices2[2] = {PROBLEM_DIMS - 4, PROBLEM_DIMS - 3};
        const float values2[2] = {1.0f, -2.0f};
        forcing_vector[0] += -0.0f;
        forcing_vector[PROBLEM_DIMS - 1] += -0.0f;
        ASSERT(mtx_res == (jmtx_matrix_crs_set_row(mtx, 0, 2, indices1, values1)));
        ASSERT(mtx_res == (jmtx_matrix_crs_set_row(mtx, PROBLEM_DIMS - 3, 2, indices2, values2)));
    }
    for (unsigned i = 1; i < PROBLEM_DIMS - 3; ++i)
    {
        const uint32_t indices[3] = {i - 1, i, i + 1};
        const float values[3] = { 1.0f, -2.0f, 1.0f };
        ASSERT(mtx_res == (jmtx_matrix_crs_set_row(mtx, i, 3, indices, values)));
    }
//    print_crs_matrix(mtx);
    uint32_t iterations = 0;
    float final_err;
    const double t0 = omp_get_wtime();
    mtx_res = jmtx_gauss_seidel_crs_parallel(
            mtx, forcing_vector, iterative_solution + 1, 1e-4f, MAX_ITERATIONS, &iterations, NULL, &final_err, aux_v1);
    const double t1 = omp_get_wtime();
    printf("Solution took %g seconds for a problem of size %d\n", t1 - t0, PROBLEM_DIMS);
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS || mtx_res == JMTX_RESULT_NOT_CONVERGED);
    iterative_solution[0] = 0;
    iterative_solution[PROBLEM_DIMS - 1] = 0;
    printf("Iterative solution had final residual ratio of %g after %u iterations\n", final_err, iterations);
    const float dx = 1.0f / (float)(PROBLEM_DIMS - 1);
#pragma omp parallel for default(none) shared(iterative_solution) shared(dx)
    for (unsigned i = 0; i < PROBLEM_DIMS; ++i)
    {
        iterative_solution[i] *= dx * dx;
    }
//    for (unsigned i = 0; i < PROBLEM_DIMS; ++i)
//    {
//        const float x = (float)i / (float)(PROBLEM_DIMS - 1);
//        printf("u_ex(%g) = %g, u_num(%g) = %g\n", x, exact_solution[i], x, iterative_solution[i]);
//    }
    MATRIX_TEST_CALL(jmtx_matrix_crs_destroy(mtx));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    return 0;
}