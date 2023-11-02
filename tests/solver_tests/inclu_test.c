//
// Created by jan on 2.11.2023.
//
#include <omp.h>
#include <stdio.h>
#include "../test_common.h"
#include "../../source/solvers/incomplete_lu_decomposition.h"
#include "../../source/matrices/sparse_multiplication.h"

enum
{
    PROBLEM_SIZE_X = 5, PROBLEM_SIZE_Y = 5,
    INTERNAL_SIZE_X = PROBLEM_SIZE_X - 2,
    INTERNAL_SIZE_Y = PROBLEM_SIZE_Y - 2,
    PROBLEM_INTERNAL_PTS = INTERNAL_SIZE_X * INTERNAL_SIZE_Y,
    WORK_DIVISIONS = 4,
};

static unsigned lexicographic_position(unsigned i, unsigned j) { return INTERNAL_SIZE_X * i + j; }

static void from_lexicographic(unsigned n, unsigned* pi, unsigned* pj)
{
    *pj = n % INTERNAL_SIZE_X;
    *pi = n / INTERNAL_SIZE_X;
}

int main()
{
    jmtx_matrix_crs* mtx;
    jmtx_result mtx_res;
    omp_set_dynamic(1);
    const int proc_count = omp_get_num_procs();
    const int max_threads = omp_get_max_threads();
    printf("OpenMP found %d processors, with a maximum of %d threads\n", proc_count, max_threads);

    const float dy = 1.0f / (PROBLEM_SIZE_Y - 1);
    const float dx = 1.0f / (PROBLEM_SIZE_X - 1);

    const float rdy2 = 1.0f / (dy * dy);
    const float rdx2 = 1.0f / (dx * dx);

    MATRIX_TEST_CALL(jmtx_matrix_crs_new(&mtx, PROBLEM_INTERNAL_PTS, PROBLEM_INTERNAL_PTS, 5 * PROBLEM_INTERNAL_PTS > PROBLEM_INTERNAL_PTS * PROBLEM_INTERNAL_PTS ?: PROBLEM_INTERNAL_PTS * PROBLEM_INTERNAL_PTS, NULL));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    //  Serial construction
    const double t0_serial = omp_get_wtime();
    for (unsigned i = 0; i < INTERNAL_SIZE_Y; ++i)
    {
        for (unsigned j = 0; j < INTERNAL_SIZE_X; ++j)
        {
            //  Point is (DX * j, DY * i)
            unsigned k = 0;
            float values[5];
            uint32_t positions[5];
            if (i != 0)
            {
                // There's a bottom boundary
                values[k] = - rdy2;
                positions[k] = lexicographic_position(i - 1, j);
                k += 1;
            }

            if (j != 0)
            {
                // There's a left boundary
                values[k] = - rdx2;
                positions[k] = lexicographic_position(i, j - 1);
                k += 1;
            }

            values[k] = 2 * (rdx2 + rdy2);
            positions[k] = lexicographic_position(i, j);
            k += 1;

            if (j != INTERNAL_SIZE_X - 1)
            {
                // There's a right boundary
                values[k] = - rdx2;
                positions[k] = lexicographic_position(i, j + 1);
                k += 1;
            }

            if (i != INTERNAL_SIZE_Y - 1)
            {
                // There's a top boundary
                values[k] = - rdy2;
                positions[k] = lexicographic_position(i + 1, j);
                k += 1;
            }

//            MATRIX_TEST_CALL(
            jmtx_matrix_crs_build_row(mtx, lexicographic_position(i, j), k, positions, values);
//                    );
//            ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
//            int beef_stat;
//            ASSERT(mtx_res == jmtx_matrix_crs_beef_check(mtx, &beef_stat));
//            ASSERT(beef_stat == 0xBeef);
        }
    }
    const double t1_serial = omp_get_wtime();

    printf("Serial construction of a %d by %d matrix took %g seconds\n", PROBLEM_INTERNAL_PTS, PROBLEM_INTERNAL_PTS, t1_serial - t0_serial);
    print_crs_matrix(mtx);

    jmtx_matrix_crs* lower = NULL;
    jmtx_matrix_ccs* upper = NULL;
    uint32_t n_iterations;
    float final_error;
    const double t0_decomp = omp_get_wtime();
    MATRIX_TEST_CALL(jmtx_incomplete_lu_crs(mtx, &lower, &upper, 1e-4f, 32, &final_error, &n_iterations, NULL));
    const double t1_decomp = omp_get_wtime();
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS || mtx_res == JMTX_RESULT_NOT_CONVERGED);

    printf("Decomposition took %g seconds for %u iterations, with final error of %g and the result: %s\n", t1_decomp -
    t0_decomp, n_iterations, final_error, jmtx_result_to_str(mtx_res));

    printf("Decomposed lower triangular matrix:\n");
    print_crs_matrix(lower);

    printf("Decomposed upper triangular matrix:\n");
    print_ccs_matrix(upper);

    printf("Product of the incomplete decomposition:\n");
    jmtx_matrix_crs* recon = NULL;
    MATRIX_TEST_CALL(jmtx_matrix_multiply_crs(lower, upper, &recon, NULL));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    print_crs_matrix(recon);


    MATRIX_TEST_CALL(jmtx_matrix_crs_destroy(recon));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    MATRIX_TEST_CALL(jmtx_matrix_crs_destroy(lower));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    MATRIX_TEST_CALL(jmtx_matrix_ccs_destroy(upper));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);

    MATRIX_TEST_CALL(jmtx_matrix_crs_destroy(mtx));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    return 0;
}
