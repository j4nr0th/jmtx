// Automatically generated from tests/float/solver_tests/incho_test.c on Thu Nov 30 20:00:39 2023
//
// Created by jan on 2.11.2023.
//
#include "../test_common.h"
#include "../../../include/jmtx/double/solvers/incomplete_cholesky_decomposition.h"
#include "../../../include/jmtx/double/matrices/sparse_multiplication.h"
#include "../../../include/jmtx/double/solvers/lu_solving.h"
#include "../../../include/jmtx/double/matrices/sparse_conversion.h"

#include <math.h>
#include <omp.h>
#include "../../../include/jmtx/double/matrices/sparse_row_compressed_safe.h"
#include "../../../include/jmtx/double/matrices/sparse_column_compressed_safe.h"

enum
{
    PROBLEM_SIZE_X = 4, PROBLEM_SIZE_Y = 4,
    INTERNAL_SIZE_X = PROBLEM_SIZE_X - 2,
    INTERNAL_SIZE_Y = PROBLEM_SIZE_Y - 2,
    PROBLEM_INTERNAL_PTS = INTERNAL_SIZE_X * INTERNAL_SIZE_Y,
    WORK_DIVISIONS = 4,
    MAXIMUM_ITERATIONS = (1 << 10),
};

static unsigned lexicographic_position(unsigned i, unsigned j) { return INTERNAL_SIZE_X * i + j; }

static void from_lexicographic(unsigned n, unsigned* pi, unsigned* pj)
{
    *pj = n % INTERNAL_SIZE_X;
    *pi = n / INTERNAL_SIZE_X;
}

int main()
{
    jmtxd_matrix_crs* mtx;
    jmtx_result mtx_res;
    omp_set_dynamic(1);
    const int proc_count = omp_get_num_procs();
    const int max_threads = omp_get_thread_num();
    printf("OpenMP found %d processors, with a maximum of %d threads\n", proc_count, max_threads);

    const double dy = 1.0f / (PROBLEM_SIZE_Y - 1);
    const double dx = 1.0f / (PROBLEM_SIZE_X - 1);

    const double rdy2 = 1.0f / (dy * dy);
    const double rdx2 = 1.0f / (dx * dx);

    MATRIX_TEST_CALL(jmtxds_matrix_crs_new(&mtx, PROBLEM_INTERNAL_PTS, PROBLEM_INTERNAL_PTS, 5 * PROBLEM_INTERNAL_PTS < PROBLEM_INTERNAL_PTS * PROBLEM_INTERNAL_PTS ? 5 * PROBLEM_INTERNAL_PTS : PROBLEM_INTERNAL_PTS * PROBLEM_INTERNAL_PTS, NULL));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    //  Serial construction
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

            jmtxd_matrix_crs_build_row(mtx, lexicographic_position(i, j), k, positions, values);
        }
    }
    jmtxd_matrix_crs* cholesky = NULL;
    const double t0_decomp = omp_get_wtime();
    MATRIX_TEST_CALL(jmtxd_incomplete_cholensk_crs(mtx, &cholesky, NULL));
    const double t1_decomp = omp_get_wtime();
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);

    print_crs_matrix(cholesky);
    printf("Decomposition took %g seconds and the result: %s\n", t1_decomp -
    t0_decomp, jmtx_result_to_str(mtx_res));

    jmtxd_matrix_crs* cpy;
    jmtxd_matrix_ccs* cho_t = NULL;
    MATRIX_TEST_CALL(jmtxd_matrix_crs_copy(cholesky, &cpy, NULL));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);

    MATRIX_TEST_CALL(jmtxd_convert_crs_to_ccs_inplace_transpose(cpy, &cho_t));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    print_ccs_matrix(cho_t);

    jmtxd_matrix_crs* approx_mtx = NULL;
    MATRIX_TEST_CALL(jmtxd_matrix_multiply_crs(cholesky, cho_t, &approx_mtx, NULL));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);

    print_crs_matrix(mtx);
    print_crs_matrix(approx_mtx);

    MATRIX_TEST_CALL(jmtxds_matrix_crs_destroy(approx_mtx));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);

    MATRIX_TEST_CALL(jmtxds_matrix_ccs_destroy(cho_t));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);

    MATRIX_TEST_CALL(jmtxds_matrix_crs_destroy(cholesky));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);

    MATRIX_TEST_CALL(jmtxds_matrix_crs_destroy(mtx));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    return 0;
}
