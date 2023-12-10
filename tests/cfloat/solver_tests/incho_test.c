// Automatically generated from tests/float/solver_tests/incho_test.c on Fri Dec  1 17:35:45 2023
//
// Created by jan on 2.11.2023.
//
#include "../test_common.h"
#include "../../../include/jmtx/cfloat/matrices/sparse_row_compressed_safe.h"
#include "../../../include/jmtx/cfloat/matrices/sparse_column_compressed_safe.h"
#include "../../../include/jmtx/cfloat/solvers/incomplete_cholesky_decomposition.h"
#include "../../../include/jmtx/cfloat/matrices/sparse_multiplication.h"
#include "../../../include/jmtx/cfloat/matrices/sparse_conversion.h"

#include <omp.h>

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
    jmtxc_matrix_crs* mtx;
    jmtx_result mtx_res;
    omp_set_dynamic(1);
    const int proc_count = omp_get_num_procs();
    const int max_threads = omp_get_thread_num();
    printf("OpenMP found %d processors, with a maximum of %d threads\n", proc_count, max_threads);

    const _Complex float dy = 1.0f / (PROBLEM_SIZE_Y - 1);
    const _Complex float dx = 1.0f / (PROBLEM_SIZE_X - 1);

    const _Complex float rdy2 = 1.0f / (dy * dy);
    const _Complex float rdx2 = 1.0f / (dx * dx);

    MATRIX_TEST_CALL(jmtxcs_matrix_crs_new(&mtx, PROBLEM_INTERNAL_PTS, PROBLEM_INTERNAL_PTS, 5 * PROBLEM_INTERNAL_PTS < PROBLEM_INTERNAL_PTS * PROBLEM_INTERNAL_PTS ? 5 * PROBLEM_INTERNAL_PTS : PROBLEM_INTERNAL_PTS * PROBLEM_INTERNAL_PTS, NULL));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    //  Serial construction
    for (unsigned i = 0; i < INTERNAL_SIZE_Y; ++i)
    {
        for (unsigned j = 0; j < INTERNAL_SIZE_X; ++j)
        {
            //  Point is (DX * j, DY * i)
            unsigned k = 0;
            _Complex float values[5];
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

            jmtxc_matrix_crs_build_row(mtx, lexicographic_position(i, j), k, positions, values);
        }
    }
    jmtxc_matrix_crs* cholesky = NULL;
    const double t0_decomp = omp_get_wtime();
    MATRIX_TEST_CALL(jmtxc_incomplete_cholesky_crs(mtx, &cholesky, NULL));
    const double t1_decomp = omp_get_wtime();
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);

    print_crs_matrix(cholesky);
    printf("Decomposition took %g seconds and the result: %s\n", t1_decomp -
    t0_decomp, jmtx_result_to_str(mtx_res));

    jmtxc_matrix_crs* cpy;
    jmtxc_matrix_ccs* cho_t = NULL;
    MATRIX_TEST_CALL(jmtxc_matrix_crs_copy(cholesky, &cpy, NULL));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);

    MATRIX_TEST_CALL(jmtxc_convert_crs_to_ccs_inplace_transpose(cpy, &cho_t));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    print_ccs_matrix(cho_t);

    jmtxc_matrix_crs* approx_mtx = NULL;
    MATRIX_TEST_CALL(jmtxc_matrix_multiply_crs(cholesky, cho_t, &approx_mtx, NULL));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);

    print_crs_matrix(mtx);
    print_crs_matrix(approx_mtx);

    MATRIX_TEST_CALL(jmtxcs_matrix_crs_destroy(approx_mtx));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);

    MATRIX_TEST_CALL(jmtxcs_matrix_ccs_destroy(cho_t));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);

    MATRIX_TEST_CALL(jmtxcs_matrix_crs_destroy(cholesky));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);

    MATRIX_TEST_CALL(jmtxcs_matrix_crs_destroy(mtx));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    return 0;
}
