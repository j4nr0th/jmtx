//
// Created by jan on 2.11.2023.
//
#include "../test_common.h"
#include "../../../include/jmtx/cdouble/matrices/sparse_diagonal_compressed_safe.h"
#include "../../../include/jmtx/cdouble/decompositions/incomplete_cholesky_decomposition.h"
#include "../../../include/jmtx/cdouble/matrices/sparse_multiplication.h"
#include "../../../include/jmtx/cdouble/solvers/lu_solving.h"
#include "../../../include/jmtx/cdouble/matrices/sparse_conversion.h"

#include <math.h>
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
    jmtxz_matrix_cds* mtx;
    jmtx_result mtx_res;
    omp_set_dynamic(1);
    const int proc_count = omp_get_num_procs();
    const int max_threads = omp_get_thread_num();
    printf("OpenMP found %d processors, with a maximum of %d threads\n", proc_count, max_threads);

    const double dy = 1.0f / (PROBLEM_SIZE_Y - 1);
    const double dx = 1.0f / (PROBLEM_SIZE_X - 1);

    const double rdy2 = 1.0f / (dy * dy);
    const double rdx2 = 1.0f / (dx * dx);

    MATRIX_TEST_CALL(jmtxzs_matrix_cds_new(&mtx, PROBLEM_INTERNAL_PTS, PROBLEM_INTERNAL_PTS, 0, (const int32_t[]){0}, NULL));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    //  Serial construction
    for (unsigned i = 0; i < INTERNAL_SIZE_Y; ++i)
    {
        for (unsigned j = 0; j < INTERNAL_SIZE_X; ++j)
        {
            //  Point is (DX * j, DY * i)
            unsigned k = 0;
            _Complex double values[5];
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

            jmtxz_matrix_cds_set_row(mtx, lexicographic_position(i, j), k, values, positions);
        }
    }
    print_cdsz_matrix(mtx);
    jmtxz_matrix_cds* cholesky = NULL;
    const double t0_decomp = omp_get_wtime();
    MATRIX_TEST_CALL(jmtxz_decompose_icho_cds(mtx, &cholesky, NULL));
    const double t1_decomp = omp_get_wtime();
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);

    print_cdsz_matrix(cholesky);
    printf("Decomposition took %g seconds and the result: %s\n", t1_decomp -
    t0_decomp, jmtx_result_to_str(mtx_res));

    jmtxz_matrix_cds* cho_t = NULL;
    MATRIX_TEST_CALL(jmtxz_matrix_cds_transpose(cholesky, &cho_t, NULL));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    print_cdsz_matrix(cho_t);

    jmtxz_matrix_cds* approx_mtx = NULL;
    MATRIX_TEST_CALL(jmtxz_multiply_matrix_cds(cholesky, cho_t, &approx_mtx, NULL));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);

    print_cdsz_matrix(mtx);
    print_cdsz_matrix(approx_mtx);

    MATRIX_TEST_CALL(jmtxzs_matrix_cds_destroy(approx_mtx));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);

    MATRIX_TEST_CALL(jmtxzs_matrix_cds_destroy(cho_t));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);

    MATRIX_TEST_CALL(jmtxzs_matrix_cds_destroy(cholesky));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);

    MATRIX_TEST_CALL(jmtxzs_matrix_cds_destroy(mtx));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    return 0;
}
