#include "../test_common.h"
#include "decompositions/incomplete_cholesky_decomposition.h"
#include "matrices/sparse_multiplication.h"

#include <math.h>
#include <omp.h>

enum
{
    PROBLEM_SIZE_X = 4,
    PROBLEM_SIZE_Y = 4,
    INTERNAL_SIZE_X = PROBLEM_SIZE_X - 2,
    INTERNAL_SIZE_Y = PROBLEM_SIZE_Y - 2,
    PROBLEM_INTERNAL_PTS = INTERNAL_SIZE_X * INTERNAL_SIZE_Y,
    WORK_DIVISIONS = 4,
    MAXIMUM_ITERATIONS = (1 << 10),
};

static unsigned lexicographic_position(const unsigned i, const unsigned j)
{
    return INTERNAL_SIZE_X * i + j;
}

int main()
{
    JMTX_NAME_TYPED(matrix_cds) * mtx;
    jmtx_result mtx_res;
    omp_set_dynamic(1);
    const int proc_count = omp_get_num_procs();
    const int max_threads = omp_get_thread_num();
    printf("OpenMP found %d processors, with a maximum of %d threads\n", proc_count, max_threads);

    double dy = 1.0f / (PROBLEM_SIZE_Y - 1);
    double dx = 1.0f / (PROBLEM_SIZE_X - 1);

    const double rdy2 = 1.0f / (dy * dy);
    const double rdx2 = 1.0f / (dx * dx);

    MATRIX_TEST_CALL(JMTX_NAME_TYPED(matrix_cds_new)(&mtx, PROBLEM_INTERNAL_PTS, PROBLEM_INTERNAL_PTS, 0,
                                                     (const int32_t[]){0}, NULL));

    //  Serial construction
    for (unsigned i = 0; i < INTERNAL_SIZE_Y; ++i)
    {
        for (unsigned j = 0; j < INTERNAL_SIZE_X; ++j)
        {
            //  Point is (DX * j, DY * i)
            unsigned k = 0;
            JMTX_SCALAR_T values[5];
            uint32_t positions[5];
            if (i != 0)
            {
                // There's a bottom boundary
                values[k] = -rdy2;
                positions[k] = lexicographic_position(i - 1, j);
                k += 1;
            }

            if (j != 0)
            {
                // There's a left boundary
                values[k] = -rdx2;
                positions[k] = lexicographic_position(i, j - 1);
                k += 1;
            }

            values[k] = 2 * (rdx2 + rdy2);
            positions[k] = lexicographic_position(i, j);
            k += 1;

            if (j != INTERNAL_SIZE_X - 1)
            {
                // There's a right boundary
                values[k] = -rdx2;
                positions[k] = lexicographic_position(i, j + 1);
                k += 1;
            }

            if (i != INTERNAL_SIZE_Y - 1)
            {
                // There's a top boundary
                values[k] = -rdy2;
                positions[k] = lexicographic_position(i + 1, j);
                k += 1;
            }

            JMTX_NAME_TYPED(matrix_cds_set_row)(mtx, lexicographic_position(i, j), k, values, positions);
        }
    }
    JMTX_NAME_TYPED(print_cds_matrix)(mtx);
    JMTX_NAME_TYPED(matrix_cds) *cholesky = NULL;
    const double t0_decomp = omp_get_wtime();
    MATRIX_TEST_CALL(JMTX_NAME_TYPED(decompose_icho_cds)(mtx, &cholesky, NULL));
    const double t1_decomp = omp_get_wtime();
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);

    JMTX_NAME_TYPED(print_cds_matrix)(cholesky);
    printf("Decomposition took %g seconds and the result: %s\n", t1_decomp - t0_decomp, jmtx_result_to_str(mtx_res));

    JMTX_NAME_TYPED(matrix_cds) *cho_t = NULL;
    MATRIX_TEST_CALL(JMTX_NAME_TYPED(matrix_cds_transpose)(cholesky, &cho_t, NULL));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    JMTX_NAME_TYPED(print_cds_matrix)(cho_t);

    JMTX_NAME_TYPED(matrix_cds) *approx_mtx = NULL;
    MATRIX_TEST_CALL(JMTX_NAME_TYPED(multiply_matrix_cds)(cholesky, cho_t, &approx_mtx, NULL));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);

    JMTX_NAME_TYPED(print_cds_matrix)(mtx);
    JMTX_NAME_TYPED(print_cds_matrix)(approx_mtx);

    JMTX_NAME_TYPED(matrix_cds_destroy)(approx_mtx);
    JMTX_NAME_TYPED(matrix_cds_destroy)(cho_t);
    JMTX_NAME_TYPED(matrix_cds_destroy)(cholesky);
    JMTX_NAME_TYPED(matrix_cds_destroy)(mtx);

    MATRIX_TEST_CALL(JMTX_NAME_TYPED(matrix_cds_new)(&mtx, 16, 16, 0, (int32_t[]){0}, NULL));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    dx = 1.0f / (double)(4 - 1);
    dy = 1.0f / (double)(4 - 1);
    {
        uint32_t dummy_size;
        ASSERT(JMTX_NAME_TYPED(matrix_cds_allocate_zero_diagonal)(mtx, 0, &dummy_size) != NULL);
        ASSERT(JMTX_NAME_TYPED(matrix_cds_allocate_zero_diagonal)(mtx, 1, &dummy_size) != NULL);
        ASSERT(JMTX_NAME_TYPED(matrix_cds_allocate_zero_diagonal)(mtx, 4, &dummy_size) != NULL);
        ASSERT(JMTX_NAME_TYPED(matrix_cds_allocate_zero_diagonal)(mtx, -1, &dummy_size) != NULL);
        ASSERT(JMTX_NAME_TYPED(matrix_cds_allocate_zero_diagonal)(mtx, -4, &dummy_size) != NULL);
    }

    //  Build the whole matrix
    for (unsigned row = 0; row < 4; ++row)
    {
        for (unsigned col = 0; col < 4; ++col)
        {
            const unsigned i = row * 4 + col;

            if (col != 0 && col != 4 - 1 && row != 0 && row != 4 - 1)
            {
                //  Entry corresponds to an interior point
                //  Left point not on boundary
                if (col > 1)
                {
                    MATRIX_TEST_CALL(JMTX_NAME_TYPED(matrix_cds_insert_entry)(mtx, i, i - 1, -1.0f / (dx * dx)));
                }
                //  Right point not on the boundary
                if (col < 4 - 2)
                {
                    MATRIX_TEST_CALL(JMTX_NAME_TYPED(matrix_cds_insert_entry)(mtx, i, i + 1, -1.0 / (dx * dx)));
                }
                //  Middle point (just exists I guess)
                MATRIX_TEST_CALL(
                    JMTX_NAME_TYPED(matrix_cds_insert_entry)(mtx, i, i, 2.0 / (dx * dx) + 2.0 / (dy * dy)));
                //  Top point not on the boundary
                if (row > 1)
                {
                    MATRIX_TEST_CALL(JMTX_NAME_TYPED(matrix_cds_insert_entry)(mtx, i, i - 4, -1.0 / (dy * dy)));
                }
                //  Bottom point not on the boundary
                if (row < 4 - 2)
                {
                    MATRIX_TEST_CALL(JMTX_NAME_TYPED(matrix_cds_insert_entry)(mtx, i, i + 4, -1.0 / (dy * dy)));
                }
            }
            else
            {
                MATRIX_TEST_CALL(JMTX_NAME_TYPED(matrix_cds_insert_entry)(mtx, i, i, 1.0));
            }
        }
    }

    JMTX_NAME_TYPED(print_cds_matrix)(mtx);

    JMTX_NAME_TYPED(matrix_cds) *cho = NULL;
    MATRIX_TEST_CALL(JMTX_NAME_TYPED(decompose_icho_cds)(mtx, &cho, NULL));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);

    JMTX_NAME_TYPED(print_cds_matrix)(cho);

    JMTX_NAME_TYPED(matrix_cds_destroy)(cho);

    JMTX_NAME_TYPED(matrix_cds_destroy)(mtx);

    return 0;
}
