//
// Created by jan on 23.11.2023.
//
#include "../../../include/jmtx/float/matrices/sparse_row_compressed.h"
#include "../../../source/float/matrices/basic_io.h"
#include "../test_common.h"
#include "../../../include/jmtx/float/matrices/sparse_row_compressed_safe.h"

enum
{
    PROBLEM_SIZE_X = 1 << 2, PROBLEM_SIZE_Y = 1 << 2,
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
    jmtx_matrix_crs* mtx;
    jmtx_result mtx_res;

    const float dy = 1.0f / (PROBLEM_SIZE_Y - 1);
    const float dx = 1.0f / (PROBLEM_SIZE_X - 1);

    const float rdy2 = 1.0f / (dy * dy);
    const float rdx2 = 1.0f / (dx * dx);

    MATRIX_TEST_CALL(jmtxs_matrix_crs_new(&mtx, PROBLEM_INTERNAL_PTS, PROBLEM_INTERNAL_PTS, 5 * PROBLEM_INTERNAL_PTS < PROBLEM_INTERNAL_PTS * PROBLEM_INTERNAL_PTS ? 5 * PROBLEM_INTERNAL_PTS : PROBLEM_INTERNAL_PTS * PROBLEM_INTERNAL_PTS, NULL));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);

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

            jmtx_matrix_crs_build_row(mtx, lexicographic_position(i, j), k, positions, values);
        }
    }

    print_crs_matrix(mtx);
    MATRIX_TEST_CALL(jmtx_matrix_crs_to_file(mtx, "test_output_matrix.txt"));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    MATRIX_TEST_CALL(jmtx_matrix_crs_to_file_explicit(mtx, "test_output_matrix_full.txt"));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);


    MATRIX_TEST_CALL(jmtxs_matrix_crs_destroy(mtx));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    return 0;
}