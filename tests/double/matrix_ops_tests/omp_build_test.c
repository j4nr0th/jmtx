// Automatically generated from tests/float/matrix_ops_tests/omp_build_test.c on Fri Dec  1 06:43:09 2023
//
// Created by jan on 30.10.2023.
//
#include <omp.h>
#include "../test_common.h"
#include <inttypes.h>
#include <stdio.h>

#include "../../../include/jmtx/double/matrices/sparse_row_compressed.h"
#include "../../../include/jmtx/double/matrices/sparse_row_compressed_safe.h"

enum
{
    PROBLEM_SIZE_X = 64, PROBLEM_SIZE_Y = 64,
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

static void construct_rows_of_matrix(jmtxd_matrix_crs* mtx, unsigned first, unsigned count, const double rdy2, const double rdx2)
{
    unsigned i, j;
    from_lexicographic(first, &i, &j);
    for (unsigned k = first; k < first + count; ++k)
    {
        unsigned l = 0;
        double values[5];
        uint32_t positions[5];
        if (i != 0)
        {
            // There's a bottom boundary
            values[l] = - rdy2;
            positions[l] = k - INTERNAL_SIZE_X;//lexicographic_position(i - 1, j);
            l += 1;
        }

        if (j != 0)
        {
            // There's a left boundary
            values[l] = - rdx2;
            positions[l] = k - 1;
            l += 1;
        }

        values[l] = 2 * (rdx2 + rdy2);
        positions[l] = k;
        l += 1;

        if (j != INTERNAL_SIZE_X - 1)
        {
            // There's a right boundary
            values[l] = - rdx2;
            positions[l] = k + 1;
            l += 1;
        }

        if (i != INTERNAL_SIZE_Y - 1)
        {
            // There's a top boundary
            values[l] = - rdy2;
            positions[l] = k + INTERNAL_SIZE_X;
            l += 1;
        }

//        MATRIX_TEST_CALL(
        jmtxd_matrix_crs_build_row(mtx, k - first, l, positions, values);
//                );
//        ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
//        int beef_stat;
//        ASSERT(mtx_res == jmtxds_matrix_crs_beef_check(mtx, &beef_stat));
//        ASSERT(beef_stat == 0xBeef);

        if ((j += 1) == INTERNAL_SIZE_X)
        {
            j = 0;
            i += 1;
        }
    }

}

int main()
{
    jmtxd_matrix_crs* mtx;
    jmtx_result mtx_res;
    omp_set_dynamic(1);
    const int proc_count = omp_get_num_procs();
    const int max_threads = omp_get_max_threads();
    printf("OpenMP found %d processors, with a maximum of %d threads\n", proc_count, max_threads);

    const double dy = 1.0f / (PROBLEM_SIZE_Y - 1);
    const double dx = 1.0f / (PROBLEM_SIZE_X - 1);

    const double rdy2 = 1.0f / (dy * dy);
    const double rdx2 = 1.0f / (dx * dx);

    MATRIX_TEST_CALL(jmtxds_matrix_crs_new(&mtx, PROBLEM_INTERNAL_PTS, PROBLEM_INTERNAL_PTS, 5 * PROBLEM_INTERNAL_PTS > PROBLEM_INTERNAL_PTS * PROBLEM_INTERNAL_PTS ?: PROBLEM_INTERNAL_PTS * PROBLEM_INTERNAL_PTS, NULL));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    //  Serial construction
    const double t0_serial = omp_get_wtime();
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

//            MATRIX_TEST_CALL(
                    jmtxd_matrix_crs_build_row(mtx, lexicographic_position(i, j), k, positions, values);
//                    );
//            ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
//            int beef_stat;
//            ASSERT(mtx_res == jmtxds_matrix_crs_beef_check(mtx, &beef_stat));
//            ASSERT(beef_stat == 0xBeef);
        }
    }
    const double t1_serial = omp_get_wtime();

    printf("Serial construction of a %d by %d matrix took %g seconds\n", PROBLEM_INTERNAL_PTS, PROBLEM_INTERNAL_PTS, t1_serial - t0_serial);

    //  Parallel construction

    const double t0_par_outer = omp_get_wtime();

    jmtxd_matrix_crs** const mtx_array = calloc(WORK_DIVISIONS, sizeof(*mtx_array));
    ASSERT(mtx_array);
    unsigned* sizes = calloc(WORK_DIVISIONS + 1, sizeof(*sizes));
    ASSERT(sizes);

    const unsigned single_size = PROBLEM_INTERNAL_PTS / WORK_DIVISIONS;
    unsigned counted = 0;
    for (unsigned i = 0; i < WORK_DIVISIONS + 1; ++i)
    {
        sizes[i] = counted;
        counted += single_size;
    }
    sizes[WORK_DIVISIONS] += PROBLEM_INTERNAL_PTS - (counted - single_size);    //  Add remaining to the last

    for (unsigned i = 0; i < WORK_DIVISIONS; ++i)
    {
        MATRIX_TEST_CALL(jmtxds_matrix_crs_new(mtx_array + i, (sizes[i + 1] - sizes[i]), PROBLEM_INTERNAL_PTS, 5 * (sizes[i + 1] - sizes[i]) < WORK_DIVISIONS * (sizes[i + 1] - sizes[i]) ?: WORK_DIVISIONS * (sizes[i + 1] - sizes[i]),  NULL));
        ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    }

    const double t0_par_inner = omp_get_wtime();
#pragma omp parallel default(none) shared(rdx2, rdy2, sizes, mtx_array)
    {
#pragma omp single
        for (unsigned i = 0; i < WORK_DIVISIONS; ++i)
        {
#pragma omp task default(none) shared(rdx2, rdy2, sizes, mtx_array) firstprivate(i)
            construct_rows_of_matrix(mtx_array[i], sizes[i], sizes[i + 1] - sizes[i], rdy2, rdx2);
        }
    }

    jmtxd_matrix_crs* joined;
    MATRIX_TEST_CALL(jmtxds_matrix_crs_join_vertically(&joined, NULL, WORK_DIVISIONS,
                                                     (const jmtxd_matrix_crs**) mtx_array));
    const double t1_par_inner = omp_get_wtime();
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);


    for (unsigned i = 0; i < WORK_DIVISIONS; ++i)
    {
        MATRIX_TEST_CALL(jmtxds_matrix_crs_destroy(mtx_array[i]));
        ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    }

    free(sizes);
    free(mtx_array);
    const double t1_par_outer = omp_get_wtime();
    const double t_construct = (t1_par_inner - t0_par_inner);
    const double t_overhead = (t1_par_outer - t0_par_outer) - t_construct;
    printf("Constructed matrix in parallel in %g seconds, with %g seconds of additional overhead\n", t_construct, t_overhead);


    //  Checking that parallel joined matrix is the same as the serial
    for (unsigned i = 0; i < PROBLEM_INTERNAL_PTS; ++i)
    {
        for (unsigned j = 0; j < PROBLEM_INTERNAL_PTS; ++j)
        {
            double v_serial, v_joined;
            ASSERT(jmtxds_matrix_crs_get_entry(mtx, i, j, &v_serial) == JMTX_RESULT_SUCCESS);
            ASSERT(jmtxds_matrix_crs_get_entry(joined, i, j, &v_joined) == JMTX_RESULT_SUCCESS);
            ASSERT(v_serial == v_joined);
        }
    }
    MATRIX_TEST_CALL(jmtxds_matrix_crs_destroy(joined));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);

    MATRIX_TEST_CALL(jmtxds_matrix_crs_destroy(mtx));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    return 0;
}
