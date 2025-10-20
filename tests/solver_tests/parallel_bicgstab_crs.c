#include <omp.h>
#include "../test_common.h"
#include <inttypes.h>
#include <math.h>
#include <stdio.h>

#include "decompositions/incomplete_lu_decomposition.h"
#include "matrices/sparse_row_compressed.h"
#include "matrices/sparse_column_compressed.h"
#include "matrices/sparse_conversion.h"
#include "solvers/bicgstab_iteration.h"

enum
{
    PROBLEM_SIZE_X = 64,
    PROBLEM_SIZE_Y = 64,
    INTERNAL_SIZE_X = PROBLEM_SIZE_X - 2,
    INTERNAL_SIZE_Y = PROBLEM_SIZE_Y - 2,
    PROBLEM_INTERNAL_PTS = INTERNAL_SIZE_X * INTERNAL_SIZE_Y,
    WORK_DIVISIONS = 4,
};

static void from_lexicographic(unsigned n, unsigned *pi, unsigned *pj)
{
    *pj = n % INTERNAL_SIZE_X;
    *pi = n / INTERNAL_SIZE_X;
}

/* Same as construct_rows_of_matrix, but inserts zeros into the matrix far away from the diagonal to make ILU
 * decomposition less sparse.
 */
static void construct_rows_of_matrix2(JMTX_NAME_TYPED(matrix_crs) * mtx, unsigned first, unsigned count,
                                      const double rdy2, const double rdx2)
{
    unsigned i, j;
    from_lexicographic(first, &i, &j);
    for (unsigned k = first; k < first + count; ++k)
    {
        jmtx_result t;
        unsigned l = 0;
        JMTX_SCALAR_T values[7];
        uint32_t positions[7];
        if (i != 0)
        {
            // There's a bottom boundary
            values[l] = -rdy2;
            positions[l] = k - INTERNAL_SIZE_X; // lexicographic_position(i - 1, j);
            l += 1;

            values[l] = 0;
            positions[l] = k - INTERNAL_SIZE_X + 1;
            l += 1;
        }

        if (j != 0)
        {
            // There's a left boundary
            values[l] = -rdx2;
            positions[l] = k - 1;
            l += 1;
        }

        values[l] = 2 * (rdx2 + rdy2);
        positions[l] = k;
        l += 1;

        if (j != INTERNAL_SIZE_X - 1)
        {
            // There's a right boundary
            values[l] = -rdx2;
            positions[l] = k + 1;
            l += 1;
        }

        if (i != INTERNAL_SIZE_Y - 1)
        {
            values[l] = 0;
            positions[l] = k + INTERNAL_SIZE_X - 1;
            l += 1;

            // There's a top boundary
            values[l] = -rdy2;
            positions[l] = k + INTERNAL_SIZE_X;
            l += 1;
        }

        //        MATRIX_TEST_CALL(
        ASSERT(t = JMTX_NAME_TYPED(matrix_crs_build_row)(mtx, k - first, l, positions, values) == JMTX_RESULT_SUCCESS);
        // jmtxd_matrix_crs_build_row(mtx, k - first, l, positions, values);
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
    jmtx_result mtx_res;
    omp_set_dynamic(1);
    // omp_set_num_threads(6);
    const int proc_count = omp_get_num_procs();
    const int max_threads = omp_get_max_threads();
    printf("OpenMP found %d processors, with a maximum of %d threads\n", proc_count, max_threads);

    const double dy = 1.0f / (PROBLEM_SIZE_Y - 1);
    const double dx = 1.0f / (PROBLEM_SIZE_X - 1);

    const double rdy2 = 1.0f / (dy * dy);
    const double rdx2 = 1.0f / (dx * dx);

    //  Parallel construction

    const double t0_par_outer = omp_get_wtime();

    JMTX_NAME_TYPED(matrix_crs) **const mtx_array = calloc(WORK_DIVISIONS, sizeof(*mtx_array));
    ASSERT(mtx_array);
    unsigned *sizes = calloc(WORK_DIVISIONS + 1, sizeof(*sizes));
    ASSERT(sizes);

    const unsigned single_size = PROBLEM_INTERNAL_PTS / WORK_DIVISIONS;
    unsigned counted = 0;
    for (unsigned i = 0; i < WORK_DIVISIONS + 1; ++i)
    {
        sizes[i] = counted;
        counted += single_size;
    }
    sizes[WORK_DIVISIONS] += PROBLEM_INTERNAL_PTS - (counted - single_size); //  Add remaining to the last

    for (unsigned i = 0; i < WORK_DIVISIONS; ++i)
    {
        MATRIX_TEST_CALL(
            JMTX_NAME_TYPED(matrix_crs_new)(mtx_array + i, sizes[i + 1] - sizes[i], PROBLEM_INTERNAL_PTS,
                                            5 * (sizes[i + 1] - sizes[i]) < WORK_DIVISIONS * (sizes[i + 1] - sizes[i])
                                                ?: WORK_DIVISIONS * (sizes[i + 1] - sizes[i]),
                                            NULL));
        ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    }

    const double t0_par_inner = omp_get_wtime();
#pragma omp parallel default(none) shared(rdx2, rdy2, sizes, mtx_array)
    {
#pragma omp single
        for (unsigned i = 0; i < WORK_DIVISIONS; ++i)
        {
#pragma omp task default(none) shared(rdx2, rdy2, sizes, mtx_array) firstprivate(i)
            construct_rows_of_matrix2(mtx_array[i], sizes[i], sizes[i + 1] - sizes[i], rdy2, rdx2);
        }
    }

    JMTX_NAME_TYPED(matrix_crs) * joined;
    MATRIX_TEST_CALL(JMTX_NAME_TYPED(matrix_crs_join_vertically)(&joined, NULL, WORK_DIVISIONS,
                                                                 (const JMTX_NAME_TYPED(matrix_crs) **)mtx_array));
    const double t1_par_inner = omp_get_wtime();
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);

    for (unsigned i = 0; i < WORK_DIVISIONS; ++i)
    {
        JMTX_NAME_TYPED(matrix_crs_destroy)(mtx_array[i]);
    }

    free(sizes);
    free(mtx_array);
    const double t1_par_outer = omp_get_wtime();
    const double t_construct = t1_par_inner - t0_par_inner;
    const double t_overhead = t1_par_outer - t0_par_outer - t_construct;
    printf("Constructed matrix in parallel in %g seconds, with %g seconds of additional overhead\n", t_construct,
           t_overhead);

    JMTX_NAME_TYPED(matrix_crs) * l, *u;
    JMTX_NAME_TYPED(matrix_ccs) * u_ccs;
    // Get the ILU decomposition
    MATRIX_TEST_CALL(JMTX_NAME_TYPED(decompose_ilu_crs)(joined, &l, &u_ccs, NULL));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    MATRIX_TEST_CALL(JMTX_NAME_TYPED(convert_ccs_to_crs)(u_ccs, &u, NULL));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    JMTX_NAME_TYPED(matrix_ccs_destroy)(u_ccs);
    u_ccs = NULL;

    JMTX_SCALAR_T *lhs = calloc(PROBLEM_INTERNAL_PTS, sizeof(*lhs));
    ASSERT(lhs);
    JMTX_SCALAR_T *rhs = calloc(PROBLEM_INTERNAL_PTS, sizeof(*rhs));
    ASSERT(rhs);
    JMTX_SCALAR_T *x = calloc(PROBLEM_INTERNAL_PTS, sizeof(*x));
    ASSERT(x);

    JMTX_SCALAR_T *aux1 = calloc(PROBLEM_INTERNAL_PTS, sizeof(*aux1));
    ASSERT(aux1);
    JMTX_SCALAR_T *aux2 = calloc(PROBLEM_INTERNAL_PTS, sizeof(*aux2));
    ASSERT(aux2);
    JMTX_SCALAR_T *aux3 = calloc(PROBLEM_INTERNAL_PTS, sizeof(*aux3));
    ASSERT(aux3);
    JMTX_SCALAR_T *aux4 = calloc(PROBLEM_INTERNAL_PTS, sizeof(*aux4));
    ASSERT(aux4);
    JMTX_SCALAR_T *aux5 = calloc(PROBLEM_INTERNAL_PTS, sizeof(*aux5));
    ASSERT(aux5);
    JMTX_SCALAR_T *aux6 = calloc(PROBLEM_INTERNAL_PTS, sizeof(*aux6));
    ASSERT(aux6);
    JMTX_SCALAR_T *aux7 = calloc(PROBLEM_INTERNAL_PTS, sizeof(*aux7));
    ASSERT(aux7);
    JMTX_SCALAR_T *aux8 = calloc(PROBLEM_INTERNAL_PTS, sizeof(*aux8));
    ASSERT(aux8);

    for (unsigned i = 0; i < PROBLEM_INTERNAL_PTS; ++i)
    {
        lhs[i] = sin((double)(i % PROBLEM_SIZE_X) / (double)(PROBLEM_SIZE_X - 1)) *
                 cos((double)(i / PROBLEM_SIZE_X) / (double)(PROBLEM_SIZE_Y - 1));
    }

    JMTX_NAME_TYPED(matrix_crs_vector_multiply)(joined, lhs, rhs);

    JMTX_NAME_TYPED(solver_arguments)
    args = {.in_convergence_criterion = 1e-9,
            .in_max_iterations = 256,
            .out_last_error = 1,
            .out_last_iteration = 0,
            .opt_error_evolution = NULL};
    const double t0 = omp_get_wtime();
    MATRIX_TEST_CALL(JMTX_NAME_TYPED(solve_iterative_pilubicgstab_crs)(joined, l, u, rhs, x, aux1, aux2, aux3, aux4,
                                                                       aux5, aux6, aux7, aux8, &args));
    const double t1 = omp_get_wtime();
    printf("Iterative solver took %g seconds for %u iterations, with the final criterion at %g\n", t1 - t0,
           args.out_last_iteration, args.out_last_error);

    free(aux8);
    free(aux7);
    free(aux6);
    free(aux5);
    free(aux4);
    free(aux3);
    free(aux2);
    free(aux1);
    free(x);
    free(rhs);
    free(lhs);
    JMTX_NAME_TYPED(matrix_crs_destroy)(u);
    JMTX_NAME_TYPED(matrix_crs_destroy)(l);
    JMTX_NAME_TYPED(matrix_crs_destroy)(joined);

    return 0;
}
