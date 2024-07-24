//
// Created by jan on 2.11.2023.
//
#include <omp.h>
#include <stdio.h>

#include <math.h>
#include <inttypes.h>
#include "../test_common.h"
#include "../../../include/jmtx/float/matrices/sparse_row_compressed_safe.h"
#include "../../../include/jmtx/float/matrices/sparse_column_compressed_safe.h"
#include "../../../include/jmtx/float/solvers/induced_dimension_reduction_iteration.h"
#include "../../../include/jmtx/float/matrices/sparse_conversion.h"
#include "../../../include/jmtx/float/matrices/dense_row_major.h"

enum
{
    PROBLEM_SIZE_X = 5, PROBLEM_SIZE_Y = 5,
    INTERNAL_SIZE_X = PROBLEM_SIZE_X - 2,
    INTERNAL_SIZE_Y = PROBLEM_SIZE_Y - 2,
    PROBLEM_INTERNAL_PTS = INTERNAL_SIZE_X * INTERNAL_SIZE_Y,
    WORK_DIVISIONS = 4,
    MAXIMUM_ITERATIONS = (PROBLEM_INTERNAL_PTS),
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

    MATRIX_TEST_CALL(jmtxs_matrix_crs_new(&mtx, PROBLEM_INTERNAL_PTS, PROBLEM_INTERNAL_PTS, 5 * PROBLEM_INTERNAL_PTS < PROBLEM_INTERNAL_PTS * PROBLEM_INTERNAL_PTS ? 5 * PROBLEM_INTERNAL_PTS : PROBLEM_INTERNAL_PTS * PROBLEM_INTERNAL_PTS, NULL));
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
//            ASSERT(mtx_res == jmtxs_matrix_crs_beef_check(mtx, &beef_stat));
//            ASSERT(beef_stat == 0xBeef);
        }
    }
    const double t1_serial = omp_get_wtime();

    printf("Serial construction of a %d by %d matrix took %g seconds\n", PROBLEM_INTERNAL_PTS, PROBLEM_INTERNAL_PTS, t1_serial - t0_serial);

    float* const initial_vector = malloc(PROBLEM_INTERNAL_PTS * sizeof(*initial_vector));
    ASSERT(initial_vector != NULL);

    float* const forcing_vector = malloc(PROBLEM_INTERNAL_PTS * sizeof(*forcing_vector));
    ASSERT(forcing_vector != NULL);

    float* const approximate_vector = malloc(PROBLEM_INTERNAL_PTS * sizeof(*approximate_vector));
    ASSERT(approximate_vector != NULL);

    float* const auxiliary_vector = malloc(PROBLEM_INTERNAL_PTS * sizeof(*auxiliary_vector));
    ASSERT(auxiliary_vector != NULL);

    float mag_y = 0;
    for (unsigned i = 0; i < INTERNAL_SIZE_Y; ++i)
    {
        for (unsigned j = 0; j < INTERNAL_SIZE_X; ++j)
        {
            const float v = ((float)(i + 1) * dx) + ((float)(j + 1) * dy);
            initial_vector[i * INTERNAL_SIZE_X + j] = v;
            mag_y += v * v;
        }
    }
    mag_y = sqrtf(mag_y);

    MATRIX_TEST_CALL(jmtxs_matrix_crs_vector_multiply(mtx, initial_vector, forcing_vector));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);


    jmtx_solver_arguments solver_arguments =
            {
            .in_convergence_criterion = 1e-5f,
            .in_max_iterations = MAXIMUM_ITERATIONS,
            };
    memset(approximate_vector, 0, sizeof(*approximate_vector) * PROBLEM_INTERNAL_PTS);
    float* aux_vec1 = calloc(PROBLEM_INTERNAL_PTS, sizeof(*aux_vec1));
    ASSERT(aux_vec1);
    float* aux_vec2 = calloc(PROBLEM_INTERNAL_PTS, sizeof(*aux_vec2));
    ASSERT(aux_vec2);
    float* aux_vec3 = calloc(PROBLEM_INTERNAL_PTS, sizeof(*aux_vec3));
    ASSERT(aux_vec3);
    float* aux_vec4 = calloc(PROBLEM_INTERNAL_PTS, sizeof(*aux_vec4));
    ASSERT(aux_vec4);
    float* aux_vec5 = calloc(PROBLEM_INTERNAL_PTS, sizeof(*aux_vec5));
    ASSERT(aux_vec5);
    float* aux_vec6 = calloc(PROBLEM_INTERNAL_PTS, sizeof(*aux_vec6));
    ASSERT(aux_vec6);
    float* aux_vec7 = calloc(PROBLEM_INTERNAL_PTS, sizeof(*aux_vec7));
    ASSERT(aux_vec7);

    jmtx_matrix_drm* p_mat, *aux_mat1, *aux_mat2, *aux_mat3, *aux_mat4, *aux_mat5;
    MATRIX_TEST_CALL(jmtx_matrix_drm_new(&p_mat, PROBLEM_INTERNAL_PTS, WORK_DIVISIONS, NULL, NULL));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    MATRIX_TEST_CALL(jmtx_matrix_drm_new(&aux_mat1, PROBLEM_INTERNAL_PTS, WORK_DIVISIONS, NULL, NULL));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    MATRIX_TEST_CALL(jmtx_matrix_drm_new(&aux_mat2, PROBLEM_INTERNAL_PTS, WORK_DIVISIONS, NULL, NULL));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    MATRIX_TEST_CALL(jmtx_matrix_drm_new(&aux_mat3, WORK_DIVISIONS, WORK_DIVISIONS, NULL, NULL));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    MATRIX_TEST_CALL(jmtx_matrix_drm_new(&aux_mat4, WORK_DIVISIONS, WORK_DIVISIONS, NULL, NULL));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    MATRIX_TEST_CALL(jmtx_matrix_drm_new(&aux_mat5, WORK_DIVISIONS, WORK_DIVISIONS, NULL, NULL));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);

    //  Initialize matrix P with random orthonormal vectors
    for (uint32_t i = 0; i < WORK_DIVISIONS; ++i)
    {
        float* elems;
        const uint_fast32_t len = jmtx_matrix_drm_get_row(p_mat, i, &elems);
        for (uint32_t j = 0; j < len; ++j)
        {
            double v = drand48();
            elems[j] = (float)(1.0 - 2.0 * v);
        }
        //  Make it orthogonal to other entries
        for (uint32_t j = 0; j < i; ++j)
        {
            float v = 0;
            float* other;
            (void)jmtx_matrix_drm_get_row(p_mat, j, &other);
            for (uint32_t k = 0; k < len; ++k)
            {
                v += elems[k] * other[k];
            }
            for (uint32_t k = 0; k < len; ++k)
            {
                elems[k] -= v * other[k];
            }
        }
        //  Normalize it
        float v = 0;
        for (uint32_t k = 0; k < len; ++k)
        {
            v += elems[k] * elems[k];
        }
        v = sqrtf(v);
        for (uint32_t k = 0; k < len; ++k)
        {
            elems[k] /= v;
        }
    }

    print_drm_matrix(p_mat);
    
    mtx_res = jmtx_solve_iterative_idrs_crs(mtx, forcing_vector, approximate_vector, aux_vec1, aux_vec2, WORK_DIVISIONS,
                                            aux_vec3, aux_vec4, aux_vec5, aux_vec6, aux_vec7, p_mat, aux_mat1, aux_mat2,
                                            aux_mat3, aux_mat4, aux_mat5, &solver_arguments);
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS || mtx_res == JMTX_RESULT_NOT_CONVERGED);

    jmtx_matrix_drm_destroy(aux_mat5);
    jmtx_matrix_drm_destroy(aux_mat4);
    jmtx_matrix_drm_destroy(aux_mat3);
    jmtx_matrix_drm_destroy(aux_mat2);
    jmtx_matrix_drm_destroy(aux_mat1);
    jmtx_matrix_drm_destroy(p_mat);

    free(aux_vec7);
    free(aux_vec6);
    free(aux_vec5);
    free(aux_vec4);
    free(aux_vec3);
    free(aux_vec2);
    free(aux_vec1);
    printf("Solving using ILU took %u iterations, with the final error of %g\n", solver_arguments.out_last_iteration, (double)solver_arguments.out_last_error);

    printf("Comparison of approximate solution vs real solution:\n");
    float rms_err = 0;
    float residual = 0;
    for (unsigned i = 0; i < PROBLEM_INTERNAL_PTS; ++i)
    {
        forcing_vector[i] -= jmtx_matrix_crs_vector_multiply_row(mtx, approximate_vector, i);
        const float err = (initial_vector[i] - approximate_vector[i]) / initial_vector[i];
        rms_err += err * err;
//        printf("Element %u, real: %g, approx: %g, err: %g, residual: %g\n", i, initial_vector[i], approximate_vector[i],
//               err, forcing_vector[i]);
        residual += forcing_vector[i] * forcing_vector[i];
    }
    rms_err = sqrtf(rms_err / PROBLEM_INTERNAL_PTS);
    residual = sqrtf(residual);
    printf("RMS Error: %g, residual/forcing: %g/%g = %g\n", rms_err, residual, mag_y, residual/mag_y);

    free(auxiliary_vector);
    free(approximate_vector);
    free(forcing_vector);
    free(initial_vector);

    MATRIX_TEST_CALL(jmtxs_matrix_crs_destroy(mtx));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    return 0;
}
