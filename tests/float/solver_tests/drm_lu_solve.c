//
// Created by jan on 2025-08-26.
//
#include "../../../include/jmtx/float/matrices/dense_row_major.h"
#include "../../../include/jmtx/float/decompositions/dense_lu_decomposition.h"
#include "../../../include/jmtx/float/solvers/lu_solving.h"

#include "../test_common.h"

#include <math.h>
#include <stdlib.h>

enum
{
    PROBLEM_DIMS_1 = 4,
    PROBLEM_DIMS_2 = 1,
};

int main()
{
    jmtxf_matrix_drm *mat, *decomp, *a, *b, *sol;
    ASSERT(jmtxf_matrix_drm_new(&mat, PROBLEM_DIMS_1, PROBLEM_DIMS_1, NULL, NULL) == JMTX_RESULT_SUCCESS);
    ASSERT(jmtxf_matrix_drm_new(&a, PROBLEM_DIMS_1, PROBLEM_DIMS_2, NULL, NULL) == JMTX_RESULT_SUCCESS);
    ASSERT(jmtxf_matrix_drm_new(&sol, PROBLEM_DIMS_1, PROBLEM_DIMS_2, NULL, NULL) == JMTX_RESULT_SUCCESS);

    srand(0);

    // fill the matrix "mat" with random garbage (probably won't be singular)
    for (uint32_t i = 0; i < PROBLEM_DIMS_1; ++i)
    {
        for (uint32_t j = 0; j < PROBLEM_DIMS_1; ++j)
        {
            jmtxf_matrix_drm_set_entry(mat, i, j, (float)rand() / (float)RAND_MAX);
        }
    }

    // fill the matrix "a" with random garbage (probably won't be singular)
    for (uint32_t i = 0; i < PROBLEM_DIMS_1; ++i)
    {
        for (uint32_t j = 0; j < PROBLEM_DIMS_2; ++j)
        {
            jmtxf_matrix_drm_set_entry(a, i, j, (float)rand() / (float)RAND_MAX);
        }
    }

    // Compute the matrix "b"
    ASSERT(jmtxf_matrix_drm_new(&b, PROBLEM_DIMS_1, PROBLEM_DIMS_2, NULL, NULL) == JMTX_RESULT_SUCCESS);
    ASSERT(jmtxf_matrix_drm_multiply_matrix(mat, a, b) == JMTX_RESULT_SUCCESS);
    ASSERT(jmtxf_matrix_drm_new(&decomp, PROBLEM_DIMS_1, PROBLEM_DIMS_1, NULL, NULL) == JMTX_RESULT_SUCCESS);

    // Compute the decomposition and solve
    // ASSERT(jmtxf_decompose_lu_pivot_drm(mat, decomp) == JMTX_RESULT_SUCCESS);
    jmtxf_decompose_lu_drm(mat, decomp);

    ASSERT(jmtx_solve_direct_lu_drm_matrix(decomp, b, sol) == JMTX_RESULT_SUCCESS);

    // Check that the computed solution "sol" is very close to the exact solution "a"
    for (uint32_t i = 0; i < PROBLEM_DIMS_1; ++i)
    {
        for (uint32_t j = 0; j < PROBLEM_DIMS_2; ++j)
        {
            ASSERT(fabsf(jmtxf_matrix_drm_get_entry(sol, i, j) - jmtxf_matrix_drm_get_entry(a, i, j)) < 1e-5f);
        }
    }

    // Compute the decomposition with pivoting now
    ASSERT(jmtxf_decompose_lu_pivot_drm(mat, decomp) == JMTX_RESULT_SUCCESS);
    ASSERT(jmtx_solve_direct_lu_drm_matrix(decomp, b, sol) == JMTX_RESULT_SUCCESS);

    // Check that the computed solution "sol" is very close to the exact solution "a"
    print_drm_matrix(a);
    print_drm_matrix(b);
    print_drm_matrix(sol);
    for (uint32_t i = 0; i < PROBLEM_DIMS_1; ++i)
    {
        for (uint32_t j = 0; j < PROBLEM_DIMS_2; ++j)
        {
            printf("Error is %g\n", fabsf(jmtxf_matrix_drm_get_entry(sol, i, j) - jmtxf_matrix_drm_get_entry(a, i, j)));
            ASSERT(fabsf(jmtxf_matrix_drm_get_entry(sol, i, j) - jmtxf_matrix_drm_get_entry(a, i, j)) < 1e-5f);
        }
    }

    return 0;
}
