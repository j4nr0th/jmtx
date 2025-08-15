// Automatically generated from source/float/solvers/incomplete_cholesky_decomposition.c on Fri Dec  1 06:43:01 2023
//
// Created by jan on 21.11.2023.
//
#include "../../../include/jmtx/double/decompositions/incomplete_cholesky_decomposition.h"
#include "../matrices/sparse_diagonal_compressed_internal.h"
#include "../matrices/sparse_row_compressed_internal.h"
#include <assert.h>
#include <math.h>

jmtx_result jmtxd_decompose_icho_crs(const jmtxd_matrix_crs *a, jmtxd_matrix_crs **p_c,
                                     const jmtx_allocator_callbacks *allocator_callbacks)
{
    if (!a)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (a->base.type != JMTXD_TYPE_CRS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (a->base.cols != a->base.rows)
    {
        //  Only doing square matrices
        return JMTX_RESULT_BAD_MATRIX;
    }
    if (!p_c)
    {
        return JMTX_RESULT_NULL_PARAM;
    }

    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }
    else if (allocator_callbacks->alloc == NULL || allocator_callbacks->free == NULL)
    {
        return JMTX_RESULT_NULL_PARAM;
    }

    //  L and U have at most this many entries (in the case that A is already triangular)
    const uint32_t n = a->base.rows;

    jmtxd_matrix_crs *c = NULL;
    jmtx_result res = jmtxd_matrix_crs_copy(a, &c, allocator_callbacks);
    if (res != JMTX_RESULT_SUCCESS)
    {
        return res;
    }

    for (uint32_t i = 0; i < n; ++i)
    {
        uint32_t *i_idx = NULL;
        double *i_val = NULL;
        const uint32_t i_cnt = jmtxd_matrix_crs_get_row(c, i, &i_idx, &i_val);
        uint32_t j = 0, p;
        for (p = 0; p < i_cnt && j <= i; ++p)
        {
            j = i_idx[p];
            //            if (j > i)
            //            {
            //                break;
            //            }
            uint32_t *j_idx;
            double *j_val;
            const uint32_t j_cnt = jmtxd_matrix_crs_get_row(c, j, &j_idx, &j_val);
            double v = 0.0f;
            uint32_t ki, kj;
            for (ki = 0, kj = 0; ki < i_cnt && kj < j_cnt && i_idx[ki] < j && j_idx[kj] < j;)
            {
                if (i_idx[ki] == j_idx[kj])
                {
                    v += i_val[ki] * j_val[kj];
                    ki += 1;
                    kj += 1;
                }
                else if (i_idx[ki] > j_idx[kj])
                {
                    kj += 1;
                }
                else // if(i_idx[ki] < j_idx[kj])
                {
                    ki += 1;
                }
            }
            assert(j_idx[kj] <= j);
            while (j_idx[kj] < j)
            {
                kj += 1;
            }
            //  Zero on diagonal should not happen because it would've been encountered by now
            //            if (i_idx[kj] != j)
            //            {
            //                jmtxd_matrix_crs_destroy(c);
            //                return JMTX_RESULT_BAD_MATRIX;
            //            }
            assert(j_idx[kj] == j);

            if (i != j)
            {
                const double l_ij = (i_val[p] - v) / j_val[kj];
                i_val[p] = l_ij;
            }
            else
            {
                const double l_ij = sqrt((i_val[p] - v));
                i_val[p] = l_ij;
                p += 1;
                break;
            }
        }
        if (j != i)
        {
            //  There was no diagonal entry!
            jmtxd_matrix_crs_destroy(c);
            return JMTX_RESULT_BAD_MATRIX;
        }

        //  Zero the rest of the row out
        while (p < i_cnt)
        {
            i_val[p] = 0;
            p += 1;
        }
    }

    jmtxd_matrix_crs_remove_zeros(c);
    *p_c = c;

    return JMTX_RESULT_SUCCESS;
}

/**
 * Uses relations for Cholesky decomposition to compute an approximate decomposition with C' such that the matrix
 * C'C'^T has the same sparsity as the starting matrix. This decomposition can be used as a preconditioner or directly.
 * As with full Cholesky decomposition, the matrix to be decomposed must be SPD. C' is computed so that it has the same
 * sparsity pattern as the matrix A. C' is an lower triangular matrix.
 *
 * @param a matrix to decompose
 * @param p_c pointer which receives the resulting C' CRS matrix
 * @param allocator_callbacks Pointer to allocators to use for allocating C', and auxiliary memory. If NULL, malloc
 * and free are used.
 * @return JMTX_RESULT_SUCCESS if successfully converged in to tolerance in under max iterations,
 * JMTX_RESULT_NOT_CONVERGED if convergence was not achieved in number of specified iterations,
 * other jmtx_result values on other failures.
 */
jmtx_result jmtxd_decompose_icho_cds(const jmtxd_matrix_cds *a, jmtxd_matrix_cds **p_c,
                                     const jmtx_allocator_callbacks *allocator_callbacks)

{
    if (!a)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (a->base.type != JMTXD_TYPE_CDS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (a->base.cols != a->base.rows)
    {
        //  Only doing square matrices
        return JMTX_RESULT_BAD_MATRIX;
    }
    if (!p_c)
    {
        return JMTX_RESULT_NULL_PARAM;
    }

    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }
    else if (allocator_callbacks->alloc == NULL || allocator_callbacks->free == NULL)
    {
        return JMTX_RESULT_NULL_PARAM;
    }

    //  L and U have at most this many entries (in the case that A is already triangular)
    const uint32_t n = a->base.rows;
    const uint_fast32_t max_per_row = jmtxd_matrix_cds_diagonal_count(a);
    uint32_t *const i_indices =
        allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*i_indices) * max_per_row);
    if (!i_indices)
    {
        return JMTX_RESULT_BAD_ALLOC;
    }
    double *const i_values = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*i_values) * max_per_row);
    if (!i_values)
    {
        allocator_callbacks->free(allocator_callbacks->state, i_indices);
        return JMTX_RESULT_BAD_ALLOC;
    }

    uint32_t *const j_indices =
        allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*j_indices) * max_per_row);
    if (!j_indices)
    {
        allocator_callbacks->free(allocator_callbacks->state, i_values);
        allocator_callbacks->free(allocator_callbacks->state, i_indices);
        return JMTX_RESULT_BAD_ALLOC;
    }
    double *const j_values = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*j_values) * max_per_row);
    if (!j_values)
    {
        allocator_callbacks->free(allocator_callbacks->state, j_indices);
        allocator_callbacks->free(allocator_callbacks->state, i_values);
        allocator_callbacks->free(allocator_callbacks->state, i_indices);
        return JMTX_RESULT_BAD_ALLOC;
    }

    jmtxd_matrix_cds *c = NULL;
    jmtx_result res = jmtxd_matrix_cds_copy(a, &c, allocator_callbacks);
    if (res != JMTX_RESULT_SUCCESS)
    {
        allocator_callbacks->free(allocator_callbacks->state, j_values);
        allocator_callbacks->free(allocator_callbacks->state, j_indices);
        allocator_callbacks->free(allocator_callbacks->state, i_values);
        allocator_callbacks->free(allocator_callbacks->state, i_indices);
        return res;
    }

    for (uint_fast32_t i = 0; i < n; ++i)
    {
        const uint_fast32_t i_cnt = jmtxd_matrix_cds_get_row(c, i, max_per_row, i_values, i_indices);
        uint_fast32_t j = 0, p;
        for (p = 0; p < i_cnt && j <= i; ++p)
        {
            if (i_values[p] == 0)
            {
                continue;
            }
            j = i_indices[p];
            const uint32_t j_cnt = jmtxd_matrix_cds_get_row(c, j, max_per_row, j_values, j_indices);
            double v = 0.0f;
            uint32_t ki, kj;
            for (ki = 0, kj = 0; ki < i_cnt && kj < j_cnt && i_indices[ki] < j && j_indices[kj] < j;)
            {
                if (i_indices[ki] == j_indices[kj])
                {
                    v += i_values[ki] * j_values[kj];
                    ki += 1;
                    kj += 1;
                }
                else if (i_indices[ki] > j_indices[kj])
                {
                    kj += 1;
                }
                else // if(i_idx[ki] < j_idx[kj])
                {
                    ki += 1;
                }
            }
            assert(j_indices[kj] <= j);
            while (j_indices[kj] < j)
            {
                kj += 1;
            }
            //  Zero on diagonal should not happen because it would've been encountered by now
            //            if (i_idx[kj] != j)
            //            {
            //                jmtxd_matrix_crs_destroy(c);
            //                return JMTX_RESULT_BAD_MATRIX;
            //            }
            assert(j_indices[kj] == j);

            if (i != j)
            {
                const double l_ij = (i_values[p] - v) / j_values[kj];
                i_values[p] = l_ij;
                jmtxd_matrix_cds_set_entry(c, i, i_indices[p], l_ij);
            }
            else
            {
                const double l_ij = sqrt((i_values[p] - v));
                i_values[p] = l_ij;
                jmtxd_matrix_cds_set_entry(c, i, i_indices[p], l_ij);
                p += 1;
                break;
            }
        }
        if (j != i)
        {
            //  There was no diagonal entry!
            allocator_callbacks->free(allocator_callbacks->state, j_values);
            allocator_callbacks->free(allocator_callbacks->state, j_indices);
            allocator_callbacks->free(allocator_callbacks->state, i_values);
            allocator_callbacks->free(allocator_callbacks->state, i_indices);
            jmtxd_matrix_cds_destroy(c);
            return JMTX_RESULT_BAD_MATRIX;
        }

        //        uint_fast32_t k = p;
        //  Zero the rest of the row out
        while (p < i_cnt)
        {
            jmtxd_matrix_cds_set_entry(c, i, i_indices[p], 0);
            i_values[p] = 0;
            p += 1;
        }

        //        res = jmtxd_matrix_cds_set_row(c, i, p, i_values, i_indices);
        //        if (res != JMTX_RESULT_SUCCESS)
        //        {
        //            allocator_callbacks->free(allocator_callbacks->state, j_values);
        //            allocator_callbacks->free(allocator_callbacks->state, j_indices);
        //            allocator_callbacks->free(allocator_callbacks->state, i_values);
        //            allocator_callbacks->free(allocator_callbacks->state, i_indices);
        //            jmtxd_matrix_cds_destroy(c);
        //            return res;
        //        }
    }
    allocator_callbacks->free(allocator_callbacks->state, j_values);
    allocator_callbacks->free(allocator_callbacks->state, j_indices);
    allocator_callbacks->free(allocator_callbacks->state, i_values);
    allocator_callbacks->free(allocator_callbacks->state, i_indices);

    //  Remove superdiagonals
    for (uint_fast32_t i = 0; i < c->super_diagonals.count; ++i)
    {
        c->base.allocator_callbacks.free(c->base.allocator_callbacks.state, c->super_diagonals.diagonals[i]);
#ifndef NDEBUG
        c->super_diagonals.diagonals[i] = (void *)0xCCCCCCCCCCCCCCCC;
#endif
    }
    c->super_diagonals.count = 0;

    *p_c = c;

    return JMTX_RESULT_SUCCESS;
}
