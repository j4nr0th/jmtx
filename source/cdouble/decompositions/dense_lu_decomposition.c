//
// Created by jan on 4.1.2024.
//
#include "../../../include/jmtx/cdouble/decompositions/dense_lu_decomposition.h"
#include "../matrices/dense_row_major_internal.h"
#include <math.h>
#include <complex.h>

/**
 * Decomposes a matrix into a lower triangular matrix L and upper triangular matrix U, storing the result in the
 * "decomposed" parameter. The diagonal of the lower triangular matrix is 1 by definition, so it is not stored.
 * Since this decomposition does not need pivoting, it means that decomposition will break down if there's a zero on
 * the diagonal or poorly conditioned.
 * @param mtx square matrix to decompose
 * @param decomposed square matrix which receives the decomposition
 */
void jmtxz_decompose_lu_drm(jmtxz_matrix_drm* mtx, jmtxz_matrix_drm* decomposed)
{
    const uint32_t n = mtx->base.rows;
    if (decomposed->rperm)
    {
        decomposed->base.allocator_callbacks.free(decomposed->base.allocator_callbacks.state, decomposed->rperm);
        decomposed->rperm = NULL;
    }
    if (decomposed->permutations)
    {
        decomposed->base.allocator_callbacks.free(decomposed->base.allocator_callbacks.state, decomposed->permutations);
        decomposed->permutations = NULL;
    }
    if (mtx->permutations)
    {
        for (uint32_t i = 0; i < n; ++i)
        {
            memcpy(decomposed->values + n * i, mtx->values + n * mtx->permutations[i], n * sizeof(*decomposed->values));
        }
    }
    else
    {
        memcpy(decomposed->values, mtx->values, n * n * sizeof(*decomposed->values));
    }
    for (uint32_t i = 0; i < n; ++i)
    {
        //  Deal with a row of L
        for (uint_fast32_t j = 0; j < i; ++j)
        {
            _Complex double v = 0;
            //  Row of L
            const _Complex double* li = decomposed->values + n * i;
            //  Column of U
            const _Complex double* uj = decomposed->values + j;
            for (uint_fast32_t k = 0; k < j; ++k)
            {
                v += li[k] * uj[k * n];
            }
            decomposed->values[n * i + j] = (decomposed->values[n * i + j] - v) / uj[n * j];
        }

        //  Deal with a column of U
        for (uint_fast32_t j = 0; j <= i; ++j)
        {
            _Complex double v = 0;
            //  Row of L
            const _Complex double* lj = decomposed->values + n * j;
            //  Column of U
            const _Complex double* ui = decomposed->values + i;
            for (uint_fast32_t k = 0; k < j; ++k)
            {
                v += lj[k] * ui[k * n];
            }
            decomposed->values[n * j + i] = (decomposed->values[n * j + i] - v);
        }
    }
}


/**
 * Decomposes a matrix into a lower triangular matrix L and upper triangular matrix U, storing the result in the
 * "decomposed" parameter. The diagonal of the lower triangular matrix is 1 by definition, so it is not stored.
 * This version uses partial pivoting, so the decomposition should not be permuted further in any way, or its
 * permutations committed. Advantage of using partial pivoting is that the decomposition won't break down if there's
 * a zero on a diagonal as well as be less sensitive to rounding errors.
 * @param mtx square matrix to decompose
 * @param decomposed square matrix which receives the decomposition
 * @return JMTX_RESULT_SUCCESS on success, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxz_decompose_lu_pivot_drm(jmtxz_matrix_drm* mtx, jmtxz_matrix_drm* decomposed)
{
    const uint32_t n = mtx->base.rows;
    if (decomposed->rperm)
    {
        for (uint_fast32_t i = 0; i < n; ++i)
        {
            decomposed->rperm[i] = i;
            decomposed->permutations[i] = i;
        }
    }
    else
    {
        decomposed->permutations = decomposed->base.allocator_callbacks.alloc(decomposed->base.allocator_callbacks.state, sizeof(*decomposed->permutations) * n);
        if (!decomposed->permutations)
        {
            return JMTX_RESULT_BAD_ALLOC;
        }
        decomposed->rperm = decomposed->base.allocator_callbacks.alloc(decomposed->base.allocator_callbacks.state, sizeof(*decomposed->rperm) * n);
        if (!decomposed->rperm)
        {
            decomposed->base.allocator_callbacks.free(decomposed->base.allocator_callbacks.state, decomposed->permutations);
            return JMTX_RESULT_BAD_ALLOC;
        }
        for (uint32_t i = 0; i < n; ++i)
        {
            decomposed->permutations[i] = i;
            decomposed->rperm[i] = i;
        }
    }
    if (mtx->permutations)
    {
        for (uint32_t i = 0; i < n; ++i)
        {
            memcpy(decomposed->values + n * i, mtx->values + n * mtx->permutations[i], n * sizeof(*decomposed->values));
        }
    }
    else
    {
        memcpy(decomposed->values, mtx->values, n * n * sizeof(*decomposed->values));
    }
    for (uint32_t i = 0; i < n; ++i)
    {
        //  Find best pivot in the matrix for the current sub-block
        uint32_t pivot = i;
        for (uint32_t j = i; j < n; ++j)
        {
            if (cabs(decomposed->values[n * decomposed->permutations[j] + i]) > cabs(decomposed->values[n * decomposed->permutations[pivot] + i]))
            {
                pivot = j;
            }
        }
        if (pivot != i && pivot != n)
        {
            (void)jmtxz_matrix_drm_swap_rows(decomposed, i, pivot);
        }

        //  Deal with a row of L
        for (uint_fast32_t j = 0; j < i; ++j)
        {
            _Complex double v = 0;
            //  Row of L
            _Complex double* li = decomposed->values + n * decomposed->permutations[i];
            //  Column of U
            const _Complex double* uj = decomposed->values + j;
            for (uint_fast32_t k = 0; k < j; ++k)
            {
                v += li[k] * uj[decomposed->permutations[k] * n];
            }
            //decomposed->values[n * decomposed->permutations[i] + j]
            li[j] = (li[j] - v) / uj[n * decomposed->permutations[j]];
        }

        //  Deal with a column of U
        for (uint_fast32_t j = 0; j <= i; ++j)
        {
            _Complex double v = 0;
            //  Row of L
            const _Complex double* lj = decomposed->values + n * decomposed->permutations[j];
            //  Column of U
            _Complex double* ui = decomposed->values + i;
            for (uint_fast32_t k = 0; k < j; ++k)
            {
                v += lj[k] * ui[decomposed->permutations[k] * n];
            }
            //decomposed->values[n * decomposed->permutations[j] + i] = (decomposed->values[n * decomposed->permutations[j] + i] - v);
            ui[n * decomposed->permutations[j]] = (ui[n * decomposed->permutations[j]] - v);
        }
    }
    return JMTX_RESULT_SUCCESS;
}
