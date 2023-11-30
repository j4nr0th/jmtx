// Automatically generated from source/float/solvers/incomplete_cholesky_decomposition.c on Thu Nov 30 20:33:08 2023
//
// Created by jan on 21.11.2023.
//
#include <assert.h>
#include <math.h>
#include "../../../include/jmtx/double/solvers/incomplete_cholesky_decomposition.h"
#include "../matrices/sparse_row_compressed_internal.h"

jmtx_result jmtxd_incomplete_cholensk_crs(
        const jmtxd_matrix_crs* a, jmtxd_matrix_crs** p_c, const jmtx_allocator_callbacks* allocator_callbacks)
{
    if (!a)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (a->base.type != JMTX_TYPE_CRS)
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

    if (!allocator_callbacks)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }
    else if (!allocator_callbacks->alloc || !allocator_callbacks->free)
    {
        return JMTX_RESULT_NULL_PARAM;
    }

    //  L and U have at most this many entries (in the case that A is already triangular)
    const uint32_t n = a->base.rows;

    jmtxd_matrix_crs* c = NULL;
    jmtx_result res = jmtxd_matrix_crs_copy(a, &c, allocator_callbacks);
    if (res != JMTX_RESULT_SUCCESS)
    {
        return res;
    }
    
    for (uint32_t i = 0; i < n; ++i)
    {
        uint32_t* i_idx = NULL;
        double* i_val = NULL;
        const uint32_t i_cnt = jmtxd_matrix_crs_get_row(c, i, &i_idx, &i_val);
        uint32_t j = 0, p;
        for (p = 0; p < i_cnt && j <= i; ++p, j = i_idx[p])
        {
//            j = i_idx[p];
//            if (j > i)
//            {
//                break;
//            }
            uint32_t* j_idx;
            double* j_val;
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
