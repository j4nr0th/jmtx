// Automatically generated from source/float/solvers/incomplete_lu_decomposition.c on Sun Dec 17 18:02:30 2023
//
// Created by jan on 2.11.2023.
//

#include <assert.h>
#include <math.h>
#include "../../../include/jmtx/cdouble/decompositions/incomplete_lu_decomposition.h"
#include "../matrices/sparse_row_compressed_internal.h"
#include "../matrices/sparse_column_compressed_internal.h"
#include "../../../tests/cdouble/test_common.h"
#include "../../../include/jmtx/cdouble/matrices/sparse_multiplication.h"

/**
 * Uses relations for LU decomposition to compute an approximate decomposition with L' and U' such that the matrix
 * L'U' has the same sparsity as the starting matrix. This decomposition can be used as a preconditioner or directly.
 *
 * For a symmetric SPD matrix, an incomplete Cholesky factorization is used instead, which exploits the symmetry of the
 * matrix to give the decomposition in the form of C' C'^T = A, where C' has the same sparsity pattern as the top of
 * the matrix A.
 *
 * @param a matrix to decompose
 * @param p_l pointer which receives the resulting L' CRS matrix
 * @param p_u pointer which receives the resulting U' CCS matrix
 * @param allocator_callbacks Pointer to allocators to use for allocating L', U', and auxiliary memory. If NULL, malloc
 * and free are used.
 * @return JMTX_RESULT_SUCCESS if successfully converged in to tolerance in under max iterations,
 * JMTX_RESULT_NOT_CONVERGED if convergence was not achieved in number of specified iterations,
 * other jmtx_result values on other failures.
 */
jmtx_result jmtxzs_decompose_ilu_crs(
        const jmtxz_matrix_crs* a, jmtxz_matrix_crs** p_l, jmtxz_matrix_ccs** p_u,
        const jmtx_allocator_callbacks* allocator_callbacks)
{
    if (!a)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (a->base.type != JMTXZ_TYPE_CRS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (a->base.cols != a->base.rows)
    {
        //  Only doing square matrices
        return JMTX_RESULT_BAD_MATRIX;
    }
    if (!p_l)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!p_u)
    {
        return JMTX_RESULT_NULL_PARAM;
    }

    if (allocator_callbacks && (!allocator_callbacks->alloc || !allocator_callbacks->free))
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    return jmtxz_decompose_ilu_crs(a, p_l, p_u, allocator_callbacks);
}


jmtx_result jmtxz_decompose_ldu_split_crs(const jmtxz_matrix_crs* a, jmtxz_matrix_crs** p_l, jmtxz_matrix_ccs** p_u,
    const jmtx_allocator_callbacks* allocator_callbacks)
{
    if (!allocator_callbacks)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }

    jmtxz_matrix_crs* l;
    jmtxz_matrix_ccs* u;



    uint32_t cols = a->base.cols;
    uint32_t* count_u = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*count_u) * cols);
    if (!count_u)
    {
        return JMTX_RESULT_BAD_ALLOC;
    }
    uint32_t rows = a->base.rows;
    uint32_t total_count_l = 0, total_count_u = 0;
    memset(count_u, 0, sizeof(*count_u) * cols);


    for (uint32_t row = 0; row < rows; ++row)
    {
        uint32_t* pcols;
        _Complex double* pvals;
        uint32_t ncols = jmtxz_matrix_crs_get_row(a, row, &pcols, &pvals);
        for (uint32_t n = 0; n < ncols; ++n)
        {
            uint32_t col = pcols[n];
            count_u[col] += (col >= row);
            total_count_u += (col >= row);
            total_count_l += (col <= row);
        }
    }

    jmtx_result res = jmtxz_matrix_crs_new(&l, rows, cols, total_count_l, allocator_callbacks);
    if (res != JMTX_RESULT_SUCCESS)
    {
        allocator_callbacks->free(allocator_callbacks->state, count_u);
        return res;
    }
    res = jmtxz_matrix_ccs_new(&u, rows, cols, total_count_u, allocator_callbacks);
    if (res != JMTX_RESULT_SUCCESS)
    {
        allocator_callbacks->free(allocator_callbacks->state, count_u);
        return res;
    }
    uint32_t* col_ends = u->end_of_column_offsets;
    //  Create end of line offsets for the upper matrix
    col_ends[0] = count_u[0];
    //  Compute cumsums for offsets
    for (uint32_t i = 1; i < cols; ++i)
    {
        col_ends[i] = count_u[i] + col_ends[i-1];
        count_u[i] = 0; //   Zero the column counts so that they can be reused later for counting bucket sizes
    }
    count_u[0] = 0;
    count_u[cols - 1] = 0;

    uint32_t lcnt = 0;
    for (uint32_t row = 0; row < rows; ++row)
    {
        uint32_t* in_cols;
        _Complex double* in_vals;
        uint32_t n_row = jmtxz_matrix_crs_get_row(a, row, &in_cols, &in_vals);

        for (uint32_t idx = 0; idx < n_row; ++idx)
        {
            const uint32_t col = in_cols[idx];
            const uint32_t ip = col > 0 ? col_ends[col-1] : 0;
            const uint32_t n_col = count_u[col];

            if (col >= row)
            {
                u->values[ip+n_col] = in_vals[idx];
                u->indices[ip+n_col] = row;
                count_u[col] += 1;
            }
            else
            {
                l->indices[lcnt] = col;
                l->values[lcnt] = in_vals[idx];
                lcnt += 1;
            }
        }
        l->indices[lcnt] = row;
        l->values[lcnt] = 1;
        lcnt += 1;
        l->end_of_row_offsets[row] = lcnt;
    }
    l->n_entries = l->end_of_row_offsets[rows - 1];
    u->n_entries = u->end_of_column_offsets[cols - 1];


    allocator_callbacks->free(allocator_callbacks->state, count_u);
    *p_l = l;
    *p_u = u;

    return JMTX_RESULT_SUCCESS;
}

jmtx_result jmtxz_decompose_ilu_crs(
        const jmtxz_matrix_crs* a, jmtxz_matrix_crs** p_l, jmtxz_matrix_ccs** p_u,
        const jmtx_allocator_callbacks* allocator_callbacks)
{
    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }

    //  L and U have at most this many entries (in the case that A is already triangular)
    const uint32_t n = a->base.rows;
    jmtxz_matrix_crs* l = NULL;
    jmtxz_matrix_ccs* u = NULL;
    jmtx_result res = jmtxz_decompose_ldu_split_crs(a, &l, &u, allocator_callbacks);
    if (res != JMTX_RESULT_SUCCESS)
    {
        return res;
    }

    for (uint32_t idx = 0; idx < n; ++idx)
    {
        //  Compute the row idx for l
        {
            uint32_t* l_idx;
            _Complex double* l_val;
            const uint32_t l_cnt = jmtxz_matrix_crs_get_row(l, idx, &l_idx, &l_val);
            for (uint32_t k = 0; k < l_cnt - 1; ++k)
            {
                uint32_t col = l_idx[k];
                uint32_t* u_idx;
                _Complex double* u_val;
                const uint32_t u_cnt = jmtxz_matrix_ccs_get_col(u, col, &u_idx, &u_val);
                const _Complex double dp = jmtxz_multiply_matrix_sparse_vectors(k, l_idx, l_val, u_cnt, u_idx, u_val);
                l_val[k] = (l_val[k] - dp) / u_val[u_cnt - 1];
                assert(u_idx[u_cnt - 1] == col);
            }
        }
        //  Compute the column idx for u
        {
            uint32_t* u_idx;
            _Complex double* u_val;
            const uint32_t u_cnt = jmtxz_matrix_ccs_get_col(u, idx, &u_idx, &u_val);
            for (uint32_t k = 0; k < u_cnt; ++k)
            {
                uint32_t row = u_idx[k];
                uint32_t* l_idx;
                _Complex double* l_val;
                const uint32_t l_cnt = jmtxz_matrix_crs_get_row(l, row, &l_idx, &l_val);
                const _Complex double dp = jmtxz_multiply_matrix_sparse_vectors(k, u_idx, u_val, l_cnt, l_idx, l_val);
                u_val[k] -= dp;
                assert(l_idx[l_cnt - 1] == row);
                assert(l_val[l_cnt - 1] == 1);
            }
        }
    }

    *p_u = u;
    *p_l = l;
    return JMTX_RESULT_SUCCESS;
}
