#include "incomplete_lu_decomposition.h"
#include <assert.h>
#include <math.h>

#include "../matrices/sparse_multiplication.h"
#include "../matrices/sparse_row_compressed.h"
#include "../matrices/sparse_column_compressed.h"
#include "../matrices/sparse_diagonal_compressed.h"
#include "../matrices/sparse_row_compressed.h"

jmtx_result JMTX_NAME_TYPED(decompose_ldu_split_crs)(const JMTX_NAME_TYPED(matrix_crs) * a,
                                                     JMTX_NAME_TYPED(matrix_crs) * *p_l,
                                                     JMTX_NAME_TYPED(matrix_ccs) * *p_u,
                                                     const jmtx_allocator_callbacks *allocator_callbacks)
{
    if (!allocator_callbacks)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }

    JMTX_NAME_TYPED(matrix_crs) * l;
    JMTX_NAME_TYPED(matrix_ccs) * u;

    JMTX_INDEX_T cols = a->base.cols;
    JMTX_INDEX_T *count_u = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*count_u) * cols);
    if (!count_u)
    {
        return JMTX_RESULT_BAD_ALLOC;
    }
    JMTX_INDEX_T rows = a->base.rows;
    JMTX_INDEX_T total_count_l = 0, total_count_u = 0;
    memset(count_u, 0, sizeof(*count_u) * cols);

    for (JMTX_INDEX_T row = 0; row < rows; ++row)
    {
        JMTX_INDEX_T *pcols;
        JMTX_SCALAR_T *pvals;
        JMTX_INDEX_T ncols = JMTX_NAME_TYPED(matrix_crs_get_row)(a, row, &pcols, &pvals);
        for (JMTX_INDEX_T n = 0; n < ncols; ++n)
        {
            JMTX_INDEX_T col = pcols[n];
            count_u[col] += (col >= row);
            total_count_u += (col >= row);
            total_count_l += (col <= row);
        }
    }

    jmtx_result res = JMTX_NAME_TYPED(matrix_crs_new)(&l, rows, cols, total_count_l, allocator_callbacks);
    if (res != JMTX_RESULT_SUCCESS)
    {
        allocator_callbacks->free(allocator_callbacks->state, count_u);
        return res;
    }
    res = JMTX_NAME_TYPED(matrix_ccs_new)(&u, rows, cols, total_count_u, allocator_callbacks);
    if (res != JMTX_RESULT_SUCCESS)
    {
        allocator_callbacks->free(allocator_callbacks->state, count_u);
        return res;
    }
    JMTX_INDEX_T *col_ends = u->end_of_column_offsets;
    //  Create end of line offsets for the upper matrix
    col_ends[0] = count_u[0];
    //  Compute cumsums for offsets
    for (JMTX_INDEX_T i = 1; i < cols; ++i)
    {
        col_ends[i] = count_u[i] + col_ends[i - 1];
        count_u[i] = 0; //   Zero the column counts so that they can be reused later for counting bucket sizes
    }
    count_u[0] = 0;
    count_u[cols - 1] = 0;

    JMTX_INDEX_T lcnt = 0;
    for (JMTX_INDEX_T row = 0; row < rows; ++row)
    {
        JMTX_INDEX_T *in_cols;
        JMTX_SCALAR_T *in_vals;
        JMTX_INDEX_T n_row = JMTX_NAME_TYPED(matrix_crs_get_row)(a, row, &in_cols, &in_vals);

        for (JMTX_INDEX_T idx = 0; idx < n_row; ++idx)
        {
            const JMTX_INDEX_T col = in_cols[idx];
            const JMTX_INDEX_T ip = col > 0 ? col_ends[col - 1] : 0;
            const JMTX_INDEX_T n_col = count_u[col];

            if (col >= row)
            {
                u->values[ip + n_col] = in_vals[idx];
                u->indices[ip + n_col] = row;
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

jmtx_result JMTX_NAME_TYPED(decompose_ilu_crs)(const JMTX_NAME_TYPED(matrix_crs) * a,
                                               JMTX_NAME_TYPED(matrix_crs) * *p_l, JMTX_NAME_TYPED(matrix_ccs) * *p_u,
                                               const jmtx_allocator_callbacks *allocator_callbacks)
{
    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }

    //  L and U have at most this many entries (in the case that A is already triangular)
    const JMTX_INDEX_T n = a->base.rows;
    JMTX_NAME_TYPED(matrix_crs) *l = NULL;
    JMTX_NAME_TYPED(matrix_ccs) *u = NULL;
    jmtx_result res = JMTX_NAME_TYPED(decompose_ldu_split_crs)(a, &l, &u, allocator_callbacks);
    if (res != JMTX_RESULT_SUCCESS)
    {
        return res;
    }

    for (JMTX_INDEX_T idx = 0; idx < n; ++idx)
    {
        //  Compute the row idx for l
        {
            JMTX_INDEX_T *l_idx;
            JMTX_SCALAR_T *l_val;
            const JMTX_INDEX_T l_cnt = JMTX_NAME_TYPED(matrix_crs_get_row)(l, idx, &l_idx, &l_val);
            for (JMTX_INDEX_T k = 0; k < l_cnt - 1; ++k)
            {
                JMTX_INDEX_T col = l_idx[k];
                JMTX_INDEX_T *u_idx;
                JMTX_SCALAR_T *u_val;
                const JMTX_INDEX_T u_cnt = JMTX_NAME_TYPED(matrix_ccs_get_col)(u, col, &u_idx, &u_val);
                const JMTX_SCALAR_T dp =
                    JMTX_NAME_TYPED(multiply_matrix_sparse_vectors)(k, l_idx, l_val, u_cnt, u_idx, u_val);
                l_val[k] = (l_val[k] - dp) / u_val[u_cnt - 1];
                assert(u_idx[u_cnt - 1] == col);
            }
        }
        //  Compute the column idx for u
        {
            JMTX_INDEX_T *u_idx;
            JMTX_SCALAR_T *u_val;
            const JMTX_INDEX_T u_cnt = JMTX_NAME_TYPED(matrix_ccs_get_col)(u, idx, &u_idx, &u_val);
            for (JMTX_INDEX_T k = 0; k < u_cnt; ++k)
            {
                JMTX_INDEX_T row = u_idx[k];
                JMTX_INDEX_T *l_idx;
                JMTX_SCALAR_T *l_val;
                const JMTX_INDEX_T l_cnt = JMTX_NAME_TYPED(matrix_crs_get_row)(l, row, &l_idx, &l_val);
                const JMTX_SCALAR_T dp =
                    JMTX_NAME_TYPED(multiply_matrix_sparse_vectors)(k, u_idx, u_val, l_cnt, l_idx, l_val);
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
