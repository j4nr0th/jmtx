//
// Created by jan on 2.11.2023.
//

#include <assert.h>
#include <math.h>
#include "incomplete_lu_decomposition.h"
#include "../matrices/sparse_row_compressed_internal.h"
#include "../matrices/sparse_column_compressed_internal.h"
#include "../matrices/sparse_conversion.h"
#include "../../tests/test_common.h"

jmtx_result jmtx_incomplete_lu_crs(
        jmtx_matrix_crs* a, jmtx_matrix_crs** p_l, jmtx_matrix_ccs** p_u, float convergence, uint32_t max_iterations,
        float* final_max_change, uint32_t* p_last_iteration, const jmtx_allocator_callbacks* allocator_callbacks)
{
    (void) convergence;
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
    if (!p_l)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!p_u)
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
    const uint32_t max_entries = a->n_entries;
    const uint32_t n = a->base.rows;
    jmtx_result res;
    uint32_t max_elements_in_direction = 0;
    jmtx_matrix_crs* l = NULL;
    jmtx_matrix_ccs* u = NULL;
    for (uint32_t i = 0; i < n; ++i)
    {
        uint32_t n_dim;
        res = jmtx_matrix_crs_entries_in_col(a, i, &n_dim);
        if (res != JMTX_RESULT_SUCCESS)
        {
            return res;
        }
        if (n_dim > max_elements_in_direction)
        {
            max_elements_in_direction = n_dim;
        }
        uint32_t* unused_idx;
        float* unused_val;
        n_dim = 0;
        res = jmtx_matrix_crs_get_row(a, i, &n_dim, &unused_idx, &unused_val);
        if (res != JMTX_RESULT_SUCCESS)
        {
            return res;
        }
        n_dim += 1;
        if (n_dim > max_elements_in_direction)
        {
            max_elements_in_direction = n_dim;
        }
    }

    uint32_t* column_indices = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*column_indices) * max_elements_in_direction);
    if (!column_indices)
    {
        return JMTX_RESULT_BAD_ALLOC;
    }
    float* column_values = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*column_values) * max_elements_in_direction);
    if (!column_values)
    {
        allocator_callbacks->free(allocator_callbacks->state, column_indices);
        return JMTX_RESULT_BAD_ALLOC;
    }


    res = jmtx_matrix_crs_new(&l, n, n, max_entries, allocator_callbacks);
    if (res != JMTX_RESULT_SUCCESS)
    {
        //  Can't make matrix :(
        allocator_callbacks->free(allocator_callbacks->state, column_values);
        allocator_callbacks->free(allocator_callbacks->state, column_indices);
        return res;
    }
    res = jmtx_matrix_ccs_new(&u, n, n, max_entries, allocator_callbacks);
    if (res != JMTX_RESULT_SUCCESS)
    {
        jmtx_matrix_crs_destroy(l);
        allocator_callbacks->free(allocator_callbacks->state, column_values);
        allocator_callbacks->free(allocator_callbacks->state, column_indices);
        //  Can't make matrix :(
        return res;
    }

    //  Make the initial guess
    for (uint32_t i = 0; i < n; ++i)
    {
        uint32_t c;
        uint32_t* indices;
        float* values;
        //  Get a row from A
        res = jmtx_matrix_crs_get_row(a, i, &c, &indices, &values);
        if (res != JMTX_RESULT_SUCCESS)
        {
            allocator_callbacks->free(allocator_callbacks->state, column_values);
            allocator_callbacks->free(allocator_callbacks->state, column_indices);
            goto failed;
        }
        //  Find where the diagonal is located (if it is there)
        uint32_t pos = jmtx_internal_find_last_leq_value(c, indices, i);
        if (indices[pos] != i)
        {
            //  There's no diagonal, meaning that LU decomposition can't be done (without pivoting, but fuck that)
            res = JMTX_RESULT_BAD_MATRIX;
            allocator_callbacks->free(allocator_callbacks->state, column_values);
            allocator_callbacks->free(allocator_callbacks->state, column_indices);
            goto failed;
        }
        memcpy(column_values, values, sizeof(*values) * c);
        memcpy(column_indices, indices, sizeof(*indices) * c);
        column_values[pos] = 1.0f;
        column_indices[pos] = i;
        pos += 1;
        //  Values bellow the diagonal go to L
        res = jmtx_matrix_crs_build_row(l, i, pos, column_indices, column_values);
        if (res != JMTX_RESULT_SUCCESS)
        {
            allocator_callbacks->free(allocator_callbacks->state, column_values);
            allocator_callbacks->free(allocator_callbacks->state, column_indices);
            goto failed;
        }

        //  Get column values from the matrix A
        res = jmtx_matrix_crs_get_col(a, i, max_elements_in_direction, &c, column_values, column_indices);
        if (res != JMTX_RESULT_SUCCESS)
        {
            allocator_callbacks->free(allocator_callbacks->state, column_values);
            allocator_callbacks->free(allocator_callbacks->state, column_indices);
            goto failed;
        }

        //  Put the column in the matrix U
        pos = jmtx_internal_find_last_leq_value(c, column_indices, i);
        if (indices[pos] != i)
        {
            //  There's no diagonal, meaning that LU decomposition can't be done (without pivoting, but fuck that)
            res = JMTX_RESULT_BAD_MATRIX;
            allocator_callbacks->free(allocator_callbacks->state, column_values);
            allocator_callbacks->free(allocator_callbacks->state, column_indices);
            goto failed;
        }
        //  Values above the diagonal go to U
        res = jmtx_matrix_ccs_build_col(u, i, pos + 1, indices, values);
        if (res != JMTX_RESULT_SUCCESS)
        {
            allocator_callbacks->free(allocator_callbacks->state, column_values);
            allocator_callbacks->free(allocator_callbacks->state, column_indices);
            goto failed;
        }
    }
    allocator_callbacks->free(allocator_callbacks->state, column_values);
    allocator_callbacks->free(allocator_callbacks->state, column_indices);
    column_values = NULL;
    column_indices = NULL;
    //  Begin to iteratively refine the matrix entries
    //  In this section it is fine to modify entries in the matrix without synchronizing access, since the entries 
    //  always get updated in-place, meaning that there is no moving of memory, just overwriting the same spot. 
    //  If the code is serial, the iterations would be like Gauss-Seidel, if it were completely parallel, with 
    //  all calculations ending at the same time, it would act like point Jacobi.
    uint32_t iteration_count = 0;
    float max_relative_change = 0;
    int converged = 0;
    while (iteration_count < max_iterations)
    {
        max_relative_change = 0;
        //  Loop over every row of A and update it
        for (int p = 0; p < (int)n; ++p)
        {
            //  Update the p-th row of L

            for (int m = 0; m < p; ++m)
            {
                float v = 0;
                float va = 0;
                jmtx_matrix_crs_get_entry(a, p, m, &va);
                if (va == 0.0f)
                {
                    continue;
                }
                uint32_t l_row_count, u_col_count;
                uint32_t* l_col_indices, * u_row_indices;
                float* l_val, *u_val;
                jmtx_matrix_crs_get_row(l, p, &l_row_count, &l_col_indices, &l_val);
                jmtx_matrix_ccs_get_col(u, m, &u_col_count, &u_row_indices, &u_val);
                for (uint32_t k_l = 0, k_u = 0; k_l < l_row_count && k_u < u_col_count && (int)l_col_indices[k_l] < m &&
                        (int)u_row_indices[k_u] < m;)
                {
                    if (l_col_indices[k_l] == u_row_indices[k_u])
                    {
                        v += l_val[k_l] * u_val[k_u];
                        k_l += 1;
                        k_u += 1;
                    }
                    else if (l_col_indices[k_l] < u_row_indices[k_u])
                    {
                        k_l += 1;
                    }
                    else
                    {
                        // (l_col_indices[k_l] < u_row_indices[k_u])
                        k_u += 1;
                    }
                }
                assert((int)u_row_indices[u_col_count - 1] == m);
                v = (v - va) / u_val[u_col_count - 1];
                const float relative_change = fabsf((v - va) / va);
                if (relative_change > max_relative_change)
                {
                    max_relative_change = relative_change;
                }
                jmtx_matrix_crs_set_entry(l, p, m, v);
            }

            for (int m = 0; m < p + 1; ++m)
            {
                float v = 0;
                float va = 0;
                uint32_t l_row_count, u_col_count;
                uint32_t* l_col_indices, * u_row_indices;
                float* l_val, *u_val;
                jmtx_matrix_crs_get_entry(a, m, p, &va);
                if (va == 0.0f)
                {
                    continue;
                }
                jmtx_matrix_crs_get_row(l, m, &l_row_count, &l_col_indices, &l_val);
                jmtx_matrix_ccs_get_col(u, p, &u_col_count, &u_row_indices, &u_val);
                for (uint32_t k_l = 0, k_u = 0; k_l < l_row_count && k_u < u_col_count && (int)l_col_indices[k_l] < m &&
                        (int)u_row_indices[k_u] < m;)
                {
                    if (l_col_indices[k_l] == u_row_indices[k_u])
                    {
                        v += l_val[k_l] * u_val[k_u];
                        k_l += 1;
                        k_u += 1;
                    }
                    else if (l_col_indices[k_l] < u_row_indices[k_u])
                    {
                        k_l += 1;
                    }
                    else
                    {
                        // (l_col_indices[k_l] < u_row_indices[k_u])
                        k_u += 1;
                    }
                }
                v = (v - va);
                const float relative_change = fabsf((v - va) / va);
                if (relative_change > max_relative_change)
                {
                    max_relative_change = relative_change;
                }
                jmtx_matrix_ccs_set_entry(u, m, p, v);
            }
        }

        //  Check convergence
        if (convergence > max_relative_change)
        {
            converged = 1;
            break;
        }
        iteration_count += 1;
    }

    *p_u = u;
    *p_l = l;
    *p_last_iteration = iteration_count;
    *final_max_change = max_relative_change;
    return converged ? JMTX_RESULT_SUCCESS : JMTX_RESULT_NOT_CONVERGED;


failed:
    if (u)
    {
        jmtx_matrix_ccs_destroy(u);
    }
    if (l)
    {
        jmtx_matrix_crs_destroy(l);
    }
    return res;
}
