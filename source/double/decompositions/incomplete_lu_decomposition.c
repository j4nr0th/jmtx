// Automatically generated from source/float/solvers/incomplete_lu_decomposition.c on Sun Dec 17 18:02:22 2023
//
// Created by jan on 2.11.2023.
//

#include <assert.h>
#include <math.h>
#include "../../../include/jmtx/double/decompositions/incomplete_lu_decomposition.h"
#include "../matrices/sparse_row_compressed_internal.h"
#include "../matrices/sparse_column_compressed_internal.h"
#include "../../../tests/double/test_common.h"

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
jmtx_result jmtxds_decompose_ilu_cds(
        const jmtxd_matrix_crs* a, jmtxd_matrix_crs** p_l, jmtxd_matrix_ccs** p_u,
        const jmtx_allocator_callbacks* allocator_callbacks)
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
    return jmtxd_decompose_ilu_cds(a, p_l, p_u, allocator_callbacks);
}

jmtx_result jmtxd_decompose_ilu_cds(
        const jmtxd_matrix_crs* a, jmtxd_matrix_crs** p_l, jmtxd_matrix_ccs** p_u,
        const jmtx_allocator_callbacks* allocator_callbacks)
{
    if (!allocator_callbacks)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }

    //  L and U have at most this many entries (in the case that A is already triangular)
    const uint32_t max_entries = a->n_entries;
    const uint32_t n = a->base.rows;
    jmtx_result res;
    uint32_t max_elements_in_direction = 0;
    jmtxd_matrix_crs* l = NULL;
    jmtxd_matrix_ccs* u = NULL;
    for (uint32_t i = 0; i < n; ++i)
    {
        uint32_t n_dim = jmtxd_matrix_crs_entries_in_col(a, i);
        if (n_dim > max_elements_in_direction)
        {
            max_elements_in_direction = n_dim;
        }
        uint32_t* unused_idx;
        double* unused_val;
        n_dim = jmtxd_matrix_crs_get_row(a, i, &unused_idx, &unused_val);
        n_dim += 1;
        if (n_dim > max_elements_in_direction)
        {
            max_elements_in_direction = n_dim;
        }
    }

    uint32_t* p_indices = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*p_indices) * 2 * max_elements_in_direction);
    if (!p_indices)
    {
        return JMTX_RESULT_BAD_ALLOC;
    }
    double* p_values = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*p_values) * 2 * max_elements_in_direction);
    if (!p_values)
    {
        allocator_callbacks->free(allocator_callbacks->state, p_indices);
        return JMTX_RESULT_BAD_ALLOC;
    }


    res = jmtxd_matrix_crs_new(&l, n, n, max_entries, allocator_callbacks);
    if (res != JMTX_RESULT_SUCCESS)
    {
        //  Can't make matrix :(
        allocator_callbacks->free(allocator_callbacks->state, p_values);
        allocator_callbacks->free(allocator_callbacks->state, p_indices);
        return res;
    }
    res = jmtxd_matrix_ccs_new(&u, n, n, max_entries, allocator_callbacks);
    if (res != JMTX_RESULT_SUCCESS)
    {
        jmtxd_matrix_crs_destroy(l);
        allocator_callbacks->free(allocator_callbacks->state, p_values);
        allocator_callbacks->free(allocator_callbacks->state, p_indices);
        //  Can't make matrix :(
        return res;
    }

    for (uint32_t i = 0; i < n; ++i)
    {
        uint32_t c;
        uint32_t* indices;
        double* values;
        //  Get a row from A
        c = jmtxd_matrix_crs_get_row(a, i, &indices, &values);
        uint32_t r;
        for (r = 0; r < c && indices[r] < i; ++r)
        {
            const uint32_t m = indices[r];
            double v = 0;
            double va = values[r];

            uint32_t u_col_count;
            uint32_t* u_row_indices;
            double* u_val;
            //  Compute the product of row p of matrix L and column m of matrix U to update the entry L_pm
            u_col_count = jmtxd_matrix_ccs_get_col(u, m, &u_row_indices, &u_val);
            for (uint32_t k_l = 0, k_u = 0; k_l < r && k_u < u_col_count && p_indices[k_l] < m &&
                                            u_row_indices[k_u] < m;)
            {
                if (p_indices[k_l] == u_row_indices[k_u])
                {
                    v += p_values[k_l] * u_val[k_u];
                    k_l += 1;
                    k_u += 1;
                }
                else if (p_indices[k_l] < u_row_indices[k_u])
                {
                    k_l += 1;
                }
                else
                {
                    // (l_col_indices[k_l] < u_row_indices[k_u])
                    k_u += 1;
                }
            }
            assert(u_row_indices[u_col_count - 1] == m);
            v = (va - v) / u_val[u_col_count - 1];
            p_values[r] = v;
            p_indices[r] = m;
        }
        p_values[r] = 1.0f;
        p_indices[r] = i;

        //  Values bellow the diagonal go to L
        res = jmtxd_matrix_crs_build_row(l, i, r + 1, p_indices, p_values);
        if (res != JMTX_RESULT_SUCCESS)
        {
            allocator_callbacks->free(allocator_callbacks->state, p_values);
            allocator_callbacks->free(allocator_callbacks->state, p_indices);
            jmtxd_matrix_ccs_destroy(u);
            jmtxd_matrix_crs_destroy(l);
            return res;
        }

        //  Get column values from the matrix A
        c = jmtxd_matrix_crs_get_col(a, i, max_elements_in_direction, p_values, p_indices);
        double* const p_vu = p_values + max_elements_in_direction;
        uint32_t* const p_iu = p_indices + max_elements_in_direction;
        //  Compute the product of row m of matrix L and column p of matrix U to update the entry U_mp
//        const uint32_t r_before = r;
        for (r = 0; r < c && p_indices[r] <= i; ++r)
        {
            const uint32_t m = p_indices[r];
            double v = 0;
            double va = p_values[r];
            assert(va != 0.0f);
            uint32_t* l_col_indices;
            double* l_val;
            uint32_t l_row_count = jmtxd_matrix_crs_get_row(l, m, &l_col_indices, &l_val);
            for (uint32_t k_l = 0, k_u = 0; k_l < l_row_count && k_u < r && l_col_indices[k_l] < m && p_indices[k_u] < m;)
            {
                if (l_col_indices[k_l] == p_indices[k_u])
                {
                    v += l_val[k_l] * p_vu[k_u];
                    k_l += 1;
                    k_u += 1;
                }
                else if (l_col_indices[k_l] < p_indices[k_u])
                {
                    k_l += 1;
                }
                else
                {
                    // (l_col_indices[k_l] < u_row_indices[k_u])
                    k_u += 1;
                }
            }
            v = (va - v);
//                jmtxd_matrix_ccs_set_entry(u, m, p, v);
            p_vu[r] = v;
            p_iu[r] = m;
        }

        assert(p_iu[r - 1] == i);
        //  Values above the diagonal go to U
        res = jmtxd_matrix_ccs_build_col(u, i, r, p_iu, p_vu);
        if (res != JMTX_RESULT_SUCCESS)
        {
            allocator_callbacks->free(allocator_callbacks->state, p_values);
            allocator_callbacks->free(allocator_callbacks->state, p_indices);
            jmtxd_matrix_ccs_destroy(u);
            jmtxd_matrix_crs_destroy(l);
            return res;
        }
    }

    allocator_callbacks->free(allocator_callbacks->state, p_values);
    allocator_callbacks->free(allocator_callbacks->state, p_indices);
    p_values = NULL;
    p_indices = NULL;

    *p_u = u;
    *p_l = l;
    return JMTX_RESULT_SUCCESS;
}
