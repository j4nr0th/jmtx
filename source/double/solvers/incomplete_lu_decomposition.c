// Automatically generated from source/float/solvers/incomplete_lu_decomposition.c on Thu Nov 30 20:33:08 2023
//
// Created by jan on 2.11.2023.
//

#include <assert.h>
#include <math.h>
#include "../../../include/jmtx/double/solvers/incomplete_lu_decomposition.h"
#include "../matrices/sparse_row_compressed_internal.h"
#include "../matrices/sparse_column_compressed_internal.h"
#include "../../../tests/double/test_common.h"

jmtx_result jmtxd_incomplete_lu_crs(
        const jmtxd_matrix_crs* a, jmtxd_matrix_crs** p_l, jmtxd_matrix_ccs** p_u,
        const jmtx_allocator_callbacks* allocator_callbacks)
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

jmtx_result jmtxd_incomplete_lu_crs_parallel(
        jmtxd_matrix_crs* a, jmtxd_matrix_crs** p_l, jmtxd_matrix_ccs** p_u,
        const jmtx_allocator_callbacks* allocator_callbacks, jmtxd_solver_arguments* args)
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
        n_dim = 0;
        n_dim = jmtxd_matrix_crs_get_row(a, i, &unused_idx, &unused_val);

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
    double* column_values = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*column_values) * max_elements_in_direction);
    if (!column_values)
    {
        allocator_callbacks->free(allocator_callbacks->state, column_indices);
        return JMTX_RESULT_BAD_ALLOC;
    }


    res = jmtxd_matrix_crs_new(&l, n, n, max_entries, allocator_callbacks);
    if (res != JMTX_RESULT_SUCCESS)
    {
        //  Can't make matrix :(
        allocator_callbacks->free(allocator_callbacks->state, column_values);
        allocator_callbacks->free(allocator_callbacks->state, column_indices);
        return res;
    }
    res = jmtxd_matrix_ccs_new(&u, n, n, max_entries, allocator_callbacks);
    if (res != JMTX_RESULT_SUCCESS)
    {
        jmtxd_matrix_crs_destroy(l);
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
        double* values;
        //  Get a row from A
        c = jmtxd_matrix_crs_get_row(a, i, &indices, &values);

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
        res = jmtxd_matrix_crs_build_row(l, i, pos, column_indices, column_values);
        if (res != JMTX_RESULT_SUCCESS)
        {
            allocator_callbacks->free(allocator_callbacks->state, column_values);
            allocator_callbacks->free(allocator_callbacks->state, column_indices);
            goto failed;
        }

        //  Get column values from the matrix A
        c = jmtxd_matrix_crs_get_col(a, i, max_elements_in_direction, column_values, column_indices);

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
        res = jmtxd_matrix_ccs_build_col(u, i, pos + 1, indices, values);
        if (res != JMTX_RESULT_SUCCESS)
        {
            allocator_callbacks->free(allocator_callbacks->state, column_values);
            allocator_callbacks->free(allocator_callbacks->state, column_indices);
            goto failed;
        }
    }
    //  Begin to iteratively refine the matrix entries
    //  In this section it is fine to modify entries in the matrix without synchronizing access, since the entries
    //  always get updated in-place, meaning that there is no moving of memory, just overwriting the same spot.
    //  If the code is serial, the iterations would be like Gauss-Seidel, if it were completely parallel, with
    //  all calculations ending at the same time, it would act like point Jacobi.
    uint32_t iteration_count = 0;
    double max_relative_change = 0;
    int converged = 0;

    allocator_callbacks->free(allocator_callbacks->state, column_values);
    allocator_callbacks->free(allocator_callbacks->state, column_indices);
    column_values = NULL;
    column_indices = NULL;
    int any_failed_allocating = 0;
    const uint32_t max_iterations = args->in_max_iterations;
    const double convergence = args->in_convergence_criterion;
    double* const p_err = args->opt_error_evolution;
#pragma omp parallel default(none) shared(iteration_count, max_elements_in_direction, max_relative_change,\
max_iterations, allocator_callbacks, a, l, u, n, convergence, converged, any_failed_allocating, p_err)
    {
        uint32_t* const col_idx = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*col_idx) * max_elements_in_direction);
        double* const col_val = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*col_val) * max_elements_in_direction);
        if (!col_idx || !col_val)
        {
            any_failed_allocating = 1;
        }
#pragma omp barrier

        if (any_failed_allocating == 0)
        {
            while (iteration_count < max_iterations && converged == 0)
            {
#pragma omp master
                {
                    max_relative_change = 0;
                }
                double local_max_relative_change = 0;
                //  Loop over every row/column pair of A and update it
#pragma omp for schedule(static)
                for (uint32_t p = 0; p < n; ++p)
                {
                    //  Update the p-th row of L
                    uint32_t n_items;
                    uint32_t* positions;
                    double* values;
                    n_items = jmtxd_matrix_crs_get_row(a, p, &positions, &values);
                    uint32_t out_items;
                    uint32_t* out_positions;
                    double* out_values;
                    out_items = jmtxd_matrix_crs_get_row(l, p, &out_positions, &out_values);
                    //  Loop over elements of L which can be updated
                    for (uint32_t r = 0, k = 0; r < n_items && (int) positions[r] < (int) p; ++r, ++k)
                    {
                        const uint32_t m = positions[r];
                        double v = 0;
                        double va = values[r];

                        uint32_t u_col_count;
                        uint32_t* u_row_indices;
                        double* u_val;
                        //  Compute the product of row p of matrix L and column m of matrix U to update the entry L_pm
                        u_col_count = jmtxd_matrix_ccs_get_col(u, m, &u_row_indices, &u_val);
                        for (uint32_t k_l = 0, k_u = 0; k_l < out_items && k_u < u_col_count && out_positions[k_l] < m &&
                                                        u_row_indices[k_u] < m;)
                        {
                            if (out_positions[k_l] == u_row_indices[k_u])
                            {
                                v += out_values[k_l] * u_val[k_u];
                                k_l += 1;
                                k_u += 1;
                            }
                            else if (out_positions[k_l] < u_row_indices[k_u])
                            {
                                k_l += 1;
                            }
                            else
                            {
                                // (l_col_indices[k_l] < u_row_indices[k_u])
                                k_u += 1;
                            }
                        }
                        v = (va - v) / u_val[u_col_count - 1];
                        const double relative_change = fabs((v - out_values[k]) / out_values[k]);
                        if (relative_change > local_max_relative_change)
                        {
                            local_max_relative_change = relative_change;
                        }
                        out_values[k] = v;
                    }

                    //  Update the p-th column of U
                    n_items =  jmtxd_matrix_crs_get_col(a, p, max_elements_in_direction, col_val, col_idx);
                    out_items = jmtxd_matrix_ccs_get_col(u, p, &out_positions, &out_values);
                    //  Compute the product of row m of matrix L and column p of matrix U to update the entry U_mp
                    for (uint32_t r = 0, k = 0; r < n_items && col_idx[r] < p + 1; ++r, ++k)
                    {
                        const uint32_t m = col_idx[r];
                        double v = 0;
                        double va = col_val[r];
                        uint32_t l_row_count;
                        uint32_t* l_col_indices;
                        double* l_val;
                        l_row_count = jmtxd_matrix_crs_get_row(l, m, &l_col_indices, &l_val);
                        for (uint32_t k_l = 0, k_u = 0; k_l < l_row_count && k_u < out_items && l_col_indices[k_l] < m &&
                                                        out_positions[k_u] < m;)
                        {
                            if (l_col_indices[k_l] == out_positions[k_u])
                            {
                                v += l_val[k_l] * out_values[k_u];
                                k_l += 1;
                                k_u += 1;
                            }
                            else if (l_col_indices[k_l] < out_positions[k_u])
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
                        const double relative_change = fabs((v - out_values[k]) / out_values[k]);
                        if (relative_change > local_max_relative_change)
                        {
                            local_max_relative_change = relative_change;
                        }
                        out_values[k] = v;
                    }
                }
#pragma omp critical
                {
                    if (local_max_relative_change > max_relative_change)
                    {
                        max_relative_change = local_max_relative_change;
                    }
                }

                //  Check convergence and record convergence
#pragma omp single
                {
                    if (convergence > max_relative_change)
                    {
                        converged = 1;
                    }
                    if (p_err)
                    {
                        p_err[iteration_count] = max_relative_change;
                    }
                    iteration_count += 1;
                }
            }
        }

        allocator_callbacks->free(allocator_callbacks->state, col_idx);
        allocator_callbacks->free(allocator_callbacks->state, col_val);
    }

    if (any_failed_allocating)
    {
        jmtxd_matrix_crs_destroy(l);
        jmtxd_matrix_ccs_destroy(u);
        return JMTX_RESULT_BAD_ALLOC;
    }

    *p_u = u;
    *p_l = l;
    args->out_last_iteration = iteration_count;
    args->out_last_error = max_relative_change;
    return converged ? JMTX_RESULT_SUCCESS : JMTX_RESULT_NOT_CONVERGED;


    failed:
    if (u)
    {
        jmtxd_matrix_ccs_destroy(u);
    }
    if (l)
    {
        jmtxd_matrix_crs_destroy(l);
    }
    return res;
}

jmtx_result jmtxd_incomplete_lu_crs_old(
        const jmtxd_matrix_crs* a, jmtxd_matrix_crs** p_l, jmtxd_matrix_ccs** p_u, double convergence, uint32_t max_iterations,
        double* final_max_change, uint32_t* p_last_iteration, const jmtx_allocator_callbacks* allocator_callbacks)
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

    uint32_t* column_indices = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*column_indices) * max_elements_in_direction);
    if (!column_indices)
    {
        return JMTX_RESULT_BAD_ALLOC;
    }
    double* column_values = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*column_values) * max_elements_in_direction);
    if (!column_values)
    {
        allocator_callbacks->free(allocator_callbacks->state, column_indices);
        return JMTX_RESULT_BAD_ALLOC;
    }


    res = jmtxd_matrix_crs_new(&l, n, n, max_entries, allocator_callbacks);
    if (res != JMTX_RESULT_SUCCESS)
    {
        //  Can't make matrix :(
        allocator_callbacks->free(allocator_callbacks->state, column_values);
        allocator_callbacks->free(allocator_callbacks->state, column_indices);
        return res;
    }
    res = jmtxd_matrix_ccs_new(&u, n, n, max_entries, allocator_callbacks);
    if (res != JMTX_RESULT_SUCCESS)
    {
        jmtxd_matrix_crs_destroy(l);
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
        double* values;
        //  Get a row from A
        c = jmtxd_matrix_crs_get_row(a, i, &indices, &values);

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
        res = jmtxd_matrix_crs_build_row(l, i, pos, column_indices, column_values);
        if (res != JMTX_RESULT_SUCCESS)
        {
            allocator_callbacks->free(allocator_callbacks->state, column_values);
            allocator_callbacks->free(allocator_callbacks->state, column_indices);
            goto failed;
        }

        //  Get column values from the matrix A
        c = jmtxd_matrix_crs_get_col(a, i, max_elements_in_direction, column_values, column_indices);

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
        res = jmtxd_matrix_ccs_build_col(u, i, pos + 1, indices, values);
        if (res != JMTX_RESULT_SUCCESS)
        {
            allocator_callbacks->free(allocator_callbacks->state, column_values);
            allocator_callbacks->free(allocator_callbacks->state, column_indices);
            goto failed;
        }
    }
    //  Begin to iteratively refine the matrix entries
    //  In this section it is fine to modify entries in the matrix without synchronizing access, since the entries
    //  always get updated in-place, meaning that there is no moving of memory, just overwriting the same spot.
    //  If the code is serial, the iterations would be like Gauss-Seidel, if it were completely parallel, with
    //  all calculations ending at the same time, it would act like point Jacobi.
    uint32_t iteration_count = 0;
    double max_relative_change = 0;
    int converged = 0;
    while (iteration_count < max_iterations)
    {
        max_relative_change = 0;
        //  Loop over every row/column pair of A and update it
        for (uint32_t p = 0; p < n; ++p)
        {
            //  Update the p-th row of L
            uint32_t* positions;
            double* values;
            uint32_t n_items = jmtxd_matrix_crs_get_row(a, p, &positions, &values);
            uint32_t* out_positions;
            double* out_values;
            uint32_t out_items = jmtxd_matrix_crs_get_row(l, p, &out_positions, &out_values);
            //  Loop over elements of L which can be updated
            for (uint32_t r = 0, k = 0; r < n_items && (int)positions[r] < (int)p; ++r, ++k)
            {
                const uint32_t m = positions[r];
                assert(out_positions[k] == m);
                double v = 0;
                double va = values[r];
                assert(va != 0.0f);

                uint32_t  u_col_count;
                uint32_t* u_row_indices;
                double *u_val;
                //  Compute the product of row p of matrix L and column m of matrix U to update the entry L_pm
                u_col_count = jmtxd_matrix_ccs_get_col(u, m, &u_row_indices, &u_val);
                for (uint32_t k_l = 0, k_u = 0; k_l < out_items && k_u < u_col_count && out_positions[k_l] < m &&
                                                u_row_indices[k_u] < m;)
                {
                    if (out_positions[k_l] == u_row_indices[k_u])
                    {
                        v += out_values[k_l] * u_val[k_u];
                        k_l += 1;
                        k_u += 1;
                    }
                    else if (out_positions[k_l] < u_row_indices[k_u])
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
                const double relative_change = fabs((v - out_values[k]) / out_values[k]);
                if (relative_change > max_relative_change)
                {
                    max_relative_change = relative_change;
                }
#ifndef NDEBUG
                double v1 = jmtxd_matrix_crs_get_entry(l, p, m);
                assert(v1 != 0.0f);
#endif
//                jmtxd_matrix_crs_set_entry(l, p, m, v);
                out_values[k] = v;
            }

            //  Update the p-th column of U
            n_items = jmtxd_matrix_crs_get_col(a, p, max_elements_in_direction, column_values, column_indices);
            out_items = jmtxd_matrix_ccs_get_col(u, p, &out_positions, &out_values);
            //  Compute the product of row m of matrix L and column p of matrix U to update the entry U_mp
            for (uint32_t r = 0, k = 0; r < n_items && column_indices[r] < p + 1; ++r, ++k)
            {
                const uint32_t m = column_indices[r];
                assert(out_positions[k] == m);
                double v = 0;
                double va = column_values[r];
                assert(va != 0.0f);
                uint32_t* l_col_indices;
                double* l_val;
                assert(va != 0.0f);
                uint32_t l_row_count = jmtxd_matrix_crs_get_row(l, m, &l_col_indices, &l_val);
                for (uint32_t k_l = 0, k_u = 0; k_l < l_row_count && k_u < out_items && l_col_indices[k_l] < m &&
                                                out_positions[k_u] < m;)
                {
                    if (l_col_indices[k_l] == out_positions[k_u])
                    {
                        v += l_val[k_l] * out_values[k_u];
                        k_l += 1;
                        k_u += 1;
                    }
                    else if (l_col_indices[k_l] < out_positions[k_u])
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
                const double relative_change = fabs((v - out_values[k]) / out_values[k]);
                if (relative_change > max_relative_change)
                {
                    max_relative_change = relative_change;
                }
#ifndef NDEBUG
                double v1 = jmtxd_matrix_ccs_get_entry(u, m, p);
                assert(v1 != 0.0f);
#endif
//                jmtxd_matrix_ccs_set_entry(u, m, p, v);
                out_values[k] = v;
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
    allocator_callbacks->free(allocator_callbacks->state, column_values);
    allocator_callbacks->free(allocator_callbacks->state, column_indices);
    column_values = NULL;
    column_indices = NULL;

    *p_u = u;
    *p_l = l;
    *p_last_iteration = iteration_count;
    *final_max_change = max_relative_change;
    return converged ? JMTX_RESULT_SUCCESS : JMTX_RESULT_NOT_CONVERGED;


    failed:
    if (u)
    {
        jmtxd_matrix_ccs_destroy(u);
    }
    if (l)
    {
        jmtxd_matrix_crs_destroy(l);
    }
    return res;
}


