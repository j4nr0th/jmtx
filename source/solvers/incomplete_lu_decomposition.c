//
// Created by jan on 2.11.2023.
//

#include "incomplete_lu_decomposition.h"
#include "../matrices/sparse_row_compressed_internal.h"
#include "../matrices/sparse_column_compressed.h"

jmtx_result jmtx_incomplete_lu_crs(
        jmtx_matrix_crs* a, jmtx_matrix_crs** p_l, jmtx_matrix_crs** p_u, float stagnation, uint32_t max_iterations,
        const jmtx_allocator_callbacks* allocator_callbacks)
{
    (void) stagnation;
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
    jmtx_result res = JMTX_RESULT_SUCCESS;
    uint32_t max_elements_in_column = 0;
    jmtx_matrix_crs* l = NULL;
    jmtx_matrix_ccs* u = NULL;
    for (uint32_t i = 0; i < n; ++i)
    {
        uint32_t n_col;
        res = jmtx_matrix_crs_entries_in_col(a, i, &n_col);
        if (res != JMTX_RESULT_SUCCESS)
        {
            return res;
        }
        if (n_col > max_elements_in_column)
        {
            max_elements_in_column = n_col;
        }
    }

    uint32_t* const column_indices = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*column_indices) * max_elements_in_column);
    if (!column_indices)
    {
        return JMTX_RESULT_BAD_ALLOC;
    }
    float* const column_values = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*column_values) * max_elements_in_column);
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
        if (pos != i)
        {
            //  There's no diagonal, meaning that LU decomposition can't be done (without pivoting, but fuck that)
            res = JMTX_RESULT_BAD_MATRIX;
            allocator_callbacks->free(allocator_callbacks->state, column_values);
            allocator_callbacks->free(allocator_callbacks->state, column_indices);
            goto failed;
        }
        //  Values bellow the diagonal go to L
        res = jmtx_matrix_crs_build_row(l, i, pos, indices, values);
        if (res != JMTX_RESULT_SUCCESS)
        {
            allocator_callbacks->free(allocator_callbacks->state, column_values);
            allocator_callbacks->free(allocator_callbacks->state, column_indices);
            goto failed;
        }

        //  Get column values from the matrix A
        res = jmtx_matrix_crs_get_col(a, i, max_elements_in_column, &c, column_values, column_indices);
        if (res != JMTX_RESULT_SUCCESS)
        {
            allocator_callbacks->free(allocator_callbacks->state, column_values);
            allocator_callbacks->free(allocator_callbacks->state, column_indices);
            goto failed;
        }

        //  Put the column in the matrix U
        pos = jmtx_internal_find_last_leq_value(c, column_indices, i);
        if (pos != i)
        {
            //  There's no diagonal, meaning that LU decomposition can't be done (without pivoting, but fuck that)
            res = JMTX_RESULT_BAD_MATRIX;
            allocator_callbacks->free(allocator_callbacks->state, column_values);
            allocator_callbacks->free(allocator_callbacks->state, column_indices);
            goto failed;
        }
        //  Values above the diagonal go to U
        res = jmtx_matrix_ccs_build_col(u, i, pos, indices, values);
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
    while (iteration_count < max_iterations)
    {
        //  Loop over every row of A and update it
        for (uint32_t p = 0; p < n; ++p)
        {
            // Column access is slower than row access, so it is best to fetch a column of U once
        }

        //  Check convergence
        iteration_count += 1;
    }


    return JMTX_RESULT_SUCCESS;


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
