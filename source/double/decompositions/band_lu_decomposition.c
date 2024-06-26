// Automatically generated from source/float/solvers/band_lu_decomposition.c on Fri Dec  1 06:43:01 2023
//
// Created by jan on 24.11.2023.
//

#include <assert.h>
#include "../../../include/jmtx/double/decompositions/band_lu_decomposition.h"
#include "../matrices/band_row_major_internal.h"


jmtx_result jmtxd_decompose_lu_brm(
        const jmtxd_matrix_brm* a, jmtxd_matrix_brm** p_l, jmtxd_matrix_brm** p_u,
        const jmtx_allocator_callbacks* allocator_callbacks)
{
    if (!a)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (a->base.type != JMTXD_TYPE_BRM)
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

    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }
    else if (allocator_callbacks->alloc == NULL || allocator_callbacks->free == NULL)
    {
        return JMTX_RESULT_NULL_PARAM;
    }

    //  L and U have the exact same bandwidth as A
    const uint_fast32_t lbw = a->lower_bandwidth;
    const uint_fast32_t ubw = a->upper_bandwidth;
    const uint32_t max_entries = lbw + ubw + 1;
    const uint32_t n = a->base.rows;
    jmtxd_matrix_brm* l = NULL;
    jmtxd_matrix_brm* u = NULL;

    double* const p_values = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*p_values) * 2 * max_entries);
    if (p_values == NULL)
    {
        return JMTX_RESULT_BAD_ALLOC;
    }

    jmtx_result res = jmtxd_matrix_brm_new(&l, n, n, 0, lbw, NULL, allocator_callbacks);
    if (res != JMTX_RESULT_SUCCESS)
    {
        //  Can't make matrix :(
        allocator_callbacks->free(allocator_callbacks->state, p_values);
        return res;
    }
    res = jmtxd_matrix_brm_new(&u, n, n, ubw, 0, NULL, allocator_callbacks);
    if (res != JMTX_RESULT_SUCCESS)
    {
        jmtxd_matrix_brm_destroy(l);
        allocator_callbacks->free(allocator_callbacks->state, p_values);
        //  Can't make matrix :(
        return res;
    }

    for (uint_fast32_t i = 0; i < n; ++i)
    {
        uint_fast32_t first_l = jmtxd_matrix_brm_first_pos_in_row(l, i);
        double* lwr_elements = NULL;
        jmtxd_matrix_brm_get_row(l, i, &lwr_elements);
        uint_fast32_t k = 0;
        double* a_elements;
        (void)jmtxd_matrix_brm_get_row(a, i, &a_elements);
        double* const upr_elements = p_values + max_entries;
        for (uint_fast32_t pl = first_l; pl < i; ++pl)
        {
            const uint_fast32_t count_upr = jmtxd_matrix_brm_get_col(u, pl, upr_elements);
            double v = 0;
            uint_fast32_t first_u = jmtxd_matrix_brm_first_pos_in_col(u, pl);
            uint_fast32_t begin, end;
            if (first_l < first_u)
            {
                begin = first_u;
            }
            else //if (first_l >= first_u)
            {
                begin = first_l;
            }
            if (pl < i)
            {
                end = pl;
            }
            else //if (pl >= i)
            {
                end = i;
            }

            for (uint_fast32_t j = begin; j < end; ++j)
            {
                v += lwr_elements[j - first_l] * upr_elements[j - first_u];
            }

            lwr_elements[k] = (a_elements[k] - v) / upr_elements[count_upr - 1];
            k += 1;
        }
        lwr_elements[k++] = 1.0f;

        const uint_fast32_t first_u = jmtxd_matrix_brm_first_pos_in_col(u, i);
        k = 0;
        (void)jmtxd_matrix_brm_get_col(a, i, p_values);
        a_elements = p_values;
        for (uint_fast32_t pu = first_u; pu <= i; ++pu)
        {
            jmtxd_matrix_brm_get_row(l, pu, &lwr_elements);
            double v = 0;
            first_l = jmtxd_matrix_brm_first_pos_in_row(l, pu);

            uint_fast32_t begin, end;
            if (first_l < first_u)
            {
                begin = first_u;
            }
            else //if (first_l >= first_u)
            {
                begin = first_l;
            }
            if (pu < i)
            {
                end = pu;
            }
            else //if (pu >= i)
            {
                end = i;
            }
            for (uint_fast32_t j = begin; j < end; ++j)
            {
                v += lwr_elements[j - first_l] * upr_elements[j - first_u];
            }
            upr_elements[k] = (a_elements[k] - v);
            k += 1;
        }
        jmtxd_matrix_brm_set_col(u, i, upr_elements);
    }

    allocator_callbacks->free(allocator_callbacks->state, p_values);

    *p_u = u;
    *p_l = l;
    return JMTX_RESULT_SUCCESS;
}
