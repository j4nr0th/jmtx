#include "band_lu_decomposition.h"
#include "../matrices/band_row_major.h"
#include <assert.h>

jmtx_result JMTX_NAME_TYPED(decompose_lu_brm)(const JMTX_NAME_TYPED(matrix_brm) * a, JMTX_NAME_TYPED(matrix_brm) * *p_l,
                                              JMTX_NAME_TYPED(matrix_brm) * *p_u,
                                              const jmtx_allocator_callbacks *allocator_callbacks)
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
    const JMTX_FAST_INT_T lbw = a->lower_bandwidth;
    const JMTX_FAST_INT_T ubw = a->upper_bandwidth;
    const JMTX_INDEX_T max_entries = lbw + ubw + 1;
    const JMTX_INDEX_T n = a->base.rows;
    JMTX_NAME_TYPED(matrix_brm) *l = NULL;
    JMTX_NAME_TYPED(matrix_brm) *u = NULL;

    JMTX_SCALAR_T *const p_values =
        allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*p_values) * 2 * max_entries);
    if (p_values == NULL)
    {
        return JMTX_RESULT_BAD_ALLOC;
    }

    jmtx_result res = JMTX_NAME_TYPED(matrix_brm_new)(&l, n, n, 0, lbw, NULL, allocator_callbacks);
    if (res != JMTX_RESULT_SUCCESS)
    {
        //  Can't make matrix :(
        allocator_callbacks->free(allocator_callbacks->state, p_values);
        return res;
    }
    res = JMTX_NAME_TYPED(matrix_brm_new)(&u, n, n, ubw, 0, NULL, allocator_callbacks);
    if (res != JMTX_RESULT_SUCCESS)
    {
        JMTX_NAME_TYPED(matrix_brm_destroy)(l);
        allocator_callbacks->free(allocator_callbacks->state, p_values);
        //  Can't make matrix :(
        return res;
    }

    for (JMTX_FAST_INT_T i = 0; i < n; ++i)
    {
        JMTX_FAST_INT_T first_l = JMTX_NAME_TYPED(matrix_brm_first_pos_in_row)(l, i);
        JMTX_SCALAR_T *lwr_elements = NULL;
        JMTX_NAME_TYPED(matrix_brm_get_row)(l, i, &lwr_elements);
        JMTX_FAST_INT_T k = 0;
        JMTX_SCALAR_T *a_elements;
        (void)JMTX_NAME_TYPED(matrix_brm_get_row)(a, i, &a_elements);
        JMTX_SCALAR_T *const upr_elements = p_values + max_entries;
        for (JMTX_FAST_INT_T pl = first_l; pl < i; ++pl)
        {
            const JMTX_FAST_INT_T count_upr = JMTX_NAME_TYPED(matrix_brm_get_col)(u, pl, upr_elements);
            JMTX_SCALAR_T v = 0;
            JMTX_FAST_INT_T first_u = JMTX_NAME_TYPED(matrix_brm_first_pos_in_col)(u, pl);
            JMTX_FAST_INT_T begin, end;
            if (first_l < first_u)
            {
                begin = first_u;
            }
            else // if (first_l >= first_u)
            {
                begin = first_l;
            }
            if (pl < i)
            {
                end = pl;
            }
            else // if (pl >= i)
            {
                end = i;
            }

            for (JMTX_FAST_INT_T j = begin; j < end; ++j)
            {
                v += lwr_elements[j - first_l] * upr_elements[j - first_u];
            }

            lwr_elements[k] = (a_elements[k] - v) / upr_elements[count_upr - 1];
            k += 1;
        }
        lwr_elements[k] = 1.0f;

        const JMTX_FAST_INT_T first_u = JMTX_NAME_TYPED(matrix_brm_first_pos_in_col)(u, i);
        k = 0;
        (void)JMTX_NAME_TYPED(matrix_brm_get_col)(a, i, p_values);
        a_elements = p_values;
        for (JMTX_FAST_INT_T pu = first_u; pu <= i; ++pu)
        {
            JMTX_NAME_TYPED(matrix_brm_get_row)(l, pu, &lwr_elements);
            JMTX_SCALAR_T v = 0;
            first_l = JMTX_NAME_TYPED(matrix_brm_first_pos_in_row)(l, pu);

            JMTX_FAST_INT_T begin, end;
            if (first_l < first_u)
            {
                begin = first_u;
            }
            else // if (first_l >= first_u)
            {
                begin = first_l;
            }
            if (pu < i)
            {
                end = pu;
            }
            else // if (pu >= i)
            {
                end = i;
            }
            for (JMTX_FAST_INT_T j = begin; j < end; ++j)
            {
                v += lwr_elements[j - first_l] * upr_elements[j - first_u];
            }
            upr_elements[k] = (a_elements[k] - v);
            k += 1;
        }
        JMTX_NAME_TYPED(matrix_brm_set_col)(u, i, upr_elements);
    }

    allocator_callbacks->free(allocator_callbacks->state, p_values);

    *p_u = u;
    *p_l = l;
    return JMTX_RESULT_SUCCESS;
}
