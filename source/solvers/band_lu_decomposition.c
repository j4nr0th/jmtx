//
// Created by jan on 24.11.2023.
//

#include <assert.h>
#include "band_lu_decomposition.h"
#include "../matrices/band_row_major_internal.h"

jmtx_result jmtx_band_lu_decomposition_brm(
        const jmtx_matrix_brm* a, jmtx_matrix_brm** p_l, jmtx_matrix_brm** p_u,
        const jmtx_allocator_callbacks* allocator_callbacks)
{
    if (!a)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (a->base.type != JMTX_TYPE_BRM)
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

    //  L and U have the exact same bandwidth as A
    const uint_fast32_t lbw = a->lower_bandwidth;
    const uint_fast32_t ubw = a->upper_bandwidth;
    const uint32_t max_entries = lbw + ubw + 1;
    const uint32_t n = a->base.rows;
    jmtx_result res;
    jmtx_matrix_brm* l = NULL;
    jmtx_matrix_brm* u = NULL;

    float* const p_values = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*p_values) * 2 * max_entries);
    if (!p_values)
    {
        return JMTX_RESULT_BAD_ALLOC;
    }
    float* const build_vals = p_values + max_entries;
    memset(p_values, 0xCC, sizeof(*p_values) * 2 * max_entries);    //TODO: remove

    res = jmtx_matrix_brm_new(&l, n, n, 0, lbw, NULL, allocator_callbacks);
    if (res != JMTX_RESULT_SUCCESS)
    {
        //  Can't make matrix :(
        allocator_callbacks->free(allocator_callbacks->state, p_values);
        return res;
    }
    res = jmtx_matrix_brm_new(&u, n, n, ubw, 0, NULL, allocator_callbacks);
    if (res != JMTX_RESULT_SUCCESS)
    {
        jmtx_matrix_brm_destroy(l);
        allocator_callbacks->free(allocator_callbacks->state, p_values);
        //  Can't make matrix :(
        return res;
    }

    for (uint_fast32_t i = 0; i < n; ++i)
    {
        uint_fast32_t c;
        float* values;
        //  Get a row from A
        c = jmtx_matrix_brm_get_row(a, i, &values);
        uint_fast32_t j;
        uint_fast32_t first = jmtx_matrix_brm_first_pos_in_row(a, i);
        float* lv = NULL;
        float* uv = p_values;
        uint_fast32_t lc = jmtx_matrix_brm_get_row(l, i, &lv);
        uint_fast32_t uc;
        for (j = first; j < i; ++j)
        {
            uc = jmtx_matrix_brm_get_col(u, j, uv);
            float v = 0;
            uint_fast32_t k;
            for (k = 0; k < j - first; ++k)
            {
                v += lv[k] * uv[k];
            }
            assert(uv[k] == jmtx_matrix_brm_get_entry(u, j, j));
            lv[j - first] = (values[j - first] - v) / uv[k];
        }
        lv[j - first] = 1.0f;
        memset(p_values, 0xCC, sizeof(*p_values) * 2 * max_entries);    //TODO: remove

        c = jmtx_matrix_brm_get_col(a, i, p_values);
        first = jmtx_matrix_brm_first_pos_in_col(u, i);
        uv = build_vals;
        for (j = first; j <= i; ++j)
        {
            float v = 0;
            lc = jmtx_matrix_brm_get_row(l, j, &lv);
            uint_fast32_t k;
            for (k = j; k < i; ++k)
            {
                v += lv[k - j] * uv[k - j];
            }
            uv[k - j] = p_values[k - j] - v;
        }
        jmtx_matrix_brm_set_col(u, i, build_vals);
        memset(p_values, 0xCC, sizeof(*p_values) * 2 * max_entries);    //TODO: remove
        (void)lc;
        (void)uc;
        (void)c;
    }

    allocator_callbacks->free(allocator_callbacks->state, p_values);

    *p_u = u;
    *p_l = l;
    return JMTX_RESULT_SUCCESS;
}
