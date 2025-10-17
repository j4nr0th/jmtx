#include "band_row_major.h"
#include <assert.h>
#include <math.h>

static inline JMTX_FAST_INT_T brm_row_offset(const JMTX_NAME_TYPED(matrix_brm) mtx[JMTX_ARRAY_ATTRIB(const static 1)],
                                             const JMTX_FAST_INT_T row)
{
    JMTX_FAST_INT_T offset = (mtx->lower_bandwidth + 1 + mtx->upper_bandwidth) * row;
    if (row < mtx->lower_bandwidth)
    {
        offset -= row * mtx->lower_bandwidth - (row - 1) * row / 2;
    }
    else
    {
        offset -= (mtx->lower_bandwidth + 1) * mtx->lower_bandwidth / 2;
    }
    if (row > mtx->base.rows - 1 - mtx->upper_bandwidth)
    {
        const JMTX_FAST_INT_T offset_row = row - (mtx->base.rows - 1 - mtx->upper_bandwidth); //  rows offset from
        offset -= offset_row * (offset_row - 1) / 2;
    }
    return offset;
}

static inline JMTX_FAST_INT_T brm_row_len(const JMTX_NAME_TYPED(matrix_brm) mtx[JMTX_ARRAY_ATTRIB(const static 1)],
                                          const JMTX_FAST_INT_T row)
{
    JMTX_FAST_INT_T width = mtx->lower_bandwidth + 1 + mtx->upper_bandwidth;
    if (row < mtx->lower_bandwidth)
    {
        width -= (mtx->lower_bandwidth - row);
    }
    if (row > mtx->base.rows - 1 - mtx->upper_bandwidth)
    {
        width -= row - (mtx->base.rows - 1 - mtx->upper_bandwidth);
    }
    return width;
}

jmtx_result JMTX_NAME_TYPED(matrix_brm_new)(JMTX_NAME_TYPED(matrix_brm) * *p_mtx, JMTX_INDEX_T rows, JMTX_INDEX_T cols,
                                            JMTX_INDEX_T ubw, JMTX_INDEX_T lbw, const JMTX_SCALAR_T *set_value,
                                            const jmtx_allocator_callbacks *allocator_callbacks)
{
    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }

    ;

    JMTX_NAME_TYPED(matrix_brm) *mtx = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*mtx));
    if (!mtx)
    {
        return JMTX_RESULT_BAD_ALLOC;
    }

    memset(mtx, 0xCC, sizeof *mtx);
    mtx->base.cols = cols;
    mtx->base.type = JMTXD_TYPE_BRM;
    mtx->base.rows = rows;
    mtx->lower_bandwidth = lbw;
    mtx->upper_bandwidth = ubw;
    mtx->base.allocator_callbacks = *allocator_callbacks;
    const JMTX_SIZE_T entry_count = brm_row_offset(mtx, rows);

    JMTX_SCALAR_T *values = allocator_callbacks->alloc(allocator_callbacks->state, entry_count * sizeof(*values));
    if (!values)
    {
        allocator_callbacks->free(allocator_callbacks->state, mtx);
        return JMTX_RESULT_BAD_ALLOC;
    }
    if (set_value)
    {
        const JMTX_SCALAR_T v = *set_value;
        if (v == 0.0f)
        {
            memset(values, 0, sizeof(*values) * entry_count);
        }
        else
        {
            for (JMTX_SIZE_T i = 0; i < entry_count; ++i)
            {
                values[i] = v;
            }
        }
    }
    mtx->values = values;

    *p_mtx = mtx;

    return JMTX_RESULT_SUCCESS;
}

void JMTX_NAME_TYPED(matrix_brm_destroy)(JMTX_NAME_TYPED(matrix_brm) * mtx)
{
    const jmtx_allocator_callbacks allocator = mtx->base.allocator_callbacks;
    allocator.free(allocator.state, mtx->values);
    allocator.free(allocator.state, mtx);
}

void JMTX_NAME_TYPED(matrix_brm_set_row)(const JMTX_NAME_TYPED(matrix_brm) * mtx, JMTX_INDEX_T row,
                                         const JMTX_SCALAR_T values[])
{
    const JMTX_FAST_INT_T offset = brm_row_offset(mtx, row);
    const JMTX_FAST_INT_T len = brm_row_len(mtx, row);
    memcpy(mtx->values + offset, values, len * sizeof(*mtx->values));
}

void JMTX_NAME_TYPED(matrix_brm_vector_multiply)(const JMTX_NAME_TYPED(matrix_brm) * mtx,
                                                 const JMTX_SCALAR_T *restrict x, JMTX_SCALAR_T *restrict y)
{
    for (JMTX_INDEX_T i = 0; i < mtx->base.rows; ++i)
    {
        const JMTX_SCALAR_T *const values = mtx->values + brm_row_offset(mtx, i);
        const JMTX_FAST_INT_T first_elm = JMTX_NAME_TYPED(matrix_brm_first_pos_in_row)(mtx, i);
        const JMTX_FAST_INT_T last_elm = JMTX_NAME_TYPED(matrix_brm_last_pos_in_row)(mtx, i);

        JMTX_SCALAR_T v = 0;
        for (JMTX_INDEX_T j = 0; j < last_elm - first_elm + 1; ++j)
        {
            v += values[j] * x[j + first_elm];
        }
        y[i] = v;
    }
}

void JMTX_NAME_TYPED(matrix_brm_set_entry)(const JMTX_NAME_TYPED(matrix_brm) * mtx, JMTX_INDEX_T i, JMTX_INDEX_T j,
                                           JMTX_SCALAR_T value)
{
    JMTX_SCALAR_T *const values = mtx->values + brm_row_offset(mtx, i);
    const JMTX_FAST_INT_T offset = j - JMTX_NAME_TYPED(matrix_brm_first_pos_in_row)(mtx, i);
    assert(j >= JMTX_NAME_TYPED(matrix_brm_first_pos_in_row)(mtx, i));
    assert(j <= JMTX_NAME_TYPED(matrix_brm_last_pos_in_row)(mtx, i));
    values[offset] = value;
}

JMTX_SCALAR_T JMTX_NAME_TYPED(matrix_brm_get_entry)(const JMTX_NAME_TYPED(matrix_brm) * mtx, JMTX_INDEX_T i,
                                                    JMTX_INDEX_T j)
{
    if ((j < JMTX_NAME_TYPED(matrix_brm_first_pos_in_row)(mtx, i)) ||
        (j > JMTX_NAME_TYPED(matrix_brm_last_pos_in_row)(mtx, i)))
    {
        return 0.0f;
    }
    const JMTX_SCALAR_T *const values = mtx->values + brm_row_offset(mtx, i);
    const JMTX_FAST_INT_T offset = j - JMTX_NAME_TYPED(matrix_brm_first_pos_in_row)(mtx, i);

    return values[offset];
}

JMTX_FAST_INT_T JMTX_NAME_TYPED(matrix_brm_get_row)(const JMTX_NAME_TYPED(matrix_brm) * mtx, JMTX_INDEX_T row,
                                                    JMTX_SCALAR_T *p_elements[1])
{
    *p_elements = mtx->values + brm_row_offset(mtx, row);
    return brm_row_len(mtx, row);
}

JMTX_INDEX_T JMTX_NAME_TYPED(matrix_brm_count_values)(const JMTX_NAME_TYPED(matrix_brm) * mtx, JMTX_SCALAR_T v)
{
    JMTX_INDEX_T r = 0;
    const JMTX_FAST_INT_T count = brm_row_offset(mtx, mtx->base.rows);
    for (JMTX_INDEX_T i = 0; i < count; ++i)
    {
        if (v == mtx->values[i])
        {
            r += 1;
        }
    }
    return r;
}

JMTX_INDEX_T JMTX_NAME_TYPED(matrix_brm_length_of_col)(const JMTX_NAME_TYPED(matrix_brm) * mtx, JMTX_INDEX_T col)
{
    return JMTX_NAME_TYPED(matrix_brm_last_pos_in_col)(mtx, col) -
           JMTX_NAME_TYPED(matrix_brm_first_pos_in_col)(mtx, col) + 1;
}

JMTX_INDEX_T JMTX_NAME_TYPED(matrix_brm_get_col)(const JMTX_NAME_TYPED(matrix_brm) * mtx, JMTX_INDEX_T col,
                                                 JMTX_SCALAR_T values[])
{
    JMTX_FAST_INT_T first_row;
    if (col < mtx->upper_bandwidth)
    {
        first_row = 0;
    }
    else
    {
        first_row = col - mtx->upper_bandwidth;
    }
    JMTX_FAST_INT_T last_row = col + 1 + mtx->lower_bandwidth;
    if (last_row > mtx->base.rows)
    {
        last_row = mtx->base.rows;
    }

    JMTX_FAST_INT_T row;
    JMTX_FAST_INT_T j, i, pos_rel;
#ifndef NDEBUG
    const JMTX_FAST_INT_T max = brm_row_offset(mtx, mtx->base.rows);
#endif
    pos_rel = (col - JMTX_NAME_TYPED(matrix_brm_first_pos_in_row)(mtx, first_row));
    i = brm_row_offset(mtx, first_row) + pos_rel;
    for (row = first_row, j = 0; row < last_row; ++row)
    {
        assert(i < max);
        values[j++] = mtx->values[i];
        const JMTX_FAST_INT_T new_pos_rel = (col - JMTX_NAME_TYPED(matrix_brm_first_pos_in_row)(mtx, row + 1));
        i += (brm_row_len(mtx, row) - (pos_rel)) + new_pos_rel;
        pos_rel = new_pos_rel;
    }

    return j;
}

jmtx_result JMTX_NAME_TYPED(matrix_brm_transpose)(const JMTX_NAME_TYPED(matrix_brm) * mtx,
                                                  JMTX_NAME_TYPED(matrix_brm) * *p_out,
                                                  const jmtx_allocator_callbacks *allocator_callbacks)
{
    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }

    const JMTX_INDEX_T n_elements = brm_row_offset(mtx, mtx->base.rows);
    const JMTX_INDEX_T new_rows = mtx->base.cols;
    const JMTX_INDEX_T new_cols = mtx->base.rows;

    JMTX_NAME_TYPED(matrix_brm) *const out = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*out));
    if (!out)
    {
        return JMTX_RESULT_BAD_ALLOC;
    }

    JMTX_SCALAR_T *const new_values =
        allocator_callbacks->alloc(allocator_callbacks->state, (n_elements) * sizeof(*new_values));
    if (!new_values)
    {
        allocator_callbacks->free(allocator_callbacks->state, out);
        return JMTX_RESULT_BAD_ALLOC;
    }

    for (JMTX_INDEX_T j = 0, n = 0; j < mtx->base.cols; ++j)
    {
        n += JMTX_NAME_TYPED(matrix_brm_get_col)(mtx, j, new_values + n);
    }

    out->base.type = JMTXD_TYPE_BRM;
    out->lower_bandwidth = mtx->upper_bandwidth;
    out->upper_bandwidth = mtx->lower_bandwidth;
    out->values = new_values;
    out->base = mtx->base;
    out->base.rows = new_rows;
    out->base.cols = new_cols;
    out->base.allocator_callbacks = *allocator_callbacks;
    *p_out = out;

    return JMTX_RESULT_SUCCESS;
}

jmtx_result JMTX_NAME_TYPED(matrix_brm_copy)(const JMTX_NAME_TYPED(matrix_brm) * mtx,
                                             JMTX_NAME_TYPED(matrix_brm) * *p_out,
                                             const jmtx_allocator_callbacks *allocator_callbacks)
{
    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }
    JMTX_NAME_TYPED(matrix_brm) *const out = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*out));
    if (!out)
    {
        return JMTX_RESULT_BAD_ALLOC;
    }
    const JMTX_FAST_INT_T n_entries = brm_row_offset(mtx, mtx->base.rows);
    JMTX_SCALAR_T *const elements =
        allocator_callbacks->alloc(allocator_callbacks->state, (n_entries) * sizeof(*elements));
    if (!elements)
    {
        allocator_callbacks->free(allocator_callbacks->state, out);
        return JMTX_RESULT_BAD_ALLOC;
    }

    memcpy(elements, mtx->values, sizeof *elements * n_entries);

    out->base.type = JMTXD_TYPE_BRM;
    out->lower_bandwidth = mtx->lower_bandwidth;
    out->upper_bandwidth = mtx->upper_bandwidth;
    out->values = elements;
    out->base = mtx->base;
    out->base.rows = mtx->base.rows;
    out->base.cols = mtx->base.cols;
    out->base.allocator_callbacks = *allocator_callbacks;
    *p_out = out;
    return JMTX_RESULT_SUCCESS;
}

JMTX_SCALAR_T JMTX_NAME_TYPED(matrix_brm_vector_multiply_row)(const JMTX_NAME_TYPED(matrix_brm) * mtx,
                                                              const JMTX_SCALAR_T *x, JMTX_INDEX_T i)
{
    const JMTX_FAST_INT_T offset = brm_row_offset(mtx, i);
    const JMTX_FAST_INT_T first_col = JMTX_NAME_TYPED(matrix_brm_first_pos_in_row)(mtx, i);
    const JMTX_FAST_INT_T len = brm_row_len(mtx, i);
    JMTX_SCALAR_T v = 0;
    const JMTX_SCALAR_T *const values = mtx->values + offset;
    for (JMTX_INDEX_T j = 0; j < len; ++j)
    {
        v += values[j] * x[first_col + j];
    }
    return v;
}

void JMTX_NAME_TYPED(matrix_brm_add_to_entry)(const JMTX_NAME_TYPED(matrix_brm) * mtx, JMTX_INDEX_T i, JMTX_INDEX_T j,
                                              JMTX_SCALAR_T value)
{
    JMTX_FAST_INT_T row_offset = brm_row_offset(mtx, i);
    JMTX_FAST_INT_T col_offset = j - JMTX_NAME_TYPED(matrix_brm_first_pos_in_row)(mtx, i);

    mtx->values[row_offset + col_offset] += value;
}

void JMTX_NAME_TYPED(matrix_brm_zero_all_entries)(const JMTX_NAME_TYPED(matrix_brm) * mtx)
{
    memset(mtx->values, 0, sizeof(*mtx->values) * brm_row_offset(mtx, mtx->base.rows));
}

void JMTX_NAME_TYPED(matrix_brm_set_all_entries)(const JMTX_NAME_TYPED(matrix_brm) * mtx, JMTX_SCALAR_T x)
{
    for (JMTX_SCALAR_T *ptr = mtx->values; ptr != mtx->values + brm_row_offset(mtx, mtx->base.rows); ++ptr)
    {
        *ptr = x;
    }
}

JMTX_HOT_FUNCTION
JMTX_FAST_INT_T JMTX_NAME_TYPED(matrix_brm_first_pos_in_row)(const JMTX_NAME_TYPED(matrix_brm) * mtx, JMTX_INDEX_T row)
{
    if (row < mtx->lower_bandwidth)
    {
        return 0;
    }
    return row - mtx->lower_bandwidth;
}

JMTX_FAST_INT_T JMTX_NAME_TYPED(matrix_brm_last_pos_in_row)(const JMTX_NAME_TYPED(matrix_brm) * mtx, JMTX_INDEX_T row)
{
    if (row > mtx->base.rows - 1 - mtx->upper_bandwidth)
    {
        return mtx->base.rows - 1;
    }
    return row + mtx->upper_bandwidth;
}

JMTX_FAST_INT_T JMTX_NAME_TYPED(matrix_brm_first_pos_in_col)(const JMTX_NAME_TYPED(matrix_brm) * mtx, JMTX_INDEX_T col)
{
    if (col < mtx->upper_bandwidth)
    {
        return 0;
    }
    return col - mtx->upper_bandwidth;
}

JMTX_FAST_INT_T JMTX_NAME_TYPED(matrix_brm_last_pos_in_col)(const JMTX_NAME_TYPED(matrix_brm) * mtx, JMTX_INDEX_T col)
{
    if (col > mtx->base.rows - 1 - mtx->lower_bandwidth)
    {
        return mtx->base.cols - 1;
    }
    return col + mtx->lower_bandwidth;
}

JMTX_FAST_INT_T JMTX_NAME_TYPED(matrix_brm_length_of_row)(const JMTX_NAME_TYPED(matrix_brm) * mtx, JMTX_INDEX_T row)
{
    return brm_row_len(mtx, row);
}

void JMTX_NAME_TYPED(matrix_brm_set_col)(const JMTX_NAME_TYPED(matrix_brm) * mtx, JMTX_INDEX_T col,
                                         const JMTX_SCALAR_T values[])
{
    JMTX_FAST_INT_T first_row;
    if (col < mtx->upper_bandwidth)
    {
        first_row = 0;
    }
    else
    {
        first_row = col - mtx->upper_bandwidth;
    }
    JMTX_FAST_INT_T last_row = col + 1 + mtx->lower_bandwidth;
    if (last_row > mtx->base.rows)
    {
        last_row = mtx->base.rows;
    }

    JMTX_FAST_INT_T row;
    JMTX_FAST_INT_T j, i, pos_rel;
#ifndef NDEBUG
    const JMTX_FAST_INT_T max = brm_row_offset(mtx, mtx->base.rows);
#endif
    pos_rel = (col - JMTX_NAME_TYPED(matrix_brm_first_pos_in_row)(mtx, first_row));
    i = brm_row_offset(mtx, first_row) + pos_rel;
    for (row = first_row, j = 0; row < last_row; ++row)
    {
        assert(i < max);
        mtx->values[i] = values[j++];
        const JMTX_FAST_INT_T new_pos_rel = (col - JMTX_NAME_TYPED(matrix_brm_first_pos_in_row)(mtx, row + 1));
        i += (brm_row_len(mtx, row) - (pos_rel)) + new_pos_rel;
        pos_rel = new_pos_rel;
    }
}

void JMTX_NAME_TYPED(matrix_brm_get_bandwidths)(const JMTX_NAME_TYPED(matrix_brm) * mtx, JMTX_INDEX_T *ubw,
                                                JMTX_INDEX_T *lbw)
{
    *ubw = mtx->upper_bandwidth;
    *lbw = mtx->lower_bandwidth;
}
