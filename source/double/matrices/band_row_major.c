// Automatically generated from source/float/matrices/band_row_major.c on Fri Dec  1 06:43:01 2023
//
// Created by jan on 13.6.2022.
//

#include <assert.h>
#include <math.h>
#include "../../../include/jmtx/double/matrices/band_row_major.h"
#include "band_row_major_internal.h"
#include "../../../include/jmtx/double/matrices/band_row_major_safe.h"


static inline uint_fast32_t brm_row_offset(const jmtxd_matrix_brm mtx[const static 1], const uint_fast32_t row)
{
    uint_fast32_t offset = (mtx->lower_bandwidth + 1 + mtx->upper_bandwidth) * row;
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
        const uint_fast32_t offset_row = row - (mtx->base.rows - 1 - mtx->upper_bandwidth);   //  rows offset from
        offset -= offset_row * (offset_row - 1) / 2;
    }
    return offset;
}

static inline uint_fast32_t brm_row_len(const jmtxd_matrix_brm mtx[const static 1], const uint_fast32_t row)
{
    uint_fast32_t width = mtx->lower_bandwidth + 1 + mtx->upper_bandwidth;
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


jmtx_result jmtxd_matrix_brm_new(
    jmtxd_matrix_brm** p_mtx, uint32_t rows, uint32_t cols, uint32_t ubw, uint32_t lbw, const double* set_value,
    const jmtx_allocator_callbacks* allocator_callbacks)
{
    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }

    ;

    jmtxd_matrix_brm* mtx = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*mtx));
    if (!mtx)
    {
        return JMTX_RESULT_BAD_ALLOC;
    }


    memset(mtx, 0xCC, sizeof*mtx);
    mtx->base.cols = cols;
    mtx->base.type = JMTXD_TYPE_BRM;
    mtx->base.rows = rows;
    mtx->lower_bandwidth = lbw;
    mtx->upper_bandwidth = ubw;
    mtx->base.allocator_callbacks = *allocator_callbacks;
    const uint64_t entry_count = brm_row_offset(mtx, rows);

    double* values = allocator_callbacks->alloc(allocator_callbacks->state, entry_count * sizeof(*values));
    if (!values)
    {
        allocator_callbacks->free(allocator_callbacks->state, mtx);
        return JMTX_RESULT_BAD_ALLOC;
    }
    if (set_value)
    {
        const double v = *set_value;
        if (v == 0.0f)
        {
            memset(values, 0, sizeof(*values) * entry_count);
        }
        else
        {
            for (uint64_t i = 0; i < entry_count; ++i)
            {
                values[i] = v;
            }
        }
    }
    mtx->values = values;

    *p_mtx = mtx;

    return JMTX_RESULT_SUCCESS;
}

jmtx_result jmtxds_matrix_brm_new(
    jmtxd_matrix_brm** p_mtx, uint32_t rows, uint32_t cols, uint32_t ubw, uint32_t lbw, const double* set_value,
    const jmtx_allocator_callbacks* allocator_callbacks)
{
    if (!p_mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!rows)
    {
        return JMTX_RESULT_BAD_PARAM;
    }
    if (!cols)
    {
        return JMTX_RESULT_BAD_PARAM;
    }
    if (ubw > cols)
    {
        return JMTX_RESULT_BAD_PARAM;
    }
    if (lbw > rows)
    {
        return JMTX_RESULT_BAD_PARAM;
    }
    if (allocator_callbacks && (!allocator_callbacks->alloc || !allocator_callbacks->realloc || !allocator_callbacks->free))
    {
        return JMTX_RESULT_BAD_PARAM;
    }
    if (set_value && !isfinite(*set_value))
    {
        return JMTX_RESULT_BAD_PARAM;
    }

    return jmtxd_matrix_brm_new(p_mtx, rows, cols, ubw, lbw, set_value, allocator_callbacks);
}

void jmtxd_matrix_brm_destroy(jmtxd_matrix_brm* mtx)
{
    jmtx_allocator_callbacks allocator = mtx->base.allocator_callbacks;
    allocator.free(allocator.state, mtx->values);
    allocator.free(allocator.state, mtx);
}

jmtx_result jmtxds_matrix_brm_destroy(jmtxd_matrix_brm* mtx)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTXD_TYPE_BRM)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    jmtxd_matrix_brm_destroy(mtx);
    return JMTX_RESULT_SUCCESS;
}

void jmtxd_matrix_brm_set_row(const jmtxd_matrix_brm* mtx, uint32_t row, double values[])
{
    const uint_fast32_t offset = brm_row_offset(mtx, row);
    const uint_fast32_t len = brm_row_len(mtx, row);
    memcpy(mtx->values + offset, values, len * sizeof(*mtx->values));
}

jmtx_result jmtxds_matrix_brm_set_row(const jmtxd_matrix_brm* mtx, uint32_t row, double values[])
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTXD_TYPE_BRM)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (mtx->base.rows <= row)
    {
        return JMTX_RESULT_INDEX_OUT_OF_BOUNDS;
    }
    if (!values)
    {
        return JMTX_RESULT_NULL_PARAM;
    }

    const uint_fast32_t len = brm_row_len(mtx, row);
    for (uint32_t i = 0; i < len; ++i)
    {
        if (!isfinite(values[i]))
        {
            return JMTX_RESULT_BAD_PARAM;
        }
    }
    
    jmtxd_matrix_brm_set_row(mtx, row, values);
    return JMTX_RESULT_SUCCESS;
}


void jmtxd_matrix_brm_vector_multiply(const jmtxd_matrix_brm* mtx, const double* restrict x, double* restrict y)
{
    for (uint32_t i = 0; i < mtx->base.rows; ++i)
    {
        const double* const values = mtx->values + brm_row_offset(mtx, i);
        const uint_fast32_t first_elm = jmtxd_matrix_brm_first_pos_in_row(mtx, i);
        const uint_fast32_t last_elm = jmtxd_matrix_brm_last_pos_in_row(mtx, i);

        double v = 0;
        for (uint32_t j = 0; j < last_elm - first_elm + 1; ++j)
        {
            v += values[j] * x[j + first_elm];
        }
        y[i] = v;
    }
}

jmtx_result jmtxds_matrix_brm_vector_multiply(const jmtxd_matrix_brm* mtx, const double* restrict x, double* restrict y)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTXD_TYPE_BRM)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (!x)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!y)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    jmtxd_matrix_brm_vector_multiply(mtx, x, y);
    return JMTX_RESULT_SUCCESS;
}

void jmtxd_matrix_brm_set_entry(const jmtxd_matrix_brm* mtx, uint32_t i, uint32_t j, double value)
{
    double* const values = mtx->values + brm_row_offset(mtx, i);
    const uint_fast32_t offset = j - jmtxd_matrix_brm_first_pos_in_row(mtx, i);
    assert(j >= jmtxd_matrix_brm_first_pos_in_row(mtx, i));
    assert(j <= jmtxd_matrix_brm_last_pos_in_row(mtx, i));
    values[offset] = value;
}

jmtx_result jmtxds_matrix_brm_set_entry(const jmtxd_matrix_brm* mtx, uint32_t i, uint32_t j, double value)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTXD_TYPE_BRM)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (j >= mtx->base.cols)
    {
        return JMTX_RESULT_INDEX_OUT_OF_BOUNDS;
    }
    if (i >= mtx->base.rows)
    {
        return JMTX_RESULT_INDEX_OUT_OF_BOUNDS;
    }
    if (!isfinite(value))
    {
        return JMTX_RESULT_BAD_PARAM;
    }
    if (j < jmtxd_matrix_brm_first_pos_in_row(mtx, i))
    {
        return JMTX_RESULT_INDEX_OUT_OF_BOUNDS;
    }
    if (j > jmtxd_matrix_brm_last_pos_in_row(mtx, i))
    {
        return JMTX_RESULT_INDEX_OUT_OF_BOUNDS;
    }
    jmtxd_matrix_brm_set_entry(mtx, i, j, value);
    return JMTX_RESULT_SUCCESS;
}

double jmtxd_matrix_brm_get_entry(const jmtxd_matrix_brm* mtx, uint32_t i, uint32_t j)
{
    if ((j < jmtxd_matrix_brm_first_pos_in_row(mtx, i)) || (j > jmtxd_matrix_brm_last_pos_in_row(mtx, i)))
    {
        return 0.0f;
    }
    double* const values = mtx->values + brm_row_offset(mtx, i);
    const uint_fast32_t offset = j - jmtxd_matrix_brm_first_pos_in_row(mtx, i);

    return values[offset];
}

jmtx_result jmtxds_matrix_brm_get_entry(const jmtxd_matrix_brm* mtx, uint32_t i, uint32_t j, double* p_value)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTXD_TYPE_BRM)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (j >= mtx->base.cols)
    {
        return JMTX_RESULT_INDEX_OUT_OF_BOUNDS;
    }
    if (i >= mtx->base.rows)
    {
        return JMTX_RESULT_INDEX_OUT_OF_BOUNDS;
    }
    *p_value = jmtxd_matrix_brm_get_entry(mtx, i, j);

    return JMTX_RESULT_SUCCESS;
}

uint_fast32_t jmtxd_matrix_brm_get_row(const jmtxd_matrix_brm* mtx, uint32_t row, double* p_elements[1])
{
    *p_elements = mtx->values + brm_row_offset(mtx, row);
    return brm_row_len(mtx, row);
}

jmtx_result jmtxds_matrix_brm_get_row(const jmtxd_matrix_brm* mtx, uint_fast32_t* n, uint32_t row, double* p_elements[1])
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTXD_TYPE_BRM)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (row >= mtx->base.rows)
    {
        return JMTX_RESULT_INDEX_OUT_OF_BOUNDS;
    }
    if (!n)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!p_elements)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    *n = jmtxd_matrix_brm_get_row(mtx, row, p_elements);
    return JMTX_RESULT_SUCCESS;
}

uint32_t jmtxd_matrix_brm_count_values(const jmtxd_matrix_brm* mtx, double v)
{
    uint32_t r = 0;
    const uint_fast32_t count = brm_row_offset(mtx, mtx->base.rows);
    for (uint32_t i = 0; i < count; ++i)
    {
        if (v == mtx->values[i])
        {
            r += 1;
        }
    }
    return r;
}

jmtx_result jmtxds_matrix_brm_count_values(const jmtxd_matrix_brm* mtx, double v, uint32_t* p_count)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTXD_TYPE_BRM)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (!p_count)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }

    *p_count = jmtxd_matrix_brm_count_values(mtx, v);
    return JMTX_RESULT_SUCCESS;
}

//static inline uint_fast32_t entries_above_diagonal(const jmtxd_matrix_brm mtx[const static 1], uint_fast32_t col)
//{
//    if (col < mtx->base.rows -  mtx->upper_bandwidth - 1)
//    {
//        return mtx->upper_bandwidth;
//    }
//    return mtx->upper_bandwidth - (mtx->base.rows -  mtx->upper_bandwidth - col);
//}
//
//static inline uint_fast32_t entries_bellow_diagonal(const jmtxd_matrix_brm mtx[const static 1], uint_fast32_t col)
//{
//    if (col > mtx->lower_bandwidth)
//    {
//        return col - 1;
//    }
//    return mtx->lower_bandwidth;
//}

uint32_t jmtxd_matrix_brm_entries_in_col(const jmtxd_matrix_brm* mtx, uint32_t col)
{
    return jmtxd_matrix_brm_last_pos_in_col(mtx, col) - jmtxd_matrix_brm_first_pos_in_col(mtx, col) + 1;
}

jmtx_result jmtxds_matrix_brm_entries_in_col(const jmtxd_matrix_brm* mtx, uint32_t col, uint32_t* p_n)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTXD_TYPE_BRM)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (col >= mtx->base.cols)
    {
        return JMTX_RESULT_INDEX_OUT_OF_BOUNDS;
    }
    if (!p_n)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    *p_n = jmtxd_matrix_brm_entries_in_col(mtx, col);
    return JMTX_RESULT_SUCCESS;
}

uint32_t
jmtxd_matrix_brm_get_col(const jmtxd_matrix_brm* mtx, uint32_t col, double values[])
{
    uint_fast32_t first_row;
    if (col < mtx->upper_bandwidth)
    {
        first_row = 0;
    }
    else
    {
        first_row = col - mtx->upper_bandwidth;
    }
    uint_fast32_t last_row = col + 1 + mtx->lower_bandwidth;
    if (last_row > mtx->base.rows)
    {
        last_row = mtx->base.rows;
    }

    uint_fast32_t row;
    uint_fast32_t j, i, pos_rel;
    const uint_fast32_t max = brm_row_offset(mtx, mtx->base.rows);
    pos_rel = (col - jmtxd_matrix_brm_first_pos_in_row(mtx, first_row));
    i = brm_row_offset(mtx, first_row) + pos_rel;
    for (row = first_row, j = 0; row < last_row; ++row)
    {
        assert(i < max);
        values[j++] = mtx->values[i];
        uint_fast32_t new_pos_rel = (col - jmtxd_matrix_brm_first_pos_in_row(mtx, row + 1));
        i += (brm_row_len(mtx, row) - (pos_rel))
                +  new_pos_rel;
        pos_rel = new_pos_rel;
    }

    return j;
}

jmtx_result jmtxds_matrix_brm_get_col(
        const jmtxd_matrix_brm* mtx, uint32_t col, uint32_t* p_count, double p_values[])
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTXD_TYPE_BRM)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (col >= mtx->base.cols)
    {
        return JMTX_RESULT_INDEX_OUT_OF_BOUNDS;
    }
    if (p_values == NULL)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!p_count)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    *p_count = jmtxd_matrix_brm_get_col(mtx, col, p_values);
    return JMTX_RESULT_SUCCESS;
}

jmtx_result jmtxd_matrix_brm_transpose(
        const jmtxd_matrix_brm* mtx, jmtxd_matrix_brm** p_out, const jmtx_allocator_callbacks* allocator_callbacks)
{
    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &mtx->base.allocator_callbacks;
    }

    const uint32_t n_elements = brm_row_offset(mtx, mtx->base.rows);
    const uint32_t new_rows = mtx->base.cols;
    const uint32_t new_cols = mtx->base.rows;

    jmtxd_matrix_brm* const out = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*out));
    if (!out)
    {
        return JMTX_RESULT_BAD_ALLOC;
    }


    double* const new_values = allocator_callbacks->alloc(allocator_callbacks->state, (n_elements) * sizeof(*new_values));
    if (!new_values)
    {
        allocator_callbacks->free(allocator_callbacks->state, out);
        return JMTX_RESULT_BAD_ALLOC;
    }


    for (uint32_t j = 0, n = 0; j < mtx->base.cols; ++j)
    {
        n += jmtxd_matrix_brm_get_col(mtx, j, new_values + n);
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

jmtx_result jmtxds_matrix_brm_transpose(
        const jmtxd_matrix_brm* mtx, jmtxd_matrix_brm** p_out, const jmtx_allocator_callbacks* allocator_callbacks)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTXD_TYPE_BRM)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (!p_out)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (allocator_callbacks && (!allocator_callbacks->alloc || !allocator_callbacks->realloc || !allocator_callbacks->free))
    {
        return JMTX_RESULT_BAD_PARAM;
    }
    return jmtxd_matrix_brm_transpose(mtx, p_out, allocator_callbacks);
}

jmtx_result jmtxd_matrix_brm_copy(const jmtxd_matrix_brm* mtx, jmtxd_matrix_brm** p_out, const jmtx_allocator_callbacks* allocator_callbacks)
{
    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &mtx->base.allocator_callbacks;
    }
    jmtxd_matrix_brm* const out = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*out));
    if (!out)
    {
        return JMTX_RESULT_BAD_ALLOC;
    }
    const uint_fast32_t n_entries = brm_row_offset(mtx, mtx->base.rows);
    double* const elements = allocator_callbacks->alloc(allocator_callbacks->state, (n_entries) * sizeof (*elements));
    if (!elements)
    {
        allocator_callbacks->free(allocator_callbacks->state, out);
        return JMTX_RESULT_BAD_ALLOC;
    }

    memcpy(elements, mtx->values, sizeof* elements * n_entries);

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

jmtx_result jmtxds_matrix_brm_copy(const jmtxd_matrix_brm* mtx, jmtxd_matrix_brm** p_out, const jmtx_allocator_callbacks* allocator_callbacks)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTXD_TYPE_BRM)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (!p_out)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (allocator_callbacks && (!allocator_callbacks->alloc || !allocator_callbacks->realloc || !allocator_callbacks->free))
    {
        return JMTX_RESULT_BAD_PARAM;
    }
    return jmtxd_matrix_brm_copy(mtx, p_out, allocator_callbacks);
}

double jmtxd_matrix_brm_vector_multiply_row(const jmtxd_matrix_brm* mtx, const double* x, uint32_t i)
{
    const uint_fast32_t offset = brm_row_offset(mtx, i);
    const uint_fast32_t first_col = jmtxd_matrix_brm_first_pos_in_row(mtx, i);
    const uint_fast32_t len = brm_row_len(mtx, i);
    double v = 0;
    const double* const values = mtx->values + offset;
    for (uint32_t j = 0; j < len; ++j)
    {
        v += values[j] * x[first_col + j];
    }
    return v;
}

jmtx_result jmtxds_matrix_brm_vector_multiply_row(const jmtxd_matrix_brm* mtx, const double* restrict x, uint32_t i, double* restrict p_r)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTXD_TYPE_BRM)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (!x)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (i >= mtx->base.rows)
    {
        return JMTX_RESULT_INDEX_OUT_OF_BOUNDS;
    }
    if (!p_r)
    {
        return JMTX_RESULT_NULL_PARAM;
    }

    *p_r = jmtxd_matrix_brm_vector_multiply_row(mtx, x, i);
    return JMTX_RESULT_SUCCESS;
}

void jmtxd_matrix_brm_add_to_entry(const jmtxd_matrix_brm* mtx, uint32_t i, uint32_t j, double value)
{
    uint_fast32_t row_offset = brm_row_offset(mtx, i);
    uint_fast32_t col_offset = j - jmtxd_matrix_brm_first_pos_in_row(mtx, i);

    mtx->values[row_offset + col_offset] += value;
}

jmtx_result jmtxds_matrix_brm_add_to_entry(const jmtxd_matrix_brm* mtx, uint32_t i, uint32_t j, double value)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTXD_TYPE_BRM)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (j >= mtx->base.cols)
    {
        return JMTX_RESULT_INDEX_OUT_OF_BOUNDS;
    }
    if (i >= mtx->base.rows)
    {
        return JMTX_RESULT_INDEX_OUT_OF_BOUNDS;
    }
    jmtxd_matrix_brm_add_to_entry(mtx, i, j, value);
    return JMTX_RESULT_SUCCESS;
}

void jmtxd_matrix_brm_zero_all_entries(const jmtxd_matrix_brm* mtx)
{
    memset(mtx->values, 0, sizeof(*mtx->values) * brm_row_offset(mtx, mtx->base.rows));
}

jmtx_result jmtxds_matrix_brm_zero_all_entries(const jmtxd_matrix_brm* mtx)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTXD_TYPE_BRM)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    jmtxd_matrix_brm_zero_all_entries(mtx);
    return JMTX_RESULT_SUCCESS;
}

void jmtxd_matrix_brm_set_all_entries(const jmtxd_matrix_brm* mtx, double x)
{
    for (double* ptr = mtx->values; ptr != mtx->values + brm_row_offset(mtx, mtx->base.rows); ++ptr)
    {
        *ptr = x;
    }
}

JMTX_HOT_FUNCTION
uint_fast32_t jmtxd_matrix_brm_first_pos_in_row(const jmtxd_matrix_brm* mtx, uint32_t row)
{
    if (row < mtx->lower_bandwidth)
    {
        return 0;
    }
    return row - mtx->lower_bandwidth;
}

uint_fast32_t jmtxd_matrix_brm_last_pos_in_row(const jmtxd_matrix_brm* mtx, uint32_t row)
{
    if (row > mtx->base.rows - 1 - mtx->upper_bandwidth)
    {
        return mtx->base.rows - 1;
    }
    return row + mtx->upper_bandwidth;
}

uint_fast32_t jmtxd_matrix_brm_first_pos_in_col(const jmtxd_matrix_brm* mtx, uint32_t col)
{
    if (col < mtx->upper_bandwidth)
    {
        return 0;
    }
    return col - mtx->upper_bandwidth;
}

uint_fast32_t jmtxd_matrix_brm_last_pos_in_col(const jmtxd_matrix_brm* mtx, uint32_t col)
{
    if (col > mtx->base.rows - 1 - mtx->lower_bandwidth)
    {
        return mtx->base.cols - 1;
    }
    return col + mtx->lower_bandwidth;
}

uint_fast32_t jmtxd_matrix_brm_length_of_row(const jmtxd_matrix_brm* mtx, uint32_t row)
{
    return brm_row_len(mtx, row);
}

void jmtxd_matrix_brm_set_col(const jmtxd_matrix_brm* mtx, uint32_t col, const double* values)
{
    uint_fast32_t first_row;
    if (col < mtx->upper_bandwidth)
    {
        first_row = 0;
    }
    else
    {
        first_row = col - mtx->upper_bandwidth;
    }
    uint_fast32_t last_row = col + 1 + mtx->lower_bandwidth;
    if (last_row > mtx->base.rows)
    {
        last_row = mtx->base.rows;
    }

    uint_fast32_t row;
    uint_fast32_t j, i, pos_rel;
    const uint_fast32_t max = brm_row_offset(mtx, mtx->base.rows);
    pos_rel = (col - jmtxd_matrix_brm_first_pos_in_row(mtx, first_row));
    i = brm_row_offset(mtx, first_row) + pos_rel;
    for (row = first_row, j = 0; row < last_row; ++row)
    {
        assert(i < max);
        mtx->values[i] = values[j++];
        uint_fast32_t new_pos_rel = (col - jmtxd_matrix_brm_first_pos_in_row(mtx, row + 1));
        i += (brm_row_len(mtx, row) - (pos_rel))
             +  new_pos_rel;
        pos_rel = new_pos_rel;
    }
}

void jmtxd_matrix_brm_get_bandwidths(const jmtxd_matrix_brm* mtx, uint32_t* ubw, uint32_t* lbw)
{
    *ubw = mtx->upper_bandwidth;
    *lbw = mtx->lower_bandwidth;
}

jmtx_result jmtxds_matrix_brm_set_all_entries(const jmtxd_matrix_brm* mtx, double x)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTXD_TYPE_BRM)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    jmtxd_matrix_brm_set_all_entries(mtx, x);

    return JMTX_RESULT_SUCCESS;
}

jmtx_result jmtxds_matrix_brm_set_col(const jmtxd_matrix_brm* mtx, uint32_t col, const double* values)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTXD_TYPE_BRM)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (mtx->base.cols <= col)
    {
        return JMTX_RESULT_INDEX_OUT_OF_BOUNDS;
    }
    if (!values)
    {
        return JMTX_RESULT_NULL_PARAM;
    }

    const uint_fast32_t len = jmtxd_matrix_brm_entries_in_col(mtx, col);
    for (uint32_t i = 0; i < len; ++i)
    {
        if (!isfinite(values[i]))
        {
            return JMTX_RESULT_BAD_PARAM;
        }
    }

    jmtxd_matrix_brm_set_col(mtx, col, values);
    return JMTX_RESULT_SUCCESS;
}

