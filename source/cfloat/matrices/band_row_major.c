// Automatically generated from source/float/matrices/band_row_major.c on Fri Dec  1 17:36:03 2023
//
// Created by jan on 13.6.2022.
//

#include <assert.h>
#include <math.h>
#include "../../../include/jmtx/cfloat/matrices/band_row_major.h"
#include "band_row_major_internal.h"
#include "../../../include/jmtx/cfloat/matrices/band_row_major_safe.h"
#include <complex.h>


static inline uint_fast32_t brm_row_offset(const jmtxc_matrix_brm mtx[JMTX_ARRAY_ATTRIB(const static 1)], const uint_fast32_t row)
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

static inline uint_fast32_t brm_row_len(const jmtxc_matrix_brm mtx[JMTX_ARRAY_ATTRIB(const static 1)], const uint_fast32_t row)
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


jmtx_result jmtxc_matrix_brm_new(
    jmtxc_matrix_brm** p_mtx, uint32_t rows, uint32_t cols, uint32_t ubw, uint32_t lbw, const _Complex float* set_value,
    const jmtx_allocator_callbacks* allocator_callbacks)
{
    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }

    ;

    jmtxc_matrix_brm* mtx = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*mtx));
    if (!mtx)
    {
        return JMTX_RESULT_BAD_ALLOC;
    }


    memset(mtx, 0xCC, sizeof*mtx);
    mtx->base.cols = cols;
    mtx->base.type = JMTXC_TYPE_BRM;
    mtx->base.rows = rows;
    mtx->lower_bandwidth = lbw;
    mtx->upper_bandwidth = ubw;
    mtx->base.allocator_callbacks = *allocator_callbacks;
    const uint64_t entry_count = brm_row_offset(mtx, rows);

    _Complex float* values = allocator_callbacks->alloc(allocator_callbacks->state, entry_count * sizeof(*values));
    if (!values)
    {
        allocator_callbacks->free(allocator_callbacks->state, mtx);
        return JMTX_RESULT_BAD_ALLOC;
    }
    if (set_value)
    {
        const _Complex float v = *set_value;
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

jmtx_result jmtxcs_matrix_brm_new(
    jmtxc_matrix_brm** p_mtx, uint32_t rows, uint32_t cols, uint32_t ubw, uint32_t lbw, const _Complex float* set_value,
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
    if (set_value && (!isfinite(crealf(*set_value)) || !isfinite(cimagf(*set_value))))
    {
        return JMTX_RESULT_BAD_PARAM;
    }

    return jmtxc_matrix_brm_new(p_mtx, rows, cols, ubw, lbw, set_value, allocator_callbacks);
}

void jmtxc_matrix_brm_destroy(jmtxc_matrix_brm* mtx)
{
    jmtx_allocator_callbacks allocator = mtx->base.allocator_callbacks;
    allocator.free(allocator.state, mtx->values);
    allocator.free(allocator.state, mtx);
}

jmtx_result jmtxcs_matrix_brm_destroy(jmtxc_matrix_brm* mtx)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTXC_TYPE_BRM)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    jmtxc_matrix_brm_destroy(mtx);
    return JMTX_RESULT_SUCCESS;
}

void jmtxc_matrix_brm_set_row(const jmtxc_matrix_brm* mtx, uint32_t row, const _Complex float values[])
{
    const uint_fast32_t offset = brm_row_offset(mtx, row);
    const uint_fast32_t len = brm_row_len(mtx, row);
    memcpy(mtx->values + offset, values, len * sizeof(*mtx->values));
}

jmtx_result jmtxcs_matrix_brm_set_row(const jmtxc_matrix_brm* mtx, uint32_t row, const _Complex float values[])
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTXC_TYPE_BRM)
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
        if ((!isfinite(crealf(values[i])) || !isfinite(cimagf(values[i]))))
        {
            return JMTX_RESULT_BAD_PARAM;
        }
    }
    
    jmtxc_matrix_brm_set_row(mtx, row, values);
    return JMTX_RESULT_SUCCESS;
}


void jmtxc_matrix_brm_vector_multiply(const jmtxc_matrix_brm* mtx, const _Complex float* restrict x, _Complex float* restrict y)
{
    for (uint32_t i = 0; i < mtx->base.rows; ++i)
    {
        const _Complex float* const values = mtx->values + brm_row_offset(mtx, i);
        const uint_fast32_t first_elm = jmtxc_matrix_brm_first_pos_in_row(mtx, i);
        const uint_fast32_t last_elm = jmtxc_matrix_brm_last_pos_in_row(mtx, i);

        _Complex float v = 0;
        for (uint32_t j = 0; j < last_elm - first_elm + 1; ++j)
        {
            v += values[j] * x[j + first_elm];
        }
        y[i] = v;
    }
}

jmtx_result jmtxcs_matrix_brm_vector_multiply(const jmtxc_matrix_brm* mtx, const _Complex float* restrict x, _Complex float* restrict y)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTXC_TYPE_BRM)
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
    jmtxc_matrix_brm_vector_multiply(mtx, x, y);
    return JMTX_RESULT_SUCCESS;
}

void jmtxc_matrix_brm_set_entry(const jmtxc_matrix_brm* mtx, uint32_t i, uint32_t j, _Complex float value)
{
    _Complex float* const values = mtx->values + brm_row_offset(mtx, i);
    const uint_fast32_t offset = j - jmtxc_matrix_brm_first_pos_in_row(mtx, i);
    assert(j >= jmtxc_matrix_brm_first_pos_in_row(mtx, i));
    assert(j <= jmtxc_matrix_brm_last_pos_in_row(mtx, i));
    values[offset] = value;
}

jmtx_result jmtxcs_matrix_brm_set_entry(const jmtxc_matrix_brm* mtx, uint32_t i, uint32_t j, _Complex float value)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTXC_TYPE_BRM)
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
    if (!isfinite(crealf(value)) || !isfinite(cimagf(value)))
    {
        return JMTX_RESULT_BAD_PARAM;
    }
    if (j < jmtxc_matrix_brm_first_pos_in_row(mtx, i))
    {
        return JMTX_RESULT_INDEX_OUT_OF_BOUNDS;
    }
    if (j > jmtxc_matrix_brm_last_pos_in_row(mtx, i))
    {
        return JMTX_RESULT_INDEX_OUT_OF_BOUNDS;
    }
    jmtxc_matrix_brm_set_entry(mtx, i, j, value);
    return JMTX_RESULT_SUCCESS;
}

_Complex float jmtxc_matrix_brm_get_entry(const jmtxc_matrix_brm* mtx, uint32_t i, uint32_t j)
{
    if ((j < jmtxc_matrix_brm_first_pos_in_row(mtx, i)) || (j > jmtxc_matrix_brm_last_pos_in_row(mtx, i)))
    {
        return 0.0f;
    }
    _Complex float* const values = mtx->values + brm_row_offset(mtx, i);
    const uint_fast32_t offset = j - jmtxc_matrix_brm_first_pos_in_row(mtx, i);

    return values[offset];
}

jmtx_result jmtxcs_matrix_brm_get_entry(const jmtxc_matrix_brm* mtx, uint32_t i, uint32_t j, _Complex float* p_value)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTXC_TYPE_BRM)
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
    *p_value = jmtxc_matrix_brm_get_entry(mtx, i, j);

    return JMTX_RESULT_SUCCESS;
}

uint_fast32_t jmtxc_matrix_brm_get_row(const jmtxc_matrix_brm* mtx, uint32_t row, _Complex float* p_elements[1])
{
    *p_elements = mtx->values + brm_row_offset(mtx, row);
    return brm_row_len(mtx, row);
}

jmtx_result jmtxcs_matrix_brm_get_row(const jmtxc_matrix_brm* mtx, uint_fast32_t* n, uint32_t row, _Complex float* p_elements[1])
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTXC_TYPE_BRM)
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
    *n = jmtxc_matrix_brm_get_row(mtx, row, p_elements);
    return JMTX_RESULT_SUCCESS;
}

uint32_t jmtxc_matrix_brm_count_values(const jmtxc_matrix_brm* mtx, _Complex float v)
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

jmtx_result jmtxcs_matrix_brm_count_values(const jmtxc_matrix_brm* mtx, _Complex float v, uint32_t* p_count)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTXC_TYPE_BRM)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (!p_count)
    {
        return JMTX_RESULT_NULL_PARAM;
    }

    *p_count = jmtxc_matrix_brm_count_values(mtx, v);
    return JMTX_RESULT_SUCCESS;
}

//static inline uint_fast32_t entries_above_diagonal(const jmtxc_matrix_brm mtx[JMTX_ARRAY_ATTRIB(const static 1)], uint_fast32_t col)
//{
//    if (col < mtx->base.rows -  mtx->upper_bandwidth - 1)
//    {
//        return mtx->upper_bandwidth;
//    }
//    return mtx->upper_bandwidth - (mtx->base.rows -  mtx->upper_bandwidth - col);
//}
//
//static inline uint_fast32_t entries_bellow_diagonal(const jmtxc_matrix_brm mtx[JMTX_ARRAY_ATTRIB(const static 1)], uint_fast32_t col)
//{
//    if (col > mtx->lower_bandwidth)
//    {
//        return col - 1;
//    }
//    return mtx->lower_bandwidth;
//}

uint32_t jmtxc_matrix_brm_length_of_col(const jmtxc_matrix_brm* mtx, uint32_t col)
{
    return jmtxc_matrix_brm_last_pos_in_col(mtx, col) - jmtxc_matrix_brm_first_pos_in_col(mtx, col) + 1;
}

jmtx_result jmtxcs_matrix_brm_length_of_col(const jmtxc_matrix_brm* mtx, uint32_t col, uint32_t* p_n)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTXC_TYPE_BRM)
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
    *p_n = jmtxc_matrix_brm_length_of_col(mtx, col);
    return JMTX_RESULT_SUCCESS;
}

uint32_t
jmtxc_matrix_brm_get_col(const jmtxc_matrix_brm* mtx, uint32_t col, _Complex float values[])
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
    pos_rel = (col - jmtxc_matrix_brm_first_pos_in_row(mtx, first_row));
    i = brm_row_offset(mtx, first_row) + pos_rel;
    for (row = first_row, j = 0; row < last_row; ++row)
    {
        assert(i < max);
        values[j++] = mtx->values[i];
        uint_fast32_t new_pos_rel = (col - jmtxc_matrix_brm_first_pos_in_row(mtx, row + 1));
        i += (brm_row_len(mtx, row) - (pos_rel))
                +  new_pos_rel;
        pos_rel = new_pos_rel;
    }

    return j;
}

jmtx_result jmtxcs_matrix_brm_get_col(
        const jmtxc_matrix_brm* mtx, uint32_t col, uint32_t* p_count, _Complex float p_values[])
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTXC_TYPE_BRM)
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
    *p_count = jmtxc_matrix_brm_get_col(mtx, col, p_values);
    return JMTX_RESULT_SUCCESS;
}

jmtx_result jmtxc_matrix_brm_transpose(
        const jmtxc_matrix_brm* mtx, jmtxc_matrix_brm** p_out, const jmtx_allocator_callbacks* allocator_callbacks)
{
    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }

    const uint32_t n_elements = brm_row_offset(mtx, mtx->base.rows);
    const uint32_t new_rows = mtx->base.cols;
    const uint32_t new_cols = mtx->base.rows;

    jmtxc_matrix_brm* const out = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*out));
    if (!out)
    {
        return JMTX_RESULT_BAD_ALLOC;
    }


    _Complex float* const new_values = allocator_callbacks->alloc(allocator_callbacks->state, (n_elements) * sizeof(*new_values));
    if (!new_values)
    {
        allocator_callbacks->free(allocator_callbacks->state, out);
        return JMTX_RESULT_BAD_ALLOC;
    }


    for (uint32_t j = 0, n = 0; j < mtx->base.cols; ++j)
    {
        n += jmtxc_matrix_brm_get_col(mtx, j, new_values + n);
    }

    out->base.type = JMTXC_TYPE_BRM;
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

jmtx_result jmtxcs_matrix_brm_transpose(
        const jmtxc_matrix_brm* mtx, jmtxc_matrix_brm** p_out, const jmtx_allocator_callbacks* allocator_callbacks)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTXC_TYPE_BRM)
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
    return jmtxc_matrix_brm_transpose(mtx, p_out, allocator_callbacks);
}

jmtx_result jmtxc_matrix_brm_copy(const jmtxc_matrix_brm* mtx, jmtxc_matrix_brm** p_out, const jmtx_allocator_callbacks* allocator_callbacks)
{
    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }
    jmtxc_matrix_brm* const out = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*out));
    if (!out)
    {
        return JMTX_RESULT_BAD_ALLOC;
    }
    const uint_fast32_t n_entries = brm_row_offset(mtx, mtx->base.rows);
    _Complex float* const elements = allocator_callbacks->alloc(allocator_callbacks->state, (n_entries) * sizeof (*elements));
    if (!elements)
    {
        allocator_callbacks->free(allocator_callbacks->state, out);
        return JMTX_RESULT_BAD_ALLOC;
    }

    memcpy(elements, mtx->values, sizeof* elements * n_entries);

    out->base.type = JMTXC_TYPE_BRM;
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

jmtx_result jmtxcs_matrix_brm_copy(const jmtxc_matrix_brm* mtx, jmtxc_matrix_brm** p_out, const jmtx_allocator_callbacks* allocator_callbacks)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTXC_TYPE_BRM)
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
    return jmtxc_matrix_brm_copy(mtx, p_out, allocator_callbacks);
}

_Complex float jmtxc_matrix_brm_vector_multiply_row(const jmtxc_matrix_brm* mtx, const _Complex float* x, uint32_t i)
{
    const uint_fast32_t offset = brm_row_offset(mtx, i);
    const uint_fast32_t first_col = jmtxc_matrix_brm_first_pos_in_row(mtx, i);
    const uint_fast32_t len = brm_row_len(mtx, i);
    _Complex float v = 0;
    const _Complex float* const values = mtx->values + offset;
    for (uint32_t j = 0; j < len; ++j)
    {
        v += values[j] * x[first_col + j];
    }
    return v;
}

jmtx_result jmtxcs_matrix_brm_vector_multiply_row(const jmtxc_matrix_brm* mtx, const _Complex float* restrict x, uint32_t i, _Complex float* restrict p_r)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTXC_TYPE_BRM)
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

    *p_r = jmtxc_matrix_brm_vector_multiply_row(mtx, x, i);
    return JMTX_RESULT_SUCCESS;
}

void jmtxc_matrix_brm_add_to_entry(const jmtxc_matrix_brm* mtx, uint32_t i, uint32_t j, _Complex float value)
{
    uint_fast32_t row_offset = brm_row_offset(mtx, i);
    uint_fast32_t col_offset = j - jmtxc_matrix_brm_first_pos_in_row(mtx, i);

    mtx->values[row_offset + col_offset] += value;
}

jmtx_result jmtxcs_matrix_brm_add_to_entry(const jmtxc_matrix_brm* mtx, uint32_t i, uint32_t j, _Complex float value)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTXC_TYPE_BRM)
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
    jmtxc_matrix_brm_add_to_entry(mtx, i, j, value);
    return JMTX_RESULT_SUCCESS;
}

void jmtxc_matrix_brm_zero_all_entries(const jmtxc_matrix_brm* mtx)
{
    memset(mtx->values, 0, sizeof(*mtx->values) * brm_row_offset(mtx, mtx->base.rows));
}

jmtx_result jmtxcs_matrix_brm_zero_all_entries(const jmtxc_matrix_brm* mtx)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTXC_TYPE_BRM)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    jmtxc_matrix_brm_zero_all_entries(mtx);
    return JMTX_RESULT_SUCCESS;
}

void jmtxc_matrix_brm_set_all_entries(const jmtxc_matrix_brm* mtx, _Complex float x)
{
    for (_Complex float* ptr = mtx->values; ptr != mtx->values + brm_row_offset(mtx, mtx->base.rows); ++ptr)
    {
        *ptr = x;
    }
}

JMTX_HOT_FUNCTION
uint_fast32_t jmtxc_matrix_brm_first_pos_in_row(const jmtxc_matrix_brm* mtx, uint32_t row)
{
    if (row < mtx->lower_bandwidth)
    {
        return 0;
    }
    return row - mtx->lower_bandwidth;
}

uint_fast32_t jmtxc_matrix_brm_last_pos_in_row(const jmtxc_matrix_brm* mtx, uint32_t row)
{
    if (row > mtx->base.rows - 1 - mtx->upper_bandwidth)
    {
        return mtx->base.rows - 1;
    }
    return row + mtx->upper_bandwidth;
}

uint_fast32_t jmtxc_matrix_brm_first_pos_in_col(const jmtxc_matrix_brm* mtx, uint32_t col)
{
    if (col < mtx->upper_bandwidth)
    {
        return 0;
    }
    return col - mtx->upper_bandwidth;
}

uint_fast32_t jmtxc_matrix_brm_last_pos_in_col(const jmtxc_matrix_brm* mtx, uint32_t col)
{
    if (col > mtx->base.rows - 1 - mtx->lower_bandwidth)
    {
        return mtx->base.cols - 1;
    }
    return col + mtx->lower_bandwidth;
}

uint_fast32_t jmtxc_matrix_brm_length_of_row(const jmtxc_matrix_brm* mtx, uint32_t row)
{
    return brm_row_len(mtx, row);
}

void jmtxc_matrix_brm_set_col(const jmtxc_matrix_brm* mtx, uint32_t col, const _Complex float values[])
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
    pos_rel = (col - jmtxc_matrix_brm_first_pos_in_row(mtx, first_row));
    i = brm_row_offset(mtx, first_row) + pos_rel;
    for (row = first_row, j = 0; row < last_row; ++row)
    {
        assert(i < max);
        mtx->values[i] = values[j++];
        uint_fast32_t new_pos_rel = (col - jmtxc_matrix_brm_first_pos_in_row(mtx, row + 1));
        i += (brm_row_len(mtx, row) - (pos_rel))
             +  new_pos_rel;
        pos_rel = new_pos_rel;
    }
}

void jmtxc_matrix_brm_get_bandwidths(const jmtxc_matrix_brm* mtx, uint32_t* ubw, uint32_t* lbw)
{
    *ubw = mtx->upper_bandwidth;
    *lbw = mtx->lower_bandwidth;
}

jmtx_result jmtxcs_matrix_brm_set_all_entries(const jmtxc_matrix_brm* mtx, _Complex float x)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTXC_TYPE_BRM)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    jmtxc_matrix_brm_set_all_entries(mtx, x);

    return JMTX_RESULT_SUCCESS;
}

jmtx_result jmtxcs_matrix_brm_set_col(const jmtxc_matrix_brm* mtx, uint32_t col, const _Complex float values[])
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTXC_TYPE_BRM)
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

    const uint_fast32_t len = jmtxc_matrix_brm_length_of_col(mtx, col);
    for (uint32_t i = 0; i < len; ++i)
    {
        if (!isfinite(crealf(values[i])) || !isfinite(cimagf(values[i])))
        {
            return JMTX_RESULT_BAD_PARAM;
        }
    }

    jmtxc_matrix_brm_set_col(mtx, col, values);
    return JMTX_RESULT_SUCCESS;
}

