// Automatically generated from source/float/matrices/sparse_column_compressed.c on Fri Dec  1 17:36:03 2023
//
// Created by jan on 15.6.2022.
//

#include <assert.h>
#include <math.h>
#include "../../../include/jmtx/cfloat/matrices/sparse_column_compressed.h"
#include "sparse_column_compressed_internal.h"
#include <complex.h>
//
// Created by jan on 13.6.2022.
//



enum {DEFAULT_RESERVED_ELEMENTS = 64};


static uint32_t ccs_get_column_entries(const jmtxc_matrix_ccs* mtx, uint32_t col, uint32_t** pp_indices, _Complex float** pp_values)
{
    uint32_t offset, len;
    if (col == 0)
    {
        offset = 0;
        len = mtx->end_of_column_offsets[0];
    }
    else
    {
        offset = mtx->end_of_column_offsets[col - 1];
        len = mtx->end_of_column_offsets[col] - offset;
    }
    *pp_indices = mtx->indices + offset;
    *pp_values = mtx->values + offset;

    return len;
}

/**
 * Inserts an entry (value-index pair) into the matrix at a specified column and a global position.
 * @param mtx matrix to which to add the entry
 * @param col what column the element is going to be in
 * @param position what is the position in the row
 * @param value what is the value of the entry
 * @param index the index of the entry
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC if memory allocation failed
 */
static jmtx_result ccs_insert_entry_at(jmtxc_matrix_ccs* mtx, uint32_t col, uint32_t position, _Complex float value, uint32_t index)
{
    const uint32_t global_position = position + (col ? mtx->end_of_column_offsets[col - 1] : 0);
    assert(position <= mtx->n_entries);
    assert(!col || (mtx->end_of_column_offsets[col - 1] == global_position || mtx->indices[global_position - 1] < index));

    if (mtx->capacity == mtx->n_entries)
    {
        //  Reallocate arrays
        const size_t new_capacity = mtx->capacity + DEFAULT_RESERVED_ELEMENTS;
        _Complex float* const new_values = mtx->base.allocator_callbacks.realloc(mtx->base.allocator_callbacks.state, mtx->values, sizeof(*mtx->values) * new_capacity);
        if (!new_values)
        {
            return JMTX_RESULT_BAD_ALLOC;
        }
        mtx->values = new_values;
        uint32_t* const new_indices = mtx->base.allocator_callbacks.realloc(mtx->base.allocator_callbacks.state, mtx->indices, sizeof(*mtx->indices) * new_capacity);
        if (!new_indices)
        {
            return JMTX_RESULT_BAD_ALLOC;
        }
        mtx->indices = new_indices;

        mtx->capacity = new_capacity;
    }

    //  Check for numer of values after the position
    const uint32_t elements_after = mtx->n_entries - global_position;
    if (elements_after)
    {
        //  Move other values out of the way
        memmove(mtx->values + global_position + 1, mtx->values + global_position, sizeof(*mtx->values) * (mtx->n_entries - global_position));
        memmove(mtx->indices + global_position + 1, mtx->indices + global_position, sizeof(*mtx->indices) * (mtx->n_entries - global_position));
    }
    //  Insert the values
    mtx->values[global_position] = value;
    mtx->indices[global_position] = index;

    //  Update offsets
    for (uint32_t i = col; i < mtx->base.cols; ++i)
    {
        mtx->end_of_column_offsets[i] += 1;
    }

    mtx->n_entries += 1;

    return JMTX_RESULT_SUCCESS;
}



jmtx_result jmtxc_matrix_ccs_new(
        jmtxc_matrix_ccs** p_mtx, uint32_t cols, uint32_t rows, uint32_t reserved_entries,
        const jmtx_allocator_callbacks* allocator_callbacks)
{
    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }

    if (reserved_entries == 0)
    {
        reserved_entries = DEFAULT_RESERVED_ELEMENTS;
        reserved_entries = reserved_entries < cols * rows ? reserved_entries : cols * rows;
    }

    uint32_t* offsets = NULL;
    uint32_t* indices = NULL;
    _Complex float* values = allocator_callbacks->alloc(allocator_callbacks->state, (reserved_entries) * sizeof(*values));
    if (!values)
    {
        return JMTX_RESULT_BAD_ALLOC;
    }
    memset(values, 0, (reserved_entries) * sizeof(*values));

    indices = allocator_callbacks->alloc(allocator_callbacks->state, (reserved_entries) * sizeof*indices);
    if (!indices)
    {
        allocator_callbacks->free(allocator_callbacks->state, values);
        return JMTX_RESULT_BAD_ALLOC;
    }
    memset(indices, 0, (reserved_entries) * sizeof(*indices));

    offsets = allocator_callbacks->alloc(allocator_callbacks->state, (cols) * sizeof*offsets);
    if (!offsets)
    {
        allocator_callbacks->free(allocator_callbacks->state, indices);
        allocator_callbacks->free(allocator_callbacks->state, values);
        return JMTX_RESULT_BAD_ALLOC;
    }
    memset(offsets, 0, (cols) * sizeof(*offsets));


    jmtxc_matrix_ccs* const this = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*this));
    if (!this)
    {
        allocator_callbacks->free(allocator_callbacks->state, offsets);
        allocator_callbacks->free(allocator_callbacks->state, indices);
        allocator_callbacks->free(allocator_callbacks->state, values);
        return JMTX_RESULT_BAD_ALLOC;
    }

    this->base.rows = rows;
    this->base.type = JMTXC_TYPE_CCS;
    this->base.cols = cols;
    this->base.allocator_callbacks = *allocator_callbacks;
    this->indices = indices;
    this->values = values;
    this->capacity = reserved_entries;
    this->n_entries = 0;
    this->end_of_column_offsets = offsets;

    *p_mtx = this;
    return JMTX_RESULT_SUCCESS;
}

jmtx_result jmtxcs_matrix_ccs_new(
        jmtxc_matrix_ccs** p_mtx, uint32_t cols, uint32_t rows, uint32_t reserved_entries,
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
    if (reserved_entries > cols * rows)
    {
        return JMTX_RESULT_BAD_PARAM;
    }
    if (allocator_callbacks && (!allocator_callbacks->alloc || !allocator_callbacks->realloc || !allocator_callbacks->free))
    {
        return JMTX_RESULT_BAD_PARAM;
    }
    return jmtxc_matrix_ccs_new(p_mtx, cols, rows, reserved_entries, allocator_callbacks);
}

void jmtxc_matrix_ccs_destroy(jmtxc_matrix_ccs* mtx)
{
    mtx->base.allocator_callbacks.free(mtx->base.allocator_callbacks.state, mtx->indices);
    mtx->base.allocator_callbacks.free(mtx->base.allocator_callbacks.state, mtx->end_of_column_offsets);
    mtx->base.allocator_callbacks.free(mtx->base.allocator_callbacks.state, mtx->values);
    mtx->base.allocator_callbacks.free(mtx->base.allocator_callbacks.state, mtx);
}

jmtx_result jmtxcs_matrix_ccs_destroy(jmtxc_matrix_ccs* mtx)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTXC_TYPE_CCS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    jmtxc_matrix_ccs_destroy(mtx);
    return JMTX_RESULT_SUCCESS;
}

jmtx_result jmtxc_matrix_ccs_shrink(jmtxc_matrix_ccs* mtx)
{
    if (mtx->n_entries == mtx->capacity)
    {
        return JMTX_RESULT_SUCCESS;
    }
    _Complex float* element_new_ptr = mtx->base.allocator_callbacks.realloc(mtx->base.allocator_callbacks.state, mtx->values, sizeof*mtx->values * (mtx->n_entries));
    if (!element_new_ptr)
    {
        return JMTX_RESULT_BAD_ALLOC;
    }
    mtx->values = element_new_ptr;
    uint32_t* new_indices_ptr = mtx->base.allocator_callbacks.realloc(mtx->base.allocator_callbacks.state, mtx->indices, sizeof*mtx->indices * (mtx->n_entries));
    if (!new_indices_ptr)
    {
        return JMTX_RESULT_BAD_ALLOC;
    }
    mtx->indices = new_indices_ptr;
    mtx->capacity = mtx->n_entries;

    return JMTX_RESULT_SUCCESS;
}

jmtx_result jmtxcs_matrix_ccs_shrink(jmtxc_matrix_ccs* mtx)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTXC_TYPE_CCS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    return jmtxc_matrix_ccs_shrink(mtx);
}

jmtx_result jmtxc_matrix_ccs_set_col(jmtxc_matrix_ccs* mtx, uint32_t col, uint32_t n, const uint32_t* indices, const _Complex float* values)
{
    const uint32_t beginning_offset = col ? mtx->end_of_column_offsets[col - 1] : 0;
    const int32_t new_elements = (int32_t)n - (int32_t)(mtx->end_of_column_offsets[col] - beginning_offset);
    const uint32_t required_capacity = (uint32_t)((int32_t)mtx->n_entries + new_elements);
    if (mtx->capacity < required_capacity)
    {
        _Complex float* new_element_ptr = mtx->base.allocator_callbacks.realloc(mtx->base.allocator_callbacks.state, mtx->values, sizeof*(mtx->values) * (required_capacity + 1));
        if (!new_element_ptr)
        {
            return JMTX_RESULT_BAD_ALLOC;
        }
        mtx->values = new_element_ptr;
        uint32_t* new_indices_ptr = mtx->base.allocator_callbacks.realloc(mtx->base.allocator_callbacks.state, mtx->indices, sizeof*(mtx->indices) * (required_capacity + 1));
        if (!new_indices_ptr)
        {
            return JMTX_RESULT_BAD_ALLOC;
        }
        mtx->indices = new_indices_ptr;
        mtx->capacity = required_capacity;
    }

    if (new_elements != 0)
    {
        const uint32_t elements_after = mtx->n_entries - mtx->end_of_column_offsets[col];
        if (elements_after)
        {
            memmove(mtx->values + mtx->end_of_column_offsets[col] + new_elements, mtx->values + mtx->end_of_column_offsets[col],
                    sizeof*mtx->values * (elements_after));
            memmove(mtx->indices + mtx->end_of_column_offsets[col] + new_elements, mtx->indices + mtx->end_of_column_offsets[col],
                    sizeof*mtx->indices * (elements_after));
        }
        memcpy(mtx->values + beginning_offset, values, sizeof*values * n);
        memcpy(mtx->indices + beginning_offset, indices, sizeof*indices * n);

        for (uint32_t i = col; i < mtx->base.cols; ++i)
        {
            mtx->end_of_column_offsets[i] += new_elements;
        }
        mtx->n_entries += new_elements;
    }
    else
    {
        memcpy(mtx->values + beginning_offset, values, sizeof*values * n);
        memcpy(mtx->indices + beginning_offset, indices, sizeof*indices * n);
    }
    return JMTX_RESULT_SUCCESS;
}

jmtx_result jmtxcs_matrix_ccs_set_col(jmtxc_matrix_ccs* mtx, uint32_t col, uint32_t n, const uint32_t* indices, const _Complex float* values)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTXC_TYPE_CCS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (mtx->base.cols <= col)
    {
        return JMTX_RESULT_INDEX_OUT_OF_BOUNDS;
    }
    if (mtx->base.rows < n)
    {
        return JMTX_RESULT_INDEX_OUT_OF_BOUNDS;
    }
    if (!indices)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!values)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (n != 0)
    {
        if (indices[0] >= mtx->base.rows)
        {
            return JMTX_RESULT_INDEX_OUT_OF_BOUNDS;
        }
        if (!isfinite(crealf(values[0])) || !isfinite(cimagf(values[0])))
        {
            return JMTX_RESULT_BAD_PARAM;
        }
        for (uint32_t i = 1; i < n; ++i)
        {
            if (indices[i - 1] >= indices[i] || (!isfinite(crealf(values[i])) || !isfinite(cimagf(values[i]))))
            {
                return JMTX_RESULT_BAD_PARAM;
            }
        }
    }
    return jmtxc_matrix_ccs_set_col(mtx, col, n, indices, values);
}

void jmtxc_matrix_ccs_vector_multiply(const jmtxc_matrix_ccs* mtx, const _Complex float* restrict x, _Complex float* restrict y)
{
    for (uint32_t i = 0; i < mtx->base.cols; ++i)
    {
        uint32_t* indices;
        _Complex float* values;
        const uint32_t n_elements = ccs_get_column_entries(mtx, i, &indices, &values);
        _Complex float v = 0;
        for (uint32_t j = 0; j < n_elements; ++j)
        {
            const uint32_t k = indices[j];
            v += values[j] * x[k];
        }
        y[i] = v;
    }
}

jmtx_result jmtxcs_matrix_ccs_vector_multiply(const jmtxc_matrix_ccs* mtx, const _Complex float* restrict x, _Complex float* restrict y)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTXC_TYPE_CCS)
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

    jmtxc_matrix_ccs_vector_multiply(mtx, x, y);
    return JMTX_RESULT_SUCCESS;
}

jmtx_result jmtxc_matrix_ccs_set_entry(jmtxc_matrix_ccs* mtx, uint32_t i, uint32_t j, _Complex float value)
{
    jmtx_result res;
    uint32_t* col_indices;
    _Complex float* col_values;
    const uint32_t n_col_elements = ccs_get_column_entries(mtx, j, &col_indices, &col_values);
    //  Check if row has any values
    if (n_col_elements != 0)
    {
        //  Find first column entry less or equal to it
        const uint32_t possible = jmtx_internal_find_last_leq_value(n_col_elements, col_indices, i);
        if (col_indices[possible] == i)
        {
            col_values[possible] = value;
            return JMTX_RESULT_SUCCESS;
        }
        else
        {
            //  Figure out where to insert
            if (col_indices[possible] <= i)
            {
                //  Insert right after the *possible*
                res = ccs_insert_entry_at(mtx, j, possible + 1, value, i);
            }
            else
            {
                //  Insert at the position of *possible*
                res = ccs_insert_entry_at(mtx, j, possible, value, i);
            }
        }
    }
    else
    {
        //  Insert the new element in an empty row
        res = ccs_insert_entry_at(mtx, j, 0, value, i);
    }

    return res;
}

jmtx_result jmtxcs_matrix_ccs_set_entry(jmtxc_matrix_ccs* mtx, uint32_t i, uint32_t j, _Complex float value)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTXC_TYPE_CCS)
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
    return jmtxc_matrix_ccs_set_entry(mtx, i, j, value);
}

_Complex float jmtxc_matrix_ccs_get_entry(const jmtxc_matrix_ccs* mtx, uint32_t i, uint32_t j)
{

    uint32_t* col_indices;
    _Complex float* col_values;
    const uint32_t n_col_elements = ccs_get_column_entries(mtx, j, &col_indices, &col_values);
    //  Check if row has any values
    if (n_col_elements != 0)
    {
        //  Find first column entry less or equal to it
        const uint32_t possible = jmtx_internal_find_last_leq_value(n_col_elements, col_indices, i);
        if (col_indices[possible] == i)
        {
            return  col_values[possible];
        }
    }
    return 0.0f;
}

jmtx_result jmtxcs_matrix_ccs_get_entry(const jmtxc_matrix_ccs* mtx, uint32_t i, uint32_t j, _Complex float* p_value)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTXC_TYPE_CCS)
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
    *p_value = jmtxc_matrix_ccs_get_entry(mtx, i, j);
    return JMTX_RESULT_SUCCESS;
}

uint32_t jmtxc_matrix_ccs_get_col(const jmtxc_matrix_ccs* mtx, uint32_t col, uint32_t** p_indices, _Complex float** p_elements)
{
    return ccs_get_column_entries(mtx, col, p_indices, p_elements);
}

jmtx_result jmtxcs_matrix_ccs_get_col(const jmtxc_matrix_ccs* mtx, uint32_t col, uint32_t* n, uint32_t** p_indices, _Complex float** p_elements)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTXC_TYPE_CCS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (col >= mtx->base.cols)
    {
        return JMTX_RESULT_INDEX_OUT_OF_BOUNDS;
    }
    if (!n)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!p_indices)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!p_elements)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    *n = jmtxc_matrix_ccs_get_col(mtx, col, p_indices, p_elements);
    return JMTX_RESULT_SUCCESS;
}


uint32_t jmtxc_matrix_ccs_count_values(const jmtxc_matrix_ccs* mtx, _Complex float v)
{
    uint32_t r = 0;
    for (uint32_t i = 0; i < mtx->n_entries; ++i)
    {
        if (v == mtx->values[i])
        {
            r += 1;
        }
    }
    return r;
}

jmtx_result jmtxcs_matrix_ccs_count_values(const jmtxc_matrix_ccs* mtx, _Complex float v, uint32_t* p_count)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTXC_TYPE_CRS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (!p_count)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }

    *p_count = jmtxc_matrix_ccs_count_values(mtx, v);
    return JMTX_RESULT_SUCCESS;
}

uint32_t jmtxc_matrix_ccs_count_indices(const jmtxc_matrix_ccs* mtx, uint32_t v)
{
    uint32_t r = 0;
    for (uint32_t i = 0; i < mtx->n_entries; ++i)
    {
        if (v == mtx->indices[i])
        {
            r += 1;
        }
    }
    return r;
}

jmtx_result jmtxcs_matrix_ccs_count_indices(const jmtxc_matrix_ccs* mtx, uint32_t v, uint32_t* p_count)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTXC_TYPE_CCS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (!p_count)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }

    *p_count = jmtxc_matrix_ccs_count_indices(mtx, v);
    return JMTX_RESULT_SUCCESS;
}

jmtx_result jmtxc_matrix_ccs_apply_unary_fn(const jmtxc_matrix_ccs* mtx, int (*unary_fn)(uint32_t i, uint32_t j, _Complex float* p_element, void* param), void* param)
{
    for (uint32_t j = 0; j < mtx->base.cols; ++j)
    {
        _Complex float* p_elements;
        uint32_t* p_indices;
        const uint32_t n_in_row = ccs_get_column_entries(mtx, j, &p_indices, &p_elements);
        for (uint32_t i = 0; i < n_in_row; ++i)
        {
            if ((unary_fn(j, p_indices[i], p_elements + i, param)))
            {
                return JMTX_RESULT_UNARY_RETURN;
            }
        }
    }
    return JMTX_RESULT_SUCCESS;
}

jmtx_result jmtxcs_matrix_ccs_apply_unary_fn(const jmtxc_matrix_ccs* mtx, int (*unary_fn)(uint32_t i, uint32_t j, _Complex float* p_element, void* param), void* param)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTXC_TYPE_CCS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (!unary_fn)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    return jmtxc_matrix_ccs_apply_unary_fn(mtx, unary_fn, param);
}

void jmtxc_matrix_ccs_remove_zeros(jmtxc_matrix_ccs* mtx)
{
//  Update offsets
    uint32_t p, c = 0, r = 0;
    while (mtx->end_of_column_offsets[r] == 0)
    {
        r += 1;
    }
    for (p = 0; p < mtx->n_entries; ++p)
    {
        assert(r < mtx->base.cols);
        if (mtx->values[p] == 0)
        {
            c += 1;
        }
        while (r < mtx->base.cols && p + 1 == mtx->end_of_column_offsets[r])
        {
            assert(mtx->end_of_column_offsets[r] >= c);
            mtx->end_of_column_offsets[r] -= c;
            r += 1;
        }
    }
    assert(r == mtx->base.cols);
    assert(c <= mtx->n_entries);
#ifndef NDEBUG
    for (r = 1; r < mtx->base.cols; ++r)
    {
        assert(mtx->end_of_column_offsets[r - 1] <= mtx->end_of_column_offsets[r]);
    }
#endif

    //  Remove the actual entries
    uint32_t p0 = mtx->n_entries;
    while (p0 != 0)
    {
        //  Check if entry must be removed
        if (mtx->values[p0 - 1] == 0)
        {
            uint32_t p1 = p0;
            while (p0 != 0 &&
                   mtx->values[p0 - 1] == 0)
            {
                p0 -= 1;
            }
            memmove(mtx->values + p0, mtx->values + p1, sizeof(*mtx->values) * (mtx->n_entries - p1));
            memmove(mtx->indices + p0, mtx->indices + p1, sizeof(*mtx->indices) * (mtx->n_entries - p1));
        }
        else
        {
            p0 -= 1;
        }
    }

    mtx->n_entries -= c;
}

jmtx_result jmtxcs_matrix_ccs_remove_zeros(jmtxc_matrix_ccs* mtx)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTXC_TYPE_CCS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }

    jmtxc_matrix_ccs_remove_zeros(mtx);
    return JMTX_RESULT_SUCCESS;
}



void jmtxc_matrix_ccs_remove_bellow(jmtxc_matrix_ccs* mtx, float v)
{
//  Update offsets
    uint32_t p, c = 0, r = 0;
    while (mtx->end_of_column_offsets[r] == 0)
    {
        r += 1;
    }
    for (p = 0; p < mtx->n_entries; ++p)
    {
        assert(r < mtx->base.cols);
        if (cabs(mtx->values[p]) < (v))
        {
            c += 1;
        }
        while (r < mtx->base.cols && p + 1 == mtx->end_of_column_offsets[r])
        {
            assert(mtx->end_of_column_offsets[r] >= c);
            mtx->end_of_column_offsets[r] -= c;
            r += 1;
        }
    }
    assert(r == mtx->base.rows);
    assert(c <= mtx->n_entries);
#ifndef NDEBUG
    for (r = 1; r < mtx->base.cols; ++r)
    {
        assert(mtx->end_of_column_offsets[r - 1] <= mtx->end_of_column_offsets[r]);
    }
#endif

    //  Remove the actual entries
    uint32_t p0 = mtx->n_entries;
    while (p0 != 0)
    {
        //  Check if entry must be removed
        if (cabsf(mtx->values[p0 - 1]) < v)
        {
            uint32_t p1 = p0;
            while (p0 != 0 && cabsf(mtx->values[p0 - 1]) < v)
            {
                p0 -= 1;
            }
            memmove(mtx->values + p0, mtx->values + p1, sizeof(*mtx->values) * (mtx->n_entries - p1));
            memmove(mtx->indices + p0, mtx->indices + p1, sizeof(*mtx->indices) * (mtx->n_entries - p1));
        }

        p0 -= 1;
    }


    //  Beef
    mtx->n_entries -= c;
}

jmtx_result jmtxcs_matrix_ccs_remove_bellow(jmtxc_matrix_ccs* mtx, float v)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTXC_TYPE_CCS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (!isfinite(v) || v < 0.0)
    {
        return JMTX_RESULT_BAD_PARAM;
    }
    jmtxc_matrix_ccs_remove_bellow(mtx, v);
    return JMTX_RESULT_SUCCESS;
}

uint32_t jmtxc_matrix_ccs_elements_in_row(const jmtxc_matrix_ccs* mtx, uint32_t row)
{
    uint32_t element_count = 0;
    for (uint32_t col = 0; col < mtx->base.cols && (!col || (mtx->end_of_column_offsets[col - 1] != mtx->n_entries)); ++col)
    {
        uint32_t* col_indices;
        _Complex float* unused_col_values;
        const uint32_t n_col_elements = ccs_get_column_entries(mtx, col, &col_indices, &unused_col_values);
        if (n_col_elements && col_indices[0] <= row && col_indices[n_col_elements - 1] >= row)
        {
            uint32_t current = jmtx_internal_find_last_leq_value(n_col_elements, col_indices, row);
            if (col_indices[current] == row)
            {
                element_count += 1;
            }
        }
    }
    return element_count;
}

jmtx_result jmtxcs_matrix_ccs_elements_in_row(const jmtxc_matrix_ccs* mtx, uint32_t row, uint32_t* p_n)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTXC_TYPE_CCS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (row >= mtx->base.rows)
    {
        return JMTX_RESULT_INDEX_OUT_OF_BOUNDS;
    }
    if (!p_n)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    *p_n = jmtxc_matrix_ccs_elements_in_row(mtx, row);
    return JMTX_RESULT_SUCCESS;
}

uint32_t jmtxc_matrix_ccs_get_row(
    const jmtxc_matrix_ccs* mtx, uint32_t row, uint32_t n, _Complex float* p_values, uint32_t* p_columns)
{
    uint32_t k = 0;
    for (uint32_t col = 0; k < n && col < mtx->base.cols && (!col || (mtx->end_of_column_offsets[col - 1] != mtx->n_entries)); ++col)
    {
        uint32_t* col_indices;
        _Complex float* unused_col_values;
        const uint32_t n_col_elements = ccs_get_column_entries(mtx, col, &col_indices, &unused_col_values);
        if (n_col_elements && col_indices[0] <= row && col_indices[n_col_elements - 1] >= row)
        {
            const uint32_t current = jmtx_internal_find_last_leq_value(n_col_elements, col_indices, row);
            if (col_indices[current] == row)
            {
                p_values[k] = mtx->values[(col ? mtx->end_of_column_offsets[col - 1] : 0) + current];
                p_columns[k] = col;
                k += 1;
            }
        }
    }

    return k;
}

jmtx_result jmtxcs_matrix_ccs_get_row(
        const jmtxc_matrix_ccs* mtx, uint32_t row, uint32_t n, _Complex float* p_values, uint32_t* p_count, uint32_t* p_columns)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTXC_TYPE_CCS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (row >= mtx->base.rows)
    {
        return JMTX_RESULT_INDEX_OUT_OF_BOUNDS;
    }
    if (p_values == NULL)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!p_columns)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!p_count)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    *p_count = jmtxc_matrix_ccs_get_row(mtx, row, n, p_values, p_columns);
    return JMTX_RESULT_SUCCESS;
}

jmtx_result jmtxc_matrix_ccs_transpose(const jmtxc_matrix_ccs* mtx, jmtxc_matrix_ccs** p_out,
                                      const jmtx_allocator_callbacks* allocator_callbacks)
{
    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &mtx->base.allocator_callbacks;
    }

    const uint32_t n_elements = mtx->n_entries;
    const uint32_t new_cols = mtx->base.rows;
    const uint32_t new_rows = mtx->base.cols;

    jmtxc_matrix_ccs* const out = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*out));
    if (!out)
    {
        return JMTX_RESULT_BAD_ALLOC;
    }

    uint32_t* const row_cum_counts = mtx->base.allocator_callbacks.alloc(mtx->base.allocator_callbacks.state, (new_cols) * sizeof(*row_cum_counts));
    if (!row_cum_counts)
    {
        allocator_callbacks->free(allocator_callbacks->state, out);
        return JMTX_RESULT_BAD_ALLOC;
    }
    uint32_t* const new_indices = mtx->base.allocator_callbacks.alloc(mtx->base.allocator_callbacks.state, (n_elements) * sizeof(*new_indices));
    if (!new_indices)
    {
        mtx->base.allocator_callbacks.free(mtx->base.allocator_callbacks.state, row_cum_counts);
        allocator_callbacks->free(allocator_callbacks->state, out);
        return JMTX_RESULT_BAD_ALLOC;
    }
    memset(new_indices, 0, (n_elements) * sizeof*new_indices);
    _Complex float* const new_elements = mtx->base.allocator_callbacks.alloc(mtx->base.allocator_callbacks.state, (n_elements) * sizeof*new_elements);
    if (!new_elements)
    {
        mtx->base.allocator_callbacks.free(mtx->base.allocator_callbacks.state, row_cum_counts);
        mtx->base.allocator_callbacks.free(mtx->base.allocator_callbacks.state, new_indices);
        allocator_callbacks->free(allocator_callbacks->state, out);
        return JMTX_RESULT_BAD_ALLOC;
    }

    *row_cum_counts = 0;

    for (uint32_t j = 0, n, p = 0; j < mtx->base.rows; ++j)
    {
        n = jmtxc_matrix_ccs_elements_in_row(mtx, j);
        uint32_t c = jmtxc_matrix_ccs_get_row(mtx, j, n, new_elements + p, new_indices + p);
        assert(c == n);
        p += n;
        row_cum_counts[j] = (j > 0 ? row_cum_counts[j - 1] : 0) + n;
    }

    out->end_of_column_offsets = row_cum_counts;
    out->values = new_elements;
    out->indices = new_indices;
    out->capacity = n_elements;
    out->n_entries = n_elements;
    out->base = mtx->base;
    out->base.rows = new_rows;
    out->base.cols = new_cols;
    out->base.allocator_callbacks = *allocator_callbacks;
    *p_out = out;
    return JMTX_RESULT_SUCCESS;
}

jmtx_result jmtxcs_matrix_ccs_transpose(const jmtxc_matrix_ccs* mtx, jmtxc_matrix_ccs** p_out,
                                      const jmtx_allocator_callbacks* allocator_callbacks)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTXC_TYPE_CCS)
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
    return jmtxc_matrix_ccs_transpose(mtx, p_out, allocator_callbacks);
}

jmtx_result jmtxc_matrix_ccs_copy(const jmtxc_matrix_ccs* mtx, jmtxc_matrix_ccs** p_out, const jmtx_allocator_callbacks* allocator_callbacks)
{
    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &mtx->base.allocator_callbacks;
    }

    jmtxc_matrix_ccs* const this = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*this));
    if (!this)
    {
        return JMTX_RESULT_BAD_ALLOC;
    }

    _Complex float* const elements = mtx->base.allocator_callbacks.alloc(mtx->base.allocator_callbacks.state, (mtx->n_entries) * sizeof *elements);
    if (!elements)
    {
        allocator_callbacks->free(allocator_callbacks->state, this);
        return JMTX_RESULT_BAD_ALLOC;
    }
    uint32_t* const indices = mtx->base.allocator_callbacks.alloc(mtx->base.allocator_callbacks.state, (mtx->n_entries) * sizeof *indices);
    if (!indices)
    {
        mtx->base.allocator_callbacks.free(mtx->base.allocator_callbacks.state, elements);
        allocator_callbacks->free(allocator_callbacks->state, this);
        return JMTX_RESULT_BAD_ALLOC;
    }
    uint32_t* const cum_sum = mtx->base.allocator_callbacks.alloc(mtx->base.allocator_callbacks.state, (mtx->base.cols) * sizeof *cum_sum);
    if (!cum_sum)
    {
        mtx->base.allocator_callbacks.free(mtx->base.allocator_callbacks.state, indices);
        mtx->base.allocator_callbacks.free(mtx->base.allocator_callbacks.state, elements);
        allocator_callbacks->free(allocator_callbacks->state, this);
        return JMTX_RESULT_BAD_ALLOC;
    }

    memcpy(elements, mtx->values, sizeof* elements * mtx->n_entries);
    memcpy(indices, mtx->indices, sizeof* indices * mtx->n_entries);
    memcpy(cum_sum, mtx->end_of_column_offsets, sizeof* cum_sum * (mtx->base.cols));
    this->base = mtx->base;
    this->values = elements;
    this->indices = indices;
    this->end_of_column_offsets = cum_sum;
    this->capacity = mtx->n_entries;
    this->n_entries = this->n_entries;
    *p_out = this;
    return JMTX_RESULT_SUCCESS;
}

jmtx_result jmtxcs_matrix_ccs_copy(const jmtxc_matrix_ccs* mtx, jmtxc_matrix_ccs** p_out, const jmtx_allocator_callbacks* allocator_callbacks)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTXC_TYPE_CCS)
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

    return jmtxc_matrix_ccs_copy(mtx, p_out, allocator_callbacks);
}

void jmtxc_matrix_ccs_zero_all_elements(const jmtxc_matrix_ccs* mtx)
{
    memset(mtx->values, 0, sizeof(*mtx->values) * mtx->n_entries);
}

jmtx_result jmtxcs_matrix_ccs_zero_all_elements(const jmtxc_matrix_ccs* mtx)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTXC_TYPE_CCS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }

    jmtxc_matrix_ccs_zero_all_elements(mtx);

    return JMTX_RESULT_SUCCESS;
}

void jmtxc_matrix_ccs_set_all_elements(const jmtxc_matrix_ccs* mtx, _Complex float x)
{
    for (_Complex float* ptr = mtx->values; ptr != mtx->values + mtx->n_entries; ++ptr)
    {
        *ptr = x;
    }

}

jmtx_result jmtxcs_matrix_ccs_set_all_elements(jmtxc_matrix_ccs* mtx, _Complex float x)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTXC_TYPE_CCS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (!isfinite(crealf(x)) || !isfinite(cimagf(x)))
    {
        return JMTX_RESULT_BAD_PARAM;
    }
    jmtxc_matrix_ccs_set_all_elements(mtx, x);
    return JMTX_RESULT_SUCCESS;
}

jmtx_result
jmtxc_matrix_ccs_build_col(jmtxc_matrix_ccs* mtx, uint32_t col, uint32_t n, const uint32_t* indices, const _Complex float* values)
{
    const uint32_t required_capacity = (uint32_t)((int32_t)mtx->n_entries + (int32_t)n);
    if (mtx->capacity < required_capacity)
    {
        _Complex float* new_element_ptr = mtx->base.allocator_callbacks.realloc(mtx->base.allocator_callbacks.state, mtx->values, sizeof*mtx->values * (required_capacity + 1));
        if (!new_element_ptr)
        {
            return JMTX_RESULT_BAD_ALLOC;
        }
        mtx->values = new_element_ptr;
        uint32_t* new_indices_ptr = mtx->base.allocator_callbacks.realloc(mtx->base.allocator_callbacks.state, mtx->indices, sizeof*mtx->indices * (required_capacity + 1));
        if (!new_indices_ptr)
        {
            return JMTX_RESULT_BAD_ALLOC;
        }
        mtx->indices = new_indices_ptr;
        mtx->capacity = required_capacity;
    }

    const uint32_t offset = col ? mtx->end_of_column_offsets[col - 1] : 0;
    memcpy(mtx->values + offset, values, sizeof*values * n);
    memcpy(mtx->indices + offset, indices, sizeof*indices * n);

    mtx->end_of_column_offsets[col] = n + offset;
    mtx->n_entries += n;

    return JMTX_RESULT_SUCCESS;
}

/**
 * Finds the upper bandwidth of the matrix; what is the furthest distance of and entry above the main diagonal
 * @param mtx matrx to find the upper bandwidth of
 * @return upper bandwidth of the matrix
 */
uint32_t jmtxc_matrix_ccs_find_upper_bandwidth(const jmtxc_matrix_ccs* mtx)
{
    //  Find the greatest distance above the main diagonal
    uint_fast32_t v_max = 0;
    for (uint_fast32_t i = 0, p = 0; i < mtx->base.cols; ++i)
    {
        uint_fast32_t j;
        for (j = p; j < mtx->end_of_column_offsets[i]; ++j)
        {
            if (mtx->indices[j] < i)
            {
                const uint_fast32_t dif = i - mtx->indices[j];
                if (dif > v_max)
                {
                    v_max = dif;
                }
            }
        }
        p += j;
    }
    return v_max;
}

/**
 * Finds the lower bandwidth of the matrix; what is the furthest distance of and entry bellow the main diagonal
 * @param mtx matrx to find the lower bandwidth of
 * @return lower bandwidth of the matrix
 */
uint32_t jmtxc_matrix_ccs_find_lower_bandwidth(const jmtxc_matrix_ccs* mtx)
{
    //  Find the greatest distance above the main diagonal
    uint_fast32_t v_max = 0;
    for (uint_fast32_t i = 0, p = 0; i < mtx->base.cols; ++i)
    {
        uint_fast32_t j;
        for (j = p; j < mtx->end_of_column_offsets[i]; ++j)
        {
            if (mtx->indices[j] > i)
            {
                const uint_fast32_t dif = mtx->indices[j] - i;
                if (dif > v_max)
                {
                    v_max = dif;
                }
            }
        }
        p += j;
    }
    return v_max;
}

