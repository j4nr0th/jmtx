// Automatically generated from source/float/matrices/sparse_row_compressed.c on Fri Dec  1 06:43:01 2023
//
// Created by jan on 13.6.2022.
//

#include <assert.h>
#include <math.h>
#include "../../../include/jmtx/double/matrices/sparse_row_compressed.h"
#include "sparse_row_compressed_internal.h"

enum{DEFAULT_RESERVED_ELEMENTS = 64};

static uint32_t crs_get_row_entries(const jmtxd_matrix_crs* mtx, uint32_t row, uint32_t* pp_indices[1], double* pp_values[1])
{
    uint32_t offset, len;
    if (row == 0)
    {
        offset = 0;
        len = mtx->end_of_row_offsets[0];
    }
    else
    {
        offset = mtx->end_of_row_offsets[row - 1];
        len = mtx->end_of_row_offsets[row] - offset;
    }
    *pp_indices = mtx->indices + offset;
    *pp_values = mtx->values + offset;

    return len;
}

/**
 * Inserts an entry (value-index pair) into the matrix at a specified row and a global position.
 * @param mtx matrix to which to add the entry
 * @param row what row the element is going to be in
 * @param position what is the position in the row
 * @param value what is the value of the entry
 * @param index the index of the entry
 * @return JMTX_RESULT_BAD_ALLOC on allocation failure
 */
static jmtx_result crs_insert_entry_at(jmtxd_matrix_crs* mtx, uint32_t row, uint32_t position, double value, uint32_t index)
{
    const uint32_t global_position = position + (row ? mtx->end_of_row_offsets[row - 1] : 0);
    assert(position <= mtx->n_entries);
    assert(!row || (mtx->end_of_row_offsets[row - 1] == global_position || mtx->indices[global_position - 1] < index));

    if (mtx->capacity == mtx->n_entries)
    {
        //  Reallocate arrays
        const size_t new_capacity = mtx->capacity + DEFAULT_RESERVED_ELEMENTS;
        double* const new_values = mtx->base.allocator_callbacks.realloc(mtx->base.allocator_callbacks.state, mtx->values, sizeof(*mtx->values) * new_capacity);
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
        mtx->capacity = (uint32_t)new_capacity;
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
    for (uint32_t i = row; i < mtx->base.rows; ++i)
    {
        mtx->end_of_row_offsets[i] += 1;
    }

    mtx->n_entries += 1;

    return JMTX_RESULT_SUCCESS;
}

jmtx_result jmtxd_matrix_crs_new(
    jmtxd_matrix_crs** p_mtx, uint32_t rows, uint32_t cols, uint32_t reserved_entries,
    const jmtx_allocator_callbacks* allocator_callbacks)
{
    if (reserved_entries == 0)
    {
        reserved_entries = DEFAULT_RESERVED_ELEMENTS;
        reserved_entries = reserved_entries < ((uint64_t)cols * (uint64_t)rows) ? reserved_entries : ((uint64_t)cols * (uint64_t)rows);
    }
    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }

    ;
    uint32_t* offsets = NULL;
    uint32_t* indices = NULL;

    jmtxd_matrix_crs* mtx = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*mtx));
    if (!mtx)
    {
        return JMTX_RESULT_BAD_ALLOC;
    }

    double* values = allocator_callbacks->alloc(allocator_callbacks->state, (reserved_entries) * sizeof(*values));
    if (!values)
    {
        allocator_callbacks->free(allocator_callbacks->state, mtx);
        return JMTX_RESULT_BAD_ALLOC;
    }
    memset(values, 0, (reserved_entries) * sizeof(*values));

    indices = allocator_callbacks->alloc(allocator_callbacks->state, (reserved_entries) * sizeof(*indices));
    if (!indices)
    {
        allocator_callbacks->free(allocator_callbacks->state, indices);
        allocator_callbacks->free(allocator_callbacks->state, mtx);
        return JMTX_RESULT_BAD_ALLOC;
    }
    memset(indices, 0, (reserved_entries) * sizeof(*indices));

    offsets = allocator_callbacks->alloc(allocator_callbacks->state, (rows) * sizeof(*offsets));
    if (!offsets)
    {
        allocator_callbacks->free(allocator_callbacks->state, offsets);
        allocator_callbacks->free(allocator_callbacks->state, indices);
        allocator_callbacks->free(allocator_callbacks->state, mtx);
        return JMTX_RESULT_BAD_ALLOC;
    }
    memset(offsets, 0, (rows) * sizeof(*offsets));

    memset(mtx, 0, sizeof*mtx);
    mtx->base.cols = cols;
    mtx->base.type = JMTXD_TYPE_CRS;
    mtx->base.rows = rows;
    mtx->base.allocator_callbacks = *allocator_callbacks;
    mtx->indices = indices;
    mtx->values = values;
    mtx->capacity = reserved_entries;
    mtx->n_entries = 0;
    mtx->end_of_row_offsets = offsets;
    *p_mtx = mtx;

    return JMTX_RESULT_SUCCESS;
}

jmtx_result jmtxds_matrix_crs_new(
    jmtxd_matrix_crs** p_mtx, uint32_t rows, uint32_t cols, uint32_t reserved_entries,
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
    if (reserved_entries > ((uint64_t)cols * (uint64_t)rows))
    {
        return JMTX_RESULT_BAD_PARAM;
    }
    if (allocator_callbacks && (!allocator_callbacks->alloc || !allocator_callbacks->realloc || !allocator_callbacks->free))
    {
        return JMTX_RESULT_BAD_PARAM;
    }

    return jmtxd_matrix_crs_new(p_mtx, rows, cols, reserved_entries, allocator_callbacks);
}

void jmtxd_matrix_crs_destroy(jmtxd_matrix_crs* mtx)
{
    jmtx_allocator_callbacks allocator = mtx->base.allocator_callbacks;
    allocator.free(allocator.state, mtx->indices);
    allocator.free(allocator.state, mtx->end_of_row_offsets);
    allocator.free(allocator.state, mtx->values);
    allocator.free(allocator.state, mtx);
}

jmtx_result jmtxds_matrix_crs_destroy(jmtxd_matrix_crs* mtx)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTXD_TYPE_CRS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    jmtxd_matrix_crs_destroy(mtx);
    return JMTX_RESULT_SUCCESS;
}

jmtx_result jmtxd_matrix_crs_shrink(jmtxd_matrix_crs* mtx)
{
    if (mtx->n_entries == mtx->capacity)
    {
        return JMTX_RESULT_SUCCESS;
    }

    double* element_new_ptr = mtx->base.allocator_callbacks.realloc(mtx->base.allocator_callbacks.state, mtx->values, sizeof*mtx->values * (mtx->n_entries));
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

jmtx_result jmtxds_matrix_crs_shrink(jmtxd_matrix_crs* mtx)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTXD_TYPE_CRS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    return jmtxd_matrix_crs_shrink(mtx);
}

jmtx_result jmtxd_matrix_crs_set_row(jmtxd_matrix_crs* mtx, uint32_t row, uint32_t n, const uint32_t indices[JMTX_ARRAY_ATTRIB(static n)], const double values[JMTX_ARRAY_ATTRIB(static n)])
{

    jmtx_result res = JMTX_RESULT_SUCCESS;
    const uint32_t beginning_offset = row ? mtx->end_of_row_offsets[row - 1] : 0;
    const int32_t new_elements = (int32_t)n - (int32_t)(mtx->end_of_row_offsets[row] - beginning_offset);
    const uint32_t required_capacity = (uint32_t)((int32_t)mtx->n_entries + new_elements);
    if (mtx->capacity < required_capacity)
    {
        double* new_element_ptr = mtx->base.allocator_callbacks.realloc(mtx->base.allocator_callbacks.state, mtx->values, sizeof*(mtx->values) * (required_capacity + 1));
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
        const uint32_t elements_after = mtx->n_entries - mtx->end_of_row_offsets[row];
        if (elements_after)
        {
            memmove(mtx->values + mtx->end_of_row_offsets[row] + new_elements, mtx->values + mtx->end_of_row_offsets[row],
                    sizeof*mtx->values * (elements_after));
            memmove(mtx->indices + mtx->end_of_row_offsets[row] + new_elements, mtx->indices + mtx->end_of_row_offsets[row],
                sizeof*mtx->indices * (elements_after));
        }
        memcpy(mtx->values + beginning_offset, values, sizeof*values * n);
        memcpy(mtx->indices + beginning_offset, indices, sizeof*indices * n);

        for (uint32_t i = row; i < mtx->base.rows; ++i)
        {
            mtx->end_of_row_offsets[i] += new_elements;
        }
        mtx->n_entries += new_elements;
    }
    else
    {
        memcpy(mtx->values + beginning_offset, values, sizeof*values * n);
        memcpy(mtx->indices + beginning_offset, indices, sizeof*indices * n);
    }

    return res;
}

jmtx_result jmtxds_matrix_crs_set_row(jmtxd_matrix_crs* mtx, uint32_t row, uint32_t n, const uint32_t* indices, const double* values)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTXD_TYPE_CRS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (mtx->base.rows <= row)
    {
        return JMTX_RESULT_INDEX_OUT_OF_BOUNDS;
    }
    if (mtx->base.cols < n)
    {
        return JMTX_RESULT_BAD_PARAM;
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
        if (indices[0] >= mtx->base.cols)
        {
            return JMTX_RESULT_INDEX_OUT_OF_BOUNDS;
        }
        if (!isfinite(values[0]))
        {
            return JMTX_RESULT_BAD_PARAM;
        }
        for (uint32_t i = 1; i < n; ++i)
        {
            if (indices[i - 1] >= indices[i] || !isfinite(values[i]))
            {
                return JMTX_RESULT_BAD_PARAM;
            }
        }
    }
    return jmtxd_matrix_crs_set_row(mtx, row, n, indices, values);
}


void jmtxd_matrix_crs_vector_multiply(const jmtxd_matrix_crs* mtx, const double* restrict x, double* restrict y)
{
    for (uint32_t i = 0; i < mtx->base.rows; ++i)
    {
        uint32_t* indices;
        double* values;
        const uint32_t row_entries = crs_get_row_entries(mtx, i, &indices, &values);
        double v = 0;
        for (uint32_t j = 0; j < row_entries; ++j)
        {
            const uint32_t k = indices[j];
            v += values[j] * x[k];
        }
        y[i] = v;
    }
}

jmtx_result jmtxds_matrix_crs_vector_multiply(const jmtxd_matrix_crs* mtx, const double* restrict x, double* restrict y)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTXD_TYPE_CRS)
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
    jmtxd_matrix_crs_vector_multiply(mtx, x, y);
    return JMTX_RESULT_SUCCESS;
}

jmtx_result jmtxd_matrix_crs_set_entry(jmtxd_matrix_crs* mtx, uint32_t i, uint32_t j, double value)
{
    jmtx_result res;
    uint32_t* row_indices;
    double* row_values;
    const uint32_t n_row_elements = crs_get_row_entries(mtx, i, &row_indices, &row_values);
    //  Check if row has any values
    if (n_row_elements != 0)
    {
        //  Find first column entry less or equal to it
        const uint32_t possible = jmtx_internal_find_last_leq_value(n_row_elements, row_indices, j);
        if (row_indices[possible] == j)
        {
            row_values[possible] = value;
            return JMTX_RESULT_SUCCESS;
        }
        else
        {
            //  Figure out where to insert
            if (row_indices[possible] <= j)
            {
                //  Insert right after the *possible*
                res = crs_insert_entry_at(mtx, i, possible + 1, value, j);
            }
            else
            {
                //  Insert at the position of *possible*
                res = crs_insert_entry_at(mtx, i, possible, value, j);
            }
        }
    }
    else
    {
        //  Insert the new element in an empty row
        res = crs_insert_entry_at(mtx, i, 0, value, j);
    }

    return res;
}

jmtx_result jmtxds_matrix_crs_set_entry(jmtxd_matrix_crs* mtx, uint32_t i, uint32_t j, double value)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTXD_TYPE_CRS)
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
    return jmtxd_matrix_crs_set_entry(mtx, i, j, value);
}

double jmtxd_matrix_crs_get_entry(const jmtxd_matrix_crs* mtx, uint32_t i, uint32_t j)
{
    uint32_t* row_indices;
    double* row_values;
    const uint32_t n_row_elements = crs_get_row_entries(mtx, i, &row_indices, &row_values);
    //  Check if row has any values
    if (n_row_elements != 0)
    {
        //  Find first column entry less or equal to it
        const uint32_t possible = jmtx_internal_find_last_leq_value(n_row_elements, row_indices, j);
        if (row_indices[possible] == j)
        {
            return row_values[possible];
        }
    }

    return 0.0f;
}

jmtx_result jmtxds_matrix_crs_get_entry(const jmtxd_matrix_crs* mtx, uint32_t i, uint32_t j, double* p_value)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTXD_TYPE_CRS)
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
    *p_value = jmtxd_matrix_crs_get_entry(mtx, i, j);

    return JMTX_RESULT_SUCCESS;
}

uint32_t jmtxd_matrix_crs_get_row(const jmtxd_matrix_crs* mtx, uint32_t row, uint32_t* p_indices[1], double* p_elements[1])
{
    return crs_get_row_entries(mtx, row, p_indices, p_elements);
}

jmtx_result jmtxds_matrix_crs_get_row(const jmtxd_matrix_crs* mtx, uint32_t row, uint32_t* n, uint32_t** p_indices, double** p_elements)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTXD_TYPE_CRS)
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
    if (!p_indices)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!p_elements)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    *n = jmtxd_matrix_crs_get_row(mtx, row, p_indices, p_elements);
    return JMTX_RESULT_SUCCESS;
}

// Maybe remove this one. It lacks and all purpose
uint32_t jmtxd_matrix_crs_count_values(const jmtxd_matrix_crs* mtx, double v)
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

jmtx_result jmtxds_matrix_crs_count_values(const jmtxd_matrix_crs* mtx, double v, uint32_t* p_count)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTXD_TYPE_CRS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (!p_count)
    {
        return JMTX_RESULT_NULL_PARAM;
    }

    *p_count = jmtxd_matrix_crs_count_values(mtx, v);
    return JMTX_RESULT_SUCCESS;
}

uint32_t jmtxd_matrix_crs_count_indices(const jmtxd_matrix_crs* mtx, uint32_t v)
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

jmtx_result jmtxds_matrix_crs_count_indices(const jmtxd_matrix_crs* mtx, uint32_t v, uint32_t* p_count)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTXD_TYPE_CRS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (!p_count)
    {
        return JMTX_RESULT_NULL_PARAM;
    }

    *p_count = jmtxd_matrix_crs_count_indices(mtx, v);
    return JMTX_RESULT_SUCCESS;
}


jmtx_result jmtxd_matrix_crs_apply_unary_fn(const jmtxd_matrix_crs* mtx, int (*unary_fn)(uint32_t i, uint32_t j, double* p_value, void* param), void* param)
{
    for (uint32_t i = 0; i < mtx->base.rows; ++i)
    {
        double* p_elements;
        uint32_t* p_indices;
        const uint32_t n_in_row = crs_get_row_entries(mtx, i, &p_indices, &p_elements);
        for (uint32_t j = 0; j < n_in_row; ++j)
        {
            if ((unary_fn(i, p_indices[j], p_elements + j, param)))
            {

                return JMTX_RESULT_UNARY_RETURN;
            }
        }
    }

    return JMTX_RESULT_SUCCESS;
}

jmtx_result jmtxds_matrix_crs_apply_unary_fn(const jmtxd_matrix_crs* mtx, int (*unary_fn)(uint32_t i, uint32_t j, double* p_value, void* param), void* param)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTXD_TYPE_CRS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (!unary_fn)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    return jmtxd_matrix_crs_apply_unary_fn(mtx, unary_fn, param);
}

void jmtxd_matrix_crs_remove_zeros(jmtxd_matrix_crs* mtx)
{
    //  Update offsets
    uint32_t p, c = 0, r = 0;
    while (mtx->end_of_row_offsets[r] == 0)
    {
        r += 1;
    }
    for (p = 0; p < mtx->n_entries; ++p)
    {
        assert(r < mtx->base.rows);
        if (mtx->values[p] == 0)
        {
            c += 1;
        }
        while (r < mtx->base.cols && p + 1 == mtx->end_of_row_offsets[r])
        {
            assert(mtx->end_of_row_offsets[r] >= c);
            mtx->end_of_row_offsets[r] -= c;
            r += 1;
        }
    }
    assert(r == mtx->base.rows);
    assert(c <= mtx->n_entries);
#ifndef NDEBUG
    for (r = 1; r < mtx->base.rows; ++r)
    {
        assert(mtx->end_of_row_offsets[r - 1] <= mtx->end_of_row_offsets[r]);
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

jmtx_result jmtxds_matrix_crs_remove_zeros(jmtxd_matrix_crs* mtx)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTXD_TYPE_CRS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }

    jmtxd_matrix_crs_remove_zeros(mtx);
    return JMTX_RESULT_SUCCESS;
}


void jmtxd_matrix_crs_remove_bellow_magnitude(jmtxd_matrix_crs* mtx, double v)
{
    //  Update offsets
    uint32_t p, c = 0, r = 0;
    while (mtx->end_of_row_offsets[r] == 0)
    {
        r += 1;
    }
    for (p = 0; p < mtx->n_entries; ++p)
    {
        assert(r < mtx->base.rows);
        if (fabs(mtx->values[p]) < v)
        {
            c += 1;
        }
        while (r < mtx->base.cols && p + 1 == mtx->end_of_row_offsets[r])
        {
            assert(mtx->end_of_row_offsets[r] >= c);
            mtx->end_of_row_offsets[r] -= c;
            r += 1;
        }
    }
    assert(r == mtx->base.rows);
    assert(c <= mtx->n_entries);
#ifndef NDEBUG
    for (r = 1; r < mtx->base.rows; ++r)
    {
        assert(mtx->end_of_row_offsets[r - 1] <= mtx->end_of_row_offsets[r]);
    }
#endif

    //  Remove the actual entries
    uint32_t p0 = mtx->n_entries;
    while (p0 != 0)
    {
        //  Check if entry must be removed
        if (fabs(mtx->values[p0 - 1]) < v)
        {
            uint32_t p1 = p0;
            while (p0 != 0 && fabs(mtx->values[p0 - 1]) < v)
            {
                p0 -= 1;
            }
            memmove(mtx->values + p0, mtx->values + p1, sizeof(*mtx->values) * (mtx->n_entries - p1));
            memmove(mtx->indices + p0, mtx->indices + p1, sizeof(*mtx->indices) * (mtx->n_entries - p1));
        }

        p0 -= 1;
    }

    mtx->n_entries -= c;
}

jmtx_result jmtxds_matrix_crs_remove_bellow_magnitude(jmtxd_matrix_crs* mtx, double v)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTXD_TYPE_CRS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (!isfinite(v))
    {
        return JMTX_RESULT_BAD_PARAM;
    }

    jmtxd_matrix_crs_remove_bellow_magnitude(mtx, v);
    return JMTX_RESULT_SUCCESS;
}

uint32_t jmtxd_matrix_crs_entries_in_col(const jmtxd_matrix_crs* mtx, uint32_t col)
{
    uint32_t element_count = 0;
    for (uint32_t row = 0; row < mtx->base.rows && (!row || (mtx->end_of_row_offsets[row - 1] != mtx->n_entries)); ++row)
    {
        uint32_t* row_indices;
        double* unused_row_values;
        const uint32_t n_row_elements = crs_get_row_entries(mtx, row, &row_indices, &unused_row_values);
        if (n_row_elements && row_indices[0] <= col && row_indices[n_row_elements - 1] >= col)
        {
            uint32_t current = jmtx_internal_find_last_leq_value(n_row_elements, row_indices, col);
            if (row_indices[current] == col)
            {
                element_count += 1;
            }
        }
    }
    return element_count;
}

jmtx_result jmtxds_matrix_crs_entries_in_col(const jmtxd_matrix_crs* mtx, uint32_t col, uint32_t* p_n)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTXD_TYPE_CRS)
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
    *p_n = jmtxd_matrix_crs_entries_in_col(mtx, col);
    return JMTX_RESULT_SUCCESS;
}

uint32_t
jmtxd_matrix_crs_get_col(const jmtxd_matrix_crs* mtx, uint32_t col, uint32_t n, double p_values[JMTX_ARRAY_ATTRIB(n)], uint32_t p_rows[JMTX_ARRAY_ATTRIB(n)])
{
    uint32_t k = 0;
    for (uint32_t row = 0; k < n && row < mtx->base.rows && (!row || (mtx->end_of_row_offsets[row - 1] != mtx->n_entries)); ++row)
    {
        uint32_t* row_indices;
        double* unused_row_values;
        const uint32_t n_row_elements = crs_get_row_entries(mtx, row, &row_indices, &unused_row_values);
        if (n_row_elements && row_indices[0] <= col && row_indices[n_row_elements - 1] >= col)
        {
            uint32_t current = jmtx_internal_find_last_leq_value(n_row_elements, row_indices, col);
            if (row_indices[current] == col)
            {
                p_values[k] = mtx->values[(row ? mtx->end_of_row_offsets[row - 1] : 0) + current];
                p_rows[k] = row;
                k += 1;
            }
        }
    }
    return k;
}

jmtx_result jmtxds_matrix_crs_get_col(
        const jmtxd_matrix_crs* mtx, uint32_t col, uint32_t n, uint32_t* p_count, double* p_values, uint32_t* p_rows)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTXD_TYPE_CRS)
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
    if (!p_rows)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    *p_count = jmtxd_matrix_crs_get_col(mtx, col, n, p_values, p_rows);
    return JMTX_RESULT_SUCCESS;
}

jmtx_result jmtxd_matrix_crs_transpose(
        const jmtxd_matrix_crs* mtx, jmtxd_matrix_crs** p_out, const jmtx_allocator_callbacks* allocator_callbacks)
{
    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }

    const uint32_t cols = mtx->base.cols;
    jmtxd_matrix_crs* out;
    jmtx_result res = jmtxd_matrix_crs_new(&out, mtx->base.rows, mtx->base.cols, mtx->n_entries, allocator_callbacks);
    if (res != JMTX_RESULT_SUCCESS)
    {
        return res;
    }
    if (!allocator_callbacks)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }

    uint32_t* col_counts = allocator_callbacks->alloc(allocator_callbacks->state, sizeof*col_counts * cols);
    if (col_counts == NULL)
    {
        jmtxd_matrix_crs_destroy(out);
        return JMTX_RESULT_BAD_ALLOC;
    }
    memset(col_counts, 0, sizeof*col_counts * cols);

    uint32_t* row_ends = out->end_of_row_offsets;
    for (uint32_t i = 0; i < mtx->n_entries; ++i)
    {
        col_counts[mtx->indices[i]] += 1;
    }
    row_ends[0] = col_counts[0];
    //  Compute cumsums for offsets
    for (uint32_t i = 1; i < cols; ++i)
    {
        row_ends[i] = col_counts[i] + row_ends[i-1];
        col_counts[i] = 0; //   Zero the row counts so that they can be reused later for counting bucket sizes
    }
    col_counts[0] = 0;
    col_counts[cols - 1] = 0;

    for (uint32_t row = 0; row < mtx->base.rows; ++row)
    {
        uint32_t* in_cols;
        double* in_vals;
        uint32_t n_col = crs_get_row_entries(mtx, row, &in_cols, &in_vals);

        for (uint32_t idx = 0; idx < n_col; ++idx)
        {
            const uint32_t col = in_cols[idx];
            const uint32_t ip = col > 0 ? row_ends[col-1] : 0;
            const uint32_t n_cols = col_counts[col];

            out->values[ip+n_cols] = in_vals[idx];
            out->indices[ip+n_cols] = row;
            col_counts[col] += 1;
        }
    }
    out->n_entries = mtx->n_entries;

    allocator_callbacks->free(allocator_callbacks->state, col_counts);
    *p_out = out;

    return JMTX_RESULT_SUCCESS;
}

jmtx_result jmtxds_matrix_crs_transpose(
        const jmtxd_matrix_crs* mtx, jmtxd_matrix_crs** p_out, const jmtx_allocator_callbacks* allocator_callbacks)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTXD_TYPE_CRS)
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
    return jmtxd_matrix_crs_transpose(mtx, p_out, allocator_callbacks);
}

jmtx_result jmtxd_matrix_crs_copy(const jmtxd_matrix_crs* mtx, jmtxd_matrix_crs** p_out, const jmtx_allocator_callbacks* allocator_callbacks)
{
    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }
    jmtxd_matrix_crs* const out = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*out));
    if (!out)
    {
        return JMTX_RESULT_BAD_ALLOC;
    }
    double* const elements = allocator_callbacks->alloc(allocator_callbacks->state, (mtx->n_entries) * sizeof (*elements));
    if (!elements)
    {
        allocator_callbacks->free(allocator_callbacks->state, out);
        return JMTX_RESULT_BAD_ALLOC;
    }
    uint32_t* const indices = allocator_callbacks->alloc(allocator_callbacks->state, (mtx->n_entries) * sizeof *indices);
    if (!indices)
    {
        allocator_callbacks->free(allocator_callbacks->state, out);
        allocator_callbacks->free(allocator_callbacks->state, elements);
        return JMTX_RESULT_BAD_ALLOC;
    }
    uint32_t* const cum_sum = allocator_callbacks->alloc(allocator_callbacks->state, (mtx->base.rows) * sizeof *cum_sum);
    if (!cum_sum)
    {
        allocator_callbacks->free(allocator_callbacks->state, out);
        allocator_callbacks->free(allocator_callbacks->state, indices);
        allocator_callbacks->free(allocator_callbacks->state, elements);
        return JMTX_RESULT_BAD_ALLOC;
    }

    memcpy(elements, mtx->values, sizeof* elements * mtx->n_entries);
    memcpy(indices, mtx->indices, sizeof* indices * mtx->n_entries);
    memcpy(cum_sum, mtx->end_of_row_offsets, sizeof* cum_sum * (mtx->base.rows));
    memcpy(out, mtx, sizeof *out);
    out->values = elements;
    out->indices = indices;
    out->end_of_row_offsets = cum_sum;
    out->base = mtx->base;
    out->capacity = out->n_entries;
    out->base.allocator_callbacks = *allocator_callbacks;
    *p_out = out;
    return JMTX_RESULT_SUCCESS;
}

jmtx_result jmtxds_matrix_crs_copy(const jmtxd_matrix_crs* mtx, jmtxd_matrix_crs** p_out, const jmtx_allocator_callbacks* allocator_callbacks)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTXD_TYPE_CRS)
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
    return jmtxd_matrix_crs_copy(mtx, p_out, allocator_callbacks);
}

jmtx_result jmtxd_matrix_crs_build_row(jmtxd_matrix_crs* mtx, uint32_t row, uint32_t n, const uint32_t indices[JMTX_ARRAY_ATTRIB(static n)], const double values[JMTX_ARRAY_ATTRIB(static n)])
{
    jmtx_result res = JMTX_RESULT_SUCCESS;
    const uint32_t required_capacity = (uint32_t)((int32_t)mtx->n_entries + (int32_t)n);
    if (mtx->capacity < required_capacity)
    {
        double* new_element_ptr = mtx->base.allocator_callbacks.realloc(mtx->base.allocator_callbacks.state, mtx->values, sizeof*mtx->values * (required_capacity + 1));
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

    const uint32_t offset = row ? mtx->end_of_row_offsets[row - 1] : 0;
    memcpy(mtx->values + offset, values, sizeof*values * n);
    memcpy(mtx->indices + offset, indices, sizeof*indices * n);

    mtx->end_of_row_offsets[row] = n + offset;
    mtx->n_entries += n;

    return res;
}

double jmtxd_matrix_crs_vector_multiply_row(const jmtxd_matrix_crs* mtx, const double* x, uint32_t i)
{
    uint32_t* indices;
    double* values;
    const uint32_t n_row = crs_get_row_entries(mtx, i, &indices, &values);
    double v = 0;
    for (uint32_t j = 0; j < n_row; ++j)
    {
        v += values[j] * x[indices[j]];
    }
    return v;
}

jmtx_result jmtxds_matrix_crs_vector_multiply_row(const jmtxd_matrix_crs* mtx, const double* restrict x, uint32_t i, double* restrict p_r)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTXD_TYPE_CRS)
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

    *p_r = jmtxd_matrix_crs_vector_multiply_row(mtx, x, i);
    return JMTX_RESULT_SUCCESS;
}

jmtx_result jmtxd_matrix_crs_add_to_entry(jmtxd_matrix_crs* mtx, uint32_t i, uint32_t j, double value)
{
    jmtx_result res;
    uint32_t* row_indices;
    double* row_values;
    const uint32_t n_row_elements = crs_get_row_entries(mtx, i, &row_indices, &row_values);
    //  Check if row has any values
    if (n_row_elements != 0)
    {
        //  Find first column entry less or equal to it
        const uint32_t possible = jmtx_internal_find_last_leq_value(n_row_elements, row_indices, j);
        if (row_indices[possible] == j)
        {
            row_values[possible] += value;
            return JMTX_RESULT_SUCCESS;
        }
        else
        {
            //  Figure out where to insert
            if (row_indices[possible] <= j)
            {
                //  Insert right after the *possible*
                res = crs_insert_entry_at(mtx, i, possible + 1, value, j);
            }
            else
            {
                //  Insert at the position of *possible*
                res = crs_insert_entry_at(mtx, i, possible, value, j);
            }
        }
    }
    else
    {
        //  Insert the new element in an empty row
        res = crs_insert_entry_at(mtx, i, 0, value, j);
    }

    return res;
}

jmtx_result jmtxds_matrix_crs_add_to_entry(jmtxd_matrix_crs* mtx, uint32_t i, uint32_t j, double value)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTXD_TYPE_CRS)
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
    return jmtxd_matrix_crs_add_to_entry(mtx, i, j, value);
}

void jmtxd_matrix_crs_zero_all_entries(const jmtxd_matrix_crs* mtx)
{
    memset(mtx->values, 0, sizeof(*mtx->values) * mtx->n_entries);
}

jmtx_result jmtxds_matrix_crs_zero_all_entries(const jmtxd_matrix_crs* mtx)
{

    if (!mtx)
    {
//        REPORT_ERROR_MESSAGE("Matrix pointer was null");
//        LEAVE_FUNCTION();
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTXD_TYPE_CRS)
    {
//        REPORT_ERROR_MESSAGE("Matrix was not compressed row sparse");
//        LEAVE_FUNCTION();
        return JMTX_RESULT_WRONG_TYPE;
    }

    jmtxd_matrix_crs_zero_all_entries(mtx);
    return JMTX_RESULT_SUCCESS;
}

void jmtxd_matrix_crs_set_all_entries(const jmtxd_matrix_crs* mtx, double x)
{
    for (double* ptr = mtx->values; ptr != mtx->values + mtx->n_entries; ++ptr)
    {
        *ptr = x;
    }
}

jmtx_result jmtxds_matrix_crs_set_all_entries(jmtxd_matrix_crs* mtx, double x)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTXD_TYPE_CRS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    jmtxd_matrix_crs_set_all_entries(mtx, x);

    return JMTX_RESULT_SUCCESS;
}

void jmtxd_matrix_crs_remove_row(jmtxd_matrix_crs* mtx, uint32_t row)
{
    uint32_t* row_indices;
    double* row_values;
    const uint32_t removed_entry_count = crs_get_row_entries(mtx, row, &row_indices, &row_values);
    const uint32_t end_offset = mtx->end_of_row_offsets[row];
    //  Check if row is not empty
    if (removed_entry_count != 0)
    {
        //  Check if there are any entries following the ones being removed
        const uint32_t following_entries = mtx->n_entries - end_offset;
        if (following_entries != 0)
        {
            //  Move other values to their new position
            memmove(mtx->values + end_offset - removed_entry_count, mtx->values + end_offset, sizeof(*mtx->values) * following_entries);
            memmove(mtx->indices + end_offset - removed_entry_count, mtx->indices + end_offset, sizeof(*mtx->indices) * following_entries);
        }
        for (uint32_t i = row; i < mtx->base.rows - 1; ++i)
        {
            mtx->end_of_row_offsets[i] = mtx->end_of_row_offsets[i + 1] - removed_entry_count;
        }
        mtx->n_entries -= removed_entry_count;
    }
    mtx->base.rows -= 1;
}

jmtx_result jmtxds_matrix_crs_remove_row(jmtxd_matrix_crs* mtx, uint32_t row)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTXD_TYPE_CRS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (mtx->base.rows <= row)
    {
        return JMTX_RESULT_INDEX_OUT_OF_BOUNDS;
    }
    jmtxd_matrix_crs_remove_row(mtx, row);
    return JMTX_RESULT_SUCCESS;
}

void jmtxd_matrix_crs_remove_column(jmtxd_matrix_crs* mtx, uint32_t col)
{
    //  i: tracks the current element to check
    //  j: how many values to be removed were found
    //  r: what row we are currently in
    uint32_t i, j, r = 0;
    while (mtx->end_of_row_offsets[r] == 0)
    {
        r += 1;
    }
    for (i = 0, j = 0; i < mtx->n_entries; ++i)
    {
        if (mtx->indices[i] == col)
        {
            j += 1;
        }
        else
        {
            if (mtx->indices[i] > col)
            {
                mtx->indices[i] -= 1;
            }
            if (j != 0)
            {
                mtx->indices[i - j] = mtx->indices[i];
                mtx->values[i - j] = mtx->values[i];
            }
        }
        while (r < mtx->base.cols && i + 1 == mtx->end_of_row_offsets[r])
        {
            assert(mtx->end_of_row_offsets[r] >= j);
            mtx->end_of_row_offsets[r] -= j;
            r += 1;
        }
    }
    assert(mtx->n_entries >= j);
    assert(r == mtx->base.rows);
    mtx->n_entries -= j;
    mtx->base.cols -= 1;
}

jmtx_result jmtxds_matrix_crs_remove_column(jmtxd_matrix_crs* mtx, uint32_t col)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTXD_TYPE_CRS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (mtx->base.cols <= col)
    {
        return JMTX_RESULT_INDEX_OUT_OF_BOUNDS;
    }
    jmtxd_matrix_crs_remove_column(mtx, col);
    return JMTX_RESULT_SUCCESS;
}

void jmtxd_matrix_crs_clear(jmtxd_matrix_crs* mtx)
{
    mtx->n_entries = 0;
    memset(mtx->end_of_row_offsets, 0, sizeof(*mtx->end_of_row_offsets) * mtx->base.rows);
}

jmtx_result jmtxds_matrix_crs_clear(jmtxd_matrix_crs* mtx)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTXD_TYPE_CRS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    jmtxd_matrix_crs_clear(mtx);
    return JMTX_RESULT_SUCCESS;
}

jmtx_result jmtxd_matrix_crs_join_vertically(jmtxd_matrix_crs** output, const jmtx_allocator_callbacks* allocator_callbacks,
                                            unsigned k, const jmtxd_matrix_crs* matrix_list[JMTX_ARRAY_ATTRIB(static k)])
{
    const uint32_t col_count = matrix_list[0]->base.cols;
    uint32_t n_rows = matrix_list[0]->base.rows;

    uint32_t element_count = matrix_list[0]->n_entries;
    //  First one is already accounted for, so no need to start at 0
    for (unsigned i = 1; i < k; ++i)
    {
        const jmtxd_matrix_crs* const e = matrix_list[i];
        if (e->base.cols != col_count)
        {
            return JMTX_RESULT_DIMS_MISMATCH;
        }
        n_rows += e->base.rows;
        element_count += e->n_entries;
    }

    jmtxd_matrix_crs* out;
    jmtx_result res = jmtxd_matrix_crs_new(&out, n_rows, col_count, element_count, allocator_callbacks);
    if (res != JMTX_RESULT_SUCCESS)
    {
        return res;
    }

    uint32_t pos = 0;
    for (unsigned i = 0; i < k; ++i)
    {
        const jmtxd_matrix_crs* const mtx = matrix_list[i];
        unsigned j;
        for (j = 0; j < mtx->base.rows; ++j)
        {
            uint32_t* indices;
            double* values;
            const uint32_t n = crs_get_row_entries(mtx, j, &indices, &values);
            //  All memory for this should be allocated in advance, so no need to check the return value
            const jmtx_result r1 = jmtxd_matrix_crs_build_row(out, pos + j, n, indices, values);
            assert(r1 == JMTX_RESULT_SUCCESS); (void)r1;
        }
        pos += j;
    }
    assert(pos == n_rows);

    *output = out;
    return JMTX_RESULT_SUCCESS;
}

jmtx_result jmtxd_matrix_crs_new_like(
        const jmtxd_matrix_crs* mtx, jmtxd_matrix_crs** p_out, const jmtx_allocator_callbacks* allocator_callbacks,
        const double* p_val)
{
    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }
    jmtxd_matrix_crs* const out = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*out));
    if (!out)
    {
        return JMTX_RESULT_BAD_ALLOC;
    }
    double* const elements = allocator_callbacks->alloc(allocator_callbacks->state, (mtx->n_entries) * sizeof (*elements));
    if (!elements)
    {
        allocator_callbacks->free(allocator_callbacks->state, out);
        return JMTX_RESULT_BAD_ALLOC;
    }
    uint32_t* const indices = allocator_callbacks->alloc(allocator_callbacks->state, (mtx->n_entries) * sizeof *indices);
    if (!indices)
    {
        allocator_callbacks->free(allocator_callbacks->state, out);
        allocator_callbacks->free(allocator_callbacks->state, elements);
        return JMTX_RESULT_BAD_ALLOC;
    }
    uint32_t* const cum_sum = allocator_callbacks->alloc(allocator_callbacks->state, (mtx->base.rows) * sizeof *cum_sum);
    if (!cum_sum)
    {
        allocator_callbacks->free(allocator_callbacks->state, out);
        allocator_callbacks->free(allocator_callbacks->state, indices);
        allocator_callbacks->free(allocator_callbacks->state, elements);
        return JMTX_RESULT_BAD_ALLOC;
    }

    if (p_val)
    {
        const double v = *p_val;
        if (v == 0)
        {
            memset(elements, 0, sizeof* elements * mtx->n_entries);
        }
        else
        {
            for (uint32_t i = 0; i < mtx->n_entries; ++i)
            {
                elements[i] = v;
            }
        }
    }

    memcpy(indices, mtx->indices, sizeof* indices * mtx->n_entries);
    memcpy(cum_sum, mtx->end_of_row_offsets, sizeof* cum_sum * (mtx->base.rows));
    memcpy(out, mtx, sizeof *out);
    out->values = elements;
    out->indices = indices;
    out->end_of_row_offsets = cum_sum;
    out->base = mtx->base;
    out->base.allocator_callbacks = *allocator_callbacks;
    *p_out = out;
    return JMTX_RESULT_SUCCESS;
}

jmtx_result jmtxds_matrix_crs_new_like(
        const jmtxd_matrix_crs* mtx, jmtxd_matrix_crs** p_out, const jmtx_allocator_callbacks* allocator_callbacks,
        const double* p_val)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTXD_TYPE_CRS)
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
    if (p_val && !isfinite(*p_val))
    {
        return JMTX_RESULT_BAD_PARAM;
    }
    return jmtxd_matrix_crs_new_like(mtx, p_out, allocator_callbacks, p_val);
}
uint32_t jmtxd_matrix_crs_find_upper_bandwidth(const jmtxd_matrix_crs* mtx)
{
    //  Find the greatest distance above the main diagonal
    uint_fast32_t v_max = 0;
    for (uint_fast32_t i = 0, p = 0; i < mtx->base.rows; ++i)
    {
        for (p = 0; p < mtx->end_of_row_offsets[i]; ++p)
        {
            if (mtx->indices[p] > i)
            {
                const uint_fast32_t dif = mtx->indices[p] - i;
                if (dif > v_max)
                {
                    v_max = dif;
                }
            }
        }
    }
    return v_max;
}

/**
 * Finds the lower bandwidth of the matrix; what is the furthest distance of and entry bellow the main diagonal
 * @param mtx matrx to find the lower bandwidth of
 * @return lower bandwidth of the matrix
 */
uint32_t jmtxd_matrix_crs_find_lower_bandwidth(const jmtxd_matrix_crs* mtx)
{
    //  Find the greatest distance above the main diagonal
    uint_fast32_t v_max = 0;
    for (uint_fast32_t i = 0, p = 0; i < mtx->base.rows; ++i)
    {
        for (; p < mtx->end_of_row_offsets[i]; ++p)
        {
            if (mtx->indices[p] < i)
            {
                const uint_fast32_t dif = i - mtx->indices[p];
                if (dif > v_max)
                {
                    v_max = dif;
                }
            }
        }
    }
    return v_max;
}

jmtx_result jmtxds_matrix_crs_join_vertically(
        jmtxd_matrix_crs** output, const jmtx_allocator_callbacks* allocator_callbacks, unsigned int k,
        const jmtxd_matrix_crs** matrix_list)
{
    if (!output)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!matrix_list)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (allocator_callbacks && (!allocator_callbacks->alloc || !allocator_callbacks->free || !allocator_callbacks->realloc))
    {
        return JMTX_RESULT_BAD_PARAM;
    }

    const uint32_t col_count = matrix_list[0]->base.cols;
    for (const jmtxd_matrix_crs** pmtx = matrix_list; pmtx != matrix_list + k; ++pmtx)
    {
        const jmtxd_matrix_crs* mtx = *pmtx;
        if (!mtx)
        {
            return JMTX_RESULT_NULL_PARAM;
        }
        if (mtx->base.type != JMTXD_TYPE_CRS)
        {
            return JMTX_RESULT_WRONG_TYPE;
        }
        if (mtx->base.cols != col_count)
        {
            return JMTX_RESULT_DIMS_MISMATCH;
        }
    }
    
    return jmtxd_matrix_crs_join_vertically(output, allocator_callbacks, k, matrix_list);
}
