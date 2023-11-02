//
// Created by jan on 15.6.2022.
//

#include <assert.h>
#include "sparse_column_compressed.h"
#include "sparse_column_compressed_internal.h"
//
// Created by jan on 13.6.2022.
//



#define DEFAULT_RESERVED_ELEMENTS 64

#define FILL_VALUE 0xDEADBEEF

static void beef_it_up(float* ptr, size_t elements)
{
    if (elements < 2)
    {
        const uint32_t BeefBuffer2[2] = {0xDEADBEEF, 0xDEADBEEF};
        memcpy(ptr, BeefBuffer2, elements * sizeof(*ptr));
    }

    while (elements & 0x3)
    {
        *(uint32_t*)ptr = 0xDEADBEEF;
        --elements;
        ++ptr;
    }

    for (uint32_t i = 0; i < elements / 4; ++i)
    {
        *((uint32_t*)(ptr + i)) = 0xDEADBEEF;
    }
    memcpy(ptr + elements / 4, ptr, elements / 4 * sizeof(*ptr));
    memcpy(ptr + elements / 2, ptr, elements / 2 * sizeof(*ptr));
    //  Beefed
}

static int beef_check(const float* ptr, size_t elements)
{
    const uint32_t* const buffer = (const uint32_t*)ptr;
    int beef_count = 0;
    for (uint32_t i = 0; i < elements; ++i)
    {
        beef_count += (buffer[i] == 0xDEADBEEF);
    }
    return beef_count;
}

static uint32_t ccs_get_column_entries(const jmtx_matrix_ccs* mtx, uint32_t col, uint32_t** pp_indices, float** pp_values)
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
 * @return
 */
static jmtx_result ccs_insert_entry_at(jmtx_matrix_ccs* mtx, uint32_t col, uint32_t position, float value, uint32_t index)
{
    const uint32_t global_position = position + (col ? mtx->end_of_column_offsets[col - 1] : 0);
    assert(position <= mtx->n_entries);
    assert(!col || (mtx->end_of_column_offsets[col - 1] == global_position || mtx->indices[global_position - 1] < index));

    if (mtx->capacity == mtx->n_entries)
    {
        //  Reallocate arrays
        const size_t new_capacity = mtx->capacity + DEFAULT_RESERVED_ELEMENTS;
        float* const new_values = mtx->base.allocator_callbacks.realloc(mtx->base.allocator_callbacks.state, mtx->values, sizeof(*mtx->values) * new_capacity);
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

#if !(defined(JMTX_NO_VERIFY_PARAMS) || defined(NDEBUG))
        //  Beef it up!
        beef_it_up(new_values + mtx->capacity, new_capacity - mtx->capacity);
        beef_it_up((float*)(new_indices + mtx->capacity), new_capacity - mtx->capacity);
#endif
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



jmtx_result jmtx_matrix_ccs_new(
        jmtx_matrix_ccs** mtx, uint32_t columns, uint32_t rows, uint32_t reserved_entries,
        const jmtx_allocator_callbacks* allocator_callbacks)
{
    if (!mtx)
    {
//        //REPORT_ERROR_MESSAGE("Matrix pointer was null");
//        //LEAVE_FUNCTION();
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!rows)
    {
//        //REPORT_ERROR_MESSAGE("Number of rows was 0");
//        //LEAVE_FUNCTION();
        return JMTX_RESULT_BAD_PARAM;
    }
    if (!columns)
    {
//        //REPORT_ERROR_MESSAGE("Number of columns was 0");
//        //LEAVE_FUNCTION();
        return JMTX_RESULT_BAD_PARAM;
    }
    if (reserved_entries > columns * rows)
    {
//        //REPORT_ERROR_MESSAGE("Number of reserved values (%u) exceeds product of columns (%u) by rows (%u)", reserved_entries, rows, columns);
//        //LEAVE_FUNCTION();
        return JMTX_RESULT_BAD_PARAM;
    }
    if (!allocator_callbacks)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }
    else if (!allocator_callbacks->alloc || !allocator_callbacks->realloc || !allocator_callbacks->free)
    {
        return JMTX_RESULT_BAD_PARAM;
    }

    if (reserved_entries == 0)
    {
        reserved_entries = DEFAULT_RESERVED_ELEMENTS;
        reserved_entries = reserved_entries < columns * rows ? reserved_entries : columns * rows;
    }

    uint32_t* offsets = NULL;
    uint32_t* indices = NULL;
    float* values = allocator_callbacks->alloc(allocator_callbacks->state, (reserved_entries) * sizeof(*values));
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

    offsets = allocator_callbacks->alloc(allocator_callbacks->state, (columns) * sizeof*offsets);
    if (!offsets)
    {
        allocator_callbacks->free(allocator_callbacks->state, indices);
        allocator_callbacks->free(allocator_callbacks->state, values);
        return JMTX_RESULT_BAD_ALLOC;
    }
    memset(offsets, 0, (columns) * sizeof(*offsets));

#if !(defined(JMTX_NO_VERIFY_PARAMS) || defined(NDEBUG))
    beef_it_up(values, reserved_entries);
    static_assert(sizeof(float) == sizeof(uint32_t), "element and index sizes must be the same");
    beef_it_up((float*)indices, reserved_entries);
#endif

    jmtx_matrix_ccs* const this = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*this));
    if (!this)
    {
        allocator_callbacks->free(allocator_callbacks->state, offsets);
        allocator_callbacks->free(allocator_callbacks->state, indices);
        allocator_callbacks->free(allocator_callbacks->state, values);
        return JMTX_RESULT_BAD_ALLOC;
    }

    this->base.rows = rows;
    this->base.type = JMTX_TYPE_CCS;
    this->base.cols = columns;
    this->base.allocator_callbacks = *allocator_callbacks;
    this->indices = indices;
    this->values = values;
    this->capacity = reserved_entries;
    this->n_entries = 0;
    this->end_of_column_offsets = offsets;

    *mtx = this;
    return JMTX_RESULT_SUCCESS;
}

jmtx_result jmtx_matrix_ccs_destroy(jmtx_matrix_ccs* mtx)
{
#ifndef JMTX_NO_VERIFY_PARAMS
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTX_TYPE_CCS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
#endif

    mtx->base.allocator_callbacks.free(mtx->base.allocator_callbacks.state, mtx->indices);
    mtx->base.allocator_callbacks.free(mtx->base.allocator_callbacks.state, mtx->end_of_column_offsets);
    mtx->base.allocator_callbacks.free(mtx->base.allocator_callbacks.state, mtx->values);
    mtx->base.allocator_callbacks.free(mtx->base.allocator_callbacks.state, mtx);

    return JMTX_RESULT_SUCCESS;
}

jmtx_result jmtx_matrix_ccs_shrink(jmtx_matrix_ccs* mtx)
{
#ifndef JMTX_NO_VERIFY_PARAMS
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTX_TYPE_CCS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
#endif
    if (mtx->n_entries == mtx->capacity)
    {
        return JMTX_RESULT_SUCCESS;
    }
    mtx->capacity = mtx->n_entries;
    float* element_new_ptr = mtx->base.allocator_callbacks.realloc(mtx->base.allocator_callbacks.state, mtx->values, sizeof*mtx->values * (mtx->n_entries));
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

    return JMTX_RESULT_SUCCESS;
}

jmtx_result jmtx_matrix_ccs_set_col(jmtx_matrix_ccs* mtx, uint32_t col, uint32_t n, const uint32_t* indices, const float* values)
{
#ifndef JMTX_NO_VERIFY_PARAMS
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTX_TYPE_CCS)
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
#endif

    const uint32_t beginning_offset = col ? mtx->end_of_column_offsets[col - 1] : 0;
    const int32_t new_elements = (int32_t)n - (int32_t)(mtx->end_of_column_offsets[col] - beginning_offset);
    const uint32_t required_capacity = (uint32_t)((int32_t)mtx->n_entries + new_elements);
    if (mtx->capacity < required_capacity)
    {
        float* new_element_ptr = mtx->base.allocator_callbacks.realloc(mtx->base.allocator_callbacks.state, mtx->values, sizeof*(mtx->values) * (required_capacity + 1));
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

jmtx_result jmtx_matrix_ccs_vector_multiply(const jmtx_matrix_ccs* mtx, const float* restrict x, float* restrict y)
{
#ifndef JMTX_NO_VERIFY_PARAMS
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTX_TYPE_CCS)
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
#endif

    for (uint32_t i = 0; i < mtx->base.cols; ++i)
    {
        uint32_t* indices;
        float* values;
        const uint32_t n_elements = ccs_get_column_entries(mtx, i, &indices, &values);
        float v = 0;
        for (uint32_t j = 0; j < n_elements; ++j)
        {
            const uint32_t k = indices[j];
            v += values[j] * x[k];
        }
        y[i] = v;
    }
    return JMTX_RESULT_SUCCESS;
}

jmtx_result jmtx_matrix_ccs_set_element(jmtx_matrix_ccs* mtx, uint32_t i, uint32_t j, float value)
{
#ifndef JMTX_NO_VERIFY_PARAMS
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTX_TYPE_CCS)
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
#endif

    jmtx_result res;
    uint32_t* col_indices;
    float* col_values;
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

jmtx_result jmtx_matrix_ccs_get_element(const jmtx_matrix_ccs* mtx, uint32_t i, uint32_t j, float* p_value)
{
#ifndef JMTX_NO_VERIFY_PARAMS
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTX_TYPE_CCS)
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
#endif

    uint32_t* col_indices;
    float* col_values;
    const uint32_t n_col_elements = ccs_get_column_entries(mtx, j, &col_indices, &col_values);
    //  Check if row has any values
    if (n_col_elements != 0)
    {
        //  Find first column entry less or equal to it
        const uint32_t possible = jmtx_internal_find_last_leq_value(n_col_elements, col_indices, i);
        if (col_indices[possible] == i)
        {
            *p_value = col_values[possible];
            return JMTX_RESULT_SUCCESS;
        }
    }
    *p_value = 0;
//    LEAVE_FUNCTION();
    return JMTX_RESULT_SUCCESS;
}

jmtx_result jmtx_matrix_ccs_get_col(const jmtx_matrix_ccs* mtx, uint32_t col, uint32_t* n, uint32_t** p_indices, float** p_elements)
{
#ifndef JMTX_NO_VERIFY_PARAMS
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTX_TYPE_CCS)
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
#endif

    uint32_t* col_indices;
    float* col_values;
    const uint32_t n_row_elements = ccs_get_column_entries(mtx, col, &col_indices, &col_values);
    *p_indices = col_indices;
    *p_elements = col_values;
    *n = n_row_elements;

    return JMTX_RESULT_SUCCESS;
}

jmtx_result jmtx_matrix_ccs_beef_check(const jmtx_matrix_ccs* mtx, int* p_beef_status)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTX_TYPE_CCS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (!p_beef_status)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }

    const int beef_status = beef_check(mtx->values, mtx->n_entries) + beef_check((void*)mtx->indices, mtx->n_entries);
    *p_beef_status = (beef_status << 16) | 0x0000Beef;

    return JMTX_RESULT_SUCCESS;
}

jmtx_result jmtx_matrix_ccs_apply_unary_fn(const jmtx_matrix_ccs* mtx, int (*unary_fn)(uint32_t i, uint32_t j, float* p_element, void* param), void* param)
{
#ifndef JMTX_NO_VERIFY_PARAMS

    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTX_TYPE_CCS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (!unary_fn)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
#endif


    for (uint32_t j = 0; j < mtx->base.cols; ++j)
    {
        float* p_elements;
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

jmtx_result jmtx_matrix_ccs_remove_zeros(jmtx_matrix_ccs* mtx)
{
#ifndef JMTX_NO_VERIFY_PARAMS
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTX_TYPE_CCS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
#endif
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


    //  Beef
    mtx->n_entries -= c;

#if !(defined(JMTX_NO_VERIFY_PARAMS) || defined(NDEBUG))
    beef_it_up(mtx->values + mtx->n_entries, c);
    static_assert(sizeof(float) == sizeof(uint32_t), "Size of index and scalar must be the same for beef");
    beef_it_up((float*)(mtx->indices + mtx->n_entries), c);
#endif
//    LEAVE_FUNCTION();
    return 0;
}

static inline bool is_bellow_magnitude(float x, float mag)
{
    //  Reinterpret both as uints
    uint32_t u_x = *(uint32_t*)&x;
    const uint32_t u_mag = *(uint32_t*)&mag;
    //  Mask out the zero bit of x
    u_x &= ~(1u << 31u);
    //  Compare as if they were ints (since float comparison works the same)
    return u_x < u_mag;
}

jmtx_result jmtx_matrix_ccs_remove_bellow(jmtx_matrix_ccs* mtx, float v)
{
#ifndef JMTX_NO_VERIFY_PARAMS
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTX_TYPE_CCS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
#endif
//  Update offsets
    uint32_t p, c = 0, r = 0;
    while (mtx->end_of_column_offsets[r] == 0)
    {
        r += 1;
    }
    for (p = 0; p < mtx->n_entries; ++p)
    {
        assert(r < mtx->base.cols);
        if (is_bellow_magnitude(mtx->values[p], v))
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
        if (is_bellow_magnitude(mtx->values[p0 - 1], v))
        {
            uint32_t p1 = p0;
            while (p0 != 0 &&
                   is_bellow_magnitude(mtx->values[p0 - 1], v))
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
#ifndef JMTX_NO_VERIFY_PARAMS
    beef_it_up(mtx->values + mtx->n_entries, c);
    static_assert(sizeof(float) == sizeof(uint32_t), "Size of index and scalar must be the same for beef");
    beef_it_up((float*)(mtx->indices + mtx->n_entries), c);
#endif
//    LEAVE_FUNCTION();
    return JMTX_RESULT_SUCCESS;
}

jmtx_result jmtx_matrix_ccs_elements_in_row(const jmtx_matrix_ccs* mtx, uint32_t row, uint32_t* p_n)
{
#ifndef JMTX_NO_VERIFY_PARAMS
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTX_TYPE_CCS)
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
#endif

    uint32_t element_count = 0;
    for (uint32_t col = 0; col < mtx->base.cols && (!col || (mtx->end_of_column_offsets[col - 1] != mtx->n_entries)); ++col)
    {
        uint32_t* col_indices;
        float* unused_col_values;
        const uint32_t n_col_elements = ccs_get_column_entries(mtx, col, &col_indices, &unused_col_values);
        if (n_col_elements)
        {
            uint32_t current = jmtx_internal_find_last_leq_value(n_col_elements, col_indices, row);
            if (col_indices[current] == row)
            {
                element_count += 1;
            }
        }
    }
    *p_n = element_count;
    return JMTX_RESULT_SUCCESS;
}

jmtx_result jmtx_matrix_ccs_get_row(
        const jmtx_matrix_ccs* mtx, uint32_t row, uint32_t n, float* p_values, uint32_t* p_count, uint32_t* p_columns)
{
#ifndef JMTX_NO_VERIFY_PARAMS
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTX_TYPE_CCS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (row >= mtx->base.rows)
    {
        return JMTX_RESULT_INDEX_OUT_OF_BOUNDS;
    }
    if (!p_values)
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
#endif
    uint32_t k = 0;
    for (uint32_t col = 0; k < n && col < mtx->base.cols && (!col || (mtx->end_of_column_offsets[col - 1] != mtx->n_entries)); ++col)
    {
        uint32_t* col_indices;
        float* unused_col_values;
        const uint32_t n_col_elements = ccs_get_column_entries(mtx, col, &col_indices, &unused_col_values);
        if (n_col_elements)
        {
            uint32_t current = jmtx_internal_find_last_leq_value(n_col_elements, col_indices, row);
            if (col_indices[current] == row)
            {
                p_values[k] = mtx->values[(col ? mtx->end_of_column_offsets[col - 1] : 0) + current];
                p_columns[k] = col;
                k += 1;
            }
        }
    }
    if (p_count)
    {
        *p_count = k;
    }
    return JMTX_RESULT_SUCCESS;
}

jmtx_result jmtx_matrix_ccs_transpose(const jmtx_matrix_ccs* mtx, jmtx_matrix_ccs** p_out,
                                      const jmtx_allocator_callbacks* allocator_callbacks)
{
//    CALL_FUNCTION(matrix_crs_transpose);
#ifndef JMTX_NO_VERIFY_PARAMS
    if (!mtx)
    {
//        REPORT_ERROR_MESSAGE("Input matrix pointer was null");
//        LEAVE_FUNCTION();
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTX_TYPE_CCS)
    {
//        REPORT_ERROR_MESSAGE("Input matrix was not compressed row sparse");
//        LEAVE_FUNCTION();
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (!p_out)
    {
//        REPORT_ERROR_MESSAGE("Output matrix pointer was null");
//        LEAVE_FUNCTION();
        return JMTX_RESULT_NULL_PARAM;
    }
#endif
    if (!allocator_callbacks)
    {
        allocator_callbacks = &mtx->base.allocator_callbacks;
    }
    else if (!allocator_callbacks->alloc || !allocator_callbacks->realloc || !allocator_callbacks->free)
    {
        return JMTX_RESULT_BAD_PARAM;
    }
    const uint32_t n_elements = mtx->n_entries;
    const uint32_t new_cols = mtx->base.rows;
    const uint32_t new_rows = mtx->base.cols;

    jmtx_matrix_ccs* const out = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*out));
    if (!out)
    {
        return JMTX_RESULT_BAD_ALLOC;
    }

    uint32_t* const row_cum_counts = mtx->base.allocator_callbacks.alloc(mtx->base.allocator_callbacks.state, (new_cols + 1) * sizeof(*row_cum_counts));
    if (!row_cum_counts)
    {
        allocator_callbacks->free(allocator_callbacks->state, out);
        return JMTX_RESULT_BAD_ALLOC;
    }
    uint32_t* const new_indices = mtx->base.allocator_callbacks.alloc(mtx->base.allocator_callbacks.state, (n_elements + 1) * sizeof(*new_indices));
    if (!new_indices)
    {
        mtx->base.allocator_callbacks.free(mtx->base.allocator_callbacks.state, row_cum_counts);
        allocator_callbacks->free(allocator_callbacks->state, out);
        return JMTX_RESULT_BAD_ALLOC;
    }
    memset(new_indices, 0, (n_elements + 1) * sizeof*new_indices);
    float* const new_elements = mtx->base.allocator_callbacks.alloc(mtx->base.allocator_callbacks.state, (n_elements + 1) * sizeof*new_elements);
    if (!new_elements)
    {
        mtx->base.allocator_callbacks.free(mtx->base.allocator_callbacks.state, row_cum_counts);
        mtx->base.allocator_callbacks.free(mtx->base.allocator_callbacks.state, new_indices);
        allocator_callbacks->free(allocator_callbacks->state, out);
        return JMTX_RESULT_BAD_ALLOC;
    }

    *row_cum_counts = 1;

    for (uint32_t j = 0, n = 0, p = 1; j < mtx->base.rows; ++j)
    {
        n = 0;
        //  This MUST NOT fail, since the parameters are all within the bounds
        jmtx_result res = jmtx_matrix_ccs_elements_in_row(mtx, j, &n);
        assert(res == JMTX_RESULT_SUCCESS);
        (void) res;
        uint32_t c;
        res = jmtx_matrix_ccs_get_row(mtx, j, n, new_elements + p, &c, new_indices + p);
        assert(res == JMTX_RESULT_SUCCESS);
        assert(c == n);
        (void) res;
        p += n;
        row_cum_counts[j + 1] = row_cum_counts[j] + n;
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

jmtx_result jmtx_matrix_ccs_copy(const jmtx_matrix_ccs* mtx, jmtx_matrix_ccs** p_out, const jmtx_allocator_callbacks* allocator_callbacks)
{
#ifndef JMTX_NO_VERIFY_PARAMS

    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTX_TYPE_CCS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (!p_out)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
#endif
    if (!allocator_callbacks)
    {
        allocator_callbacks = &mtx->base.allocator_callbacks;
    }
    else if (!allocator_callbacks->alloc || !allocator_callbacks->realloc || !allocator_callbacks->free)
    {
        return JMTX_RESULT_BAD_PARAM;
    }
    jmtx_matrix_ccs* const this = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*this));
    if (!this)
    {
        return JMTX_RESULT_BAD_ALLOC;
    }

    float* const elements = mtx->base.allocator_callbacks.alloc(mtx->base.allocator_callbacks.state, (mtx->n_entries) * sizeof *elements);
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

    memcpy(elements + 1, mtx->values + 1, sizeof* elements * mtx->n_entries);
    memcpy(indices + 1, mtx->indices + 1, sizeof* indices * mtx->n_entries);
    memcpy(cum_sum, mtx->end_of_column_offsets, sizeof* cum_sum * (mtx->base.cols));
    this->base = mtx->base;
    this->values = elements;
    this->indices = indices;
    this->end_of_column_offsets = cum_sum;
    *p_out = this;
    return JMTX_RESULT_SUCCESS;
}

jmtx_result jmtx_matrix_ccs_zero_all_elements(jmtx_matrix_ccs* mtx)
{
#ifndef JMTX_NO_VERIFY_PARAMS
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTX_TYPE_CCS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
#endif
    memset(mtx->values, 0, sizeof(*mtx->values) * mtx->n_entries);

    return JMTX_RESULT_SUCCESS;
}

jmtx_result jmtx_matrix_ccs_set_all_elements(jmtx_matrix_ccs* mtx, float x)
{
#ifndef JMTX_NO_VERIFY_PARAMS
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTX_TYPE_CCS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
#endif
    for (float* ptr = mtx->values; ptr != mtx->values + mtx->n_entries; ++ptr)
    {
        *ptr = x;
    }

    return JMTX_RESULT_SUCCESS;
}

jmtx_result
jmtx_matrix_ccs_build_col(jmtx_matrix_ccs* mtx, uint32_t col, uint32_t n, const uint32_t* indices, const float* values)
{
    const uint32_t required_capacity = (uint32_t)((int32_t)mtx->n_entries + (int32_t)n);
    if (mtx->capacity < required_capacity)
    {
        float* new_element_ptr = mtx->base.allocator_callbacks.realloc(mtx->base.allocator_callbacks.state, mtx->values, sizeof*mtx->values * (required_capacity + 1));
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
    }

    const uint32_t offset = col ? mtx->end_of_column_offsets[col - 1] : 0;
    memcpy(mtx->values + offset, values, sizeof*values * n);
    memcpy(mtx->indices + offset, indices, sizeof*indices * n);

    mtx->end_of_column_offsets[col] = n + offset;
    mtx->n_entries += n;

    return JMTX_RESULT_SUCCESS;
}
