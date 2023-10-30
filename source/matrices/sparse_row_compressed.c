//
// Created by jan on 13.6.2022.
//

#include <assert.h>
#include "sparse_row_compressed.h"
#include "sparse_row_compressed_internal.h"

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

static uint32_t crs_get_row_entries(const jmtx_matrix_crs* mtx, uint32_t row, uint32_t** pp_indices, float** pp_values)
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
 * @return
 */
static jmtx_result crs_insert_entry_at(jmtx_matrix_crs* mtx, uint32_t row, uint32_t position, float value, uint32_t index)
{
    const uint32_t global_position = position + (row ? mtx->end_of_row_offsets[row - 1] : 0);
    assert(position <= mtx->n_entries);
    assert(!row || (mtx->end_of_row_offsets[row - 1] == global_position || mtx->indices[global_position - 1] < index));

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

#ifndef JMTX_NO_VERIFY_PARAMS
        //  Beef it up!
        beef_it_up(new_values + mtx->capacity, new_capacity - mtx->capacity);
        beef_it_up((float*)(new_indices + mtx->capacity), new_capacity - mtx->capacity);
#endif
        mtx->capacity = new_capacity;
    }

    //  Check for numer of elements after the position
    const uint32_t elements_after = mtx->n_entries - global_position;
    if (elements_after)
    {
        //  Move other elements out of the way
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

jmtx_result jmtx_matrix_crs_new(
        jmtx_matrix_crs** p_mtx, uint32_t cols, uint32_t rows, uint32_t reserved_entries,
        const jmtx_allocator_callbacks* allocator_callbacks)
{
#ifndef JMTX_NO_VERIFY_PARAMS
    if (!p_mtx)
    {
//        REPORT_ERROR_MESSAGE("Matrix pointer was null");
//        LEAVE_FUNCTION();
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!rows)
    {
//        REPORT_ERROR_MESSAGE("Number of rows was 0");
//        LEAVE_FUNCTION();
        return JMTX_RESULT_BAD_PARAM;
    }
    if (!cols)
    {
//        REPORT_ERROR_MESSAGE("Number of columns was 0");
//        LEAVE_FUNCTION();
        return JMTX_RESULT_BAD_PARAM;
    }
    if (reserved_entries > cols * rows)
    {
//        REPORT_ERROR_MESSAGE("Number of reserved elements (%u) exceeds product of columns (%u) by rows (%u)", reserved_entries, rows, columns);
//        LEAVE_FUNCTION();
        return JMTX_RESULT_BAD_PARAM;
    }
#endif
    if (reserved_entries == 0)
    {
        reserved_entries = DEFAULT_RESERVED_ELEMENTS;
        reserved_entries = reserved_entries < cols * rows ? reserved_entries : cols * rows;
    }
    if (!allocator_callbacks)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }
    else if (!allocator_callbacks->alloc || !allocator_callbacks->realloc || !allocator_callbacks->free)
    {
        return JMTX_RESULT_BAD_PARAM;
    }
    jmtx_result mtx_res = 0;
    uint32_t* offsets = NULL;
    uint32_t* indices = NULL;

    jmtx_matrix_crs* mtx = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*mtx));
    if (!mtx)
    {
        return JMTX_RESULT_BAD_ALLOC;
    }

    float* values = allocator_callbacks->alloc(allocator_callbacks->state, (reserved_entries) * sizeof(*values));
    if (!values)
    {
        mtx_res = JMTX_RESULT_BAD_ALLOC;
//        CALLOC_FAILED((1 + reserved_entries) * sizeof*p_elements);
        goto fail1;
    }
    memset(values, 0, (reserved_entries) * sizeof(*values));

    indices = allocator_callbacks->alloc(allocator_callbacks->state, (reserved_entries) * sizeof(*indices));
    if (!indices)
    {
        mtx_res = JMTX_RESULT_BAD_ALLOC;
//        CALLOC_FAILED((1 + reserved_entries) * sizeof*indices);
        goto fail2;
    }
    memset(indices, 0, (reserved_entries) * sizeof(*indices));

    offsets = allocator_callbacks->alloc(allocator_callbacks->state, (rows) * sizeof(*offsets));
    if (!offsets)
    {
        mtx_res = JMTX_RESULT_BAD_ALLOC;
//        CALLOC_FAILED((columns + 1) * sizeof*offsets);
        goto fail3;
    }
    memset(offsets, 0, (rows) * sizeof(*offsets));

#ifndef JMTX_NO_VERIFY_PARAMS
    beef_it_up(values, reserved_entries);
    static_assert(sizeof(float) == sizeof(uint32_t), "element and index sizes must be the same");
    beef_it_up((float*)indices, reserved_entries);
#endif
    memset(mtx, 0, sizeof*mtx);
    mtx->base.cols = cols;
    mtx->base.type = JMTX_TYPE_CRS;
    mtx->base.rows = rows;
    mtx->base.allocator_callbacks = *allocator_callbacks;
    mtx->indices = indices;
    mtx->values = values;
    mtx->capacity = reserved_entries;
    mtx->n_entries = 0;
    mtx->end_of_row_offsets = offsets;
    *p_mtx = mtx;
//    LEAVE_FUNCTION();
    return mtx_res;
    fail3: allocator_callbacks->free(allocator_callbacks->state, indices);
    fail2: allocator_callbacks->free(allocator_callbacks->state, offsets);
    fail1: allocator_callbacks->free(allocator_callbacks->state, values);
    allocator_callbacks->free(allocator_callbacks->state, mtx);
//    LEAVE_FUNCTION();
    return mtx_res;
}

jmtx_result jmtx_matrix_crs_destroy(jmtx_matrix_crs* mtx)
{
#ifndef JMTX_NO_VERIFY_PARAMS
    if (!mtx)
    {
//        REPORT_ERROR_MESSAGE("Matrix pointer was null");
//        LEAVE_FUNCTION();
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTX_TYPE_CRS)
    {
//        REPORT_ERROR_MESSAGE("Matrix was not compressed row sparse");
//        LEAVE_FUNCTION();
        return JMTX_RESULT_WRONG_TYPE;
    }
#endif
    mtx->base.allocator_callbacks.free(mtx->base.allocator_callbacks.state, mtx->indices);
    mtx->base.allocator_callbacks.free(mtx->base.allocator_callbacks.state, mtx->end_of_row_offsets);
    mtx->base.allocator_callbacks.free(mtx->base.allocator_callbacks.state, mtx->values);
    mtx->base.allocator_callbacks.free(mtx->base.allocator_callbacks.state, mtx);
//    LEAVE_FUNCTION();
    return JMTX_RESULT_SUCCESS;
}

jmtx_result jmtx_matrix_crs_shrink(jmtx_matrix_crs* mtx)
{
#ifndef JMTX_NO_VERIFY_PARAMS
//    CALL_FUNCTION(matrix_crs_shrink);
    if (!mtx)
    {
//        REPORT_ERROR_MESSAGE("Matrix pointer was null");
//        LEAVE_FUNCTION();
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTX_TYPE_CRS)
    {
//        REPORT_ERROR_MESSAGE("Matrix was not compressed row sparse");
//        LEAVE_FUNCTION();
        return JMTX_RESULT_WRONG_TYPE;
    }
#endif
    jmtx_result res = JMTX_RESULT_SUCCESS;

    if (mtx->n_entries == mtx->capacity)
    {
//        LEAVE_FUNCTION();
        return JMTX_RESULT_SUCCESS;
    }

    float* element_new_ptr = mtx->base.allocator_callbacks.realloc(mtx->base.allocator_callbacks.state, mtx->values, sizeof*mtx->values * (mtx->n_entries));
    if (!element_new_ptr)
    {
        res = JMTX_RESULT_BAD_ALLOC;
//        REALLOC_FAILED(sizeof*mtx->elements * (mtx->n_elements));
        goto end;
    }
    mtx->values = element_new_ptr;
    uint32_t* new_indices_ptr = mtx->base.allocator_callbacks.realloc(mtx->base.allocator_callbacks.state, mtx->indices, sizeof*mtx->indices * (mtx->n_entries));
    if (!new_indices_ptr)
    {
        res = JMTX_RESULT_BAD_ALLOC;
//        REALLOC_FAILED(sizeof*mtx->indices * (mtx->n_elements));
        goto end;
    }
    mtx->indices = new_indices_ptr;
    mtx->capacity = mtx->n_entries;
end:
//    LEAVE_FUNCTION();
    return res;
}

jmtx_result jmtx_matrix_crs_set_row(jmtx_matrix_crs* mtx, uint32_t row, uint32_t n, const uint32_t* indices, const float* values)
{
#ifndef JMTX_NO_VERIFY_PARAMS
//    CALL_FUNCTION(matrix_crs_set_row);
    if (!mtx)
    {
//        REPORT_ERROR_MESSAGE("Matrix pointer was null");
//        LEAVE_FUNCTION();
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTX_TYPE_CRS)
    {
//        REPORT_ERROR_MESSAGE("Matrix was not compressed row sparse");
//        LEAVE_FUNCTION();
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (mtx->base.rows <= row)
    {
//        REPORT_ERROR_MESSAGE("Matrix has %u rows, but row %u was requested", mtx->base.rows, row);
//        LEAVE_FUNCTION();
        return JMTX_RESULT_INDEX_OUT_OF_BOUNDS;
    }
    if (mtx->base.cols < n)
    {
//        REPORT_ERROR_MESSAGE("Matrix has %u columns, but %u elements were specified to be set in row %u", mtx->base.cols, n, row);
//        LEAVE_FUNCTION();
        return JMTX_RESULT_INDEX_OUT_OF_BOUNDS;
    }
    if (!indices)
    {
//        REPORT_ERROR_MESSAGE("Indices pointer was null");
//        LEAVE_FUNCTION();
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!values)
    {
//        REPORT_ERROR_MESSAGE("Elements pointer was null");
//        LEAVE_FUNCTION();
        return JMTX_RESULT_NULL_PARAM;
    }
#endif
    jmtx_result res = JMTX_RESULT_SUCCESS;
    const uint32_t beginning_offset = row ? mtx->end_of_row_offsets[row - 1] : 0;
    const int32_t new_elements = (int32_t)n - (int32_t)(mtx->end_of_row_offsets[row] - beginning_offset);
    const uint32_t required_capacity = (uint32_t)((int32_t)mtx->n_entries + new_elements);
    if (mtx->capacity < required_capacity)
    {
        float* new_element_ptr = mtx->base.allocator_callbacks.realloc(mtx->base.allocator_callbacks.state, mtx->values, sizeof*(mtx->values) * (required_capacity + 1));
        if (!new_element_ptr)
        {
            res = JMTX_RESULT_BAD_ALLOC;
//            REALLOC_FAILED(sizeof*mtx->elements * (required_capacity + 1));
            goto end;
        }
        mtx->values = new_element_ptr;
        uint32_t* new_indices_ptr = mtx->base.allocator_callbacks.realloc(mtx->base.allocator_callbacks.state, mtx->indices, sizeof*(mtx->indices) * (required_capacity + 1));
        if (!new_indices_ptr)
        {
            res = JMTX_RESULT_BAD_ALLOC;
//            REALLOC_FAILED(sizeof*mtx->indices * (required_capacity + 1));
            goto end;
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
end:
//    LEAVE_FUNCTION();
    return res;
}

jmtx_result jmtx_matrix_crs_vector_multiply(const jmtx_matrix_crs* mtx, const float* restrict x, float* restrict y)
{
//    CALL_FUNCTION(matrix_crs_vector_multiply);
#ifndef JMTX_NO_VERIFY_PARAMS
    if (!mtx)
    {
//        REPORT_ERROR_MESSAGE("Matrix pointer was null");
//        LEAVE_FUNCTION();
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTX_TYPE_CRS)
    {
//        REPORT_ERROR_MESSAGE("Matrix was not compressed row sparse");
//        LEAVE_FUNCTION();
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (!x)
    {
//        REPORT_ERROR_MESSAGE("Vector x was null");
//        LEAVE_FUNCTION();
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!y)
    {
//        REPORT_ERROR_MESSAGE("Vector y was null");
//        LEAVE_FUNCTION();
        return JMTX_RESULT_NULL_PARAM;
    }
#endif


    jmtx_result res = JMTX_RESULT_SUCCESS;

    for (uint32_t i = 0; i < mtx->base.rows; ++i)
    {
        uint32_t* indices;
        float* values;
        const uint32_t row_entries = crs_get_row_entries(mtx, i, &indices, &values);
        float v = 0;
        for (uint32_t j = 0; j < row_entries; ++j)
        {
            const uint32_t k = indices[j];
            v += values[j] * x[k];
        }
        y[i] = v;
    }

    return res;
}

jmtx_result jmtx_matrix_crs_set_entry(jmtx_matrix_crs* mtx, uint32_t i, uint32_t j, float value)
{
//    CALL_FUNCTION(matrix_crs_set_element);
#ifndef JMTX_NO_VERIFY_PARAMS
    if (!mtx)
    {
//        REPORT_ERROR_MESSAGE("Matrix pointer was null");
//        LEAVE_FUNCTION();
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTX_TYPE_CRS)
    {
//        REPORT_ERROR_MESSAGE("Matrix was not compressed row sparse");
//        LEAVE_FUNCTION();
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (j >= mtx->base.cols)
    {
//        REPORT_ERROR_MESSAGE("Matrix has %u columns but column %u was requested", mtx->base.cols, j);
//        LEAVE_FUNCTION();
        return JMTX_RESULT_INDEX_OUT_OF_BOUNDS;
    }
    if (i >= mtx->base.rows)
    {
//        REPORT_ERROR_MESSAGE("Matrix has %u rows but row %u was requested", mtx->base.rows, i);
//        LEAVE_FUNCTION();
        return JMTX_RESULT_INDEX_OUT_OF_BOUNDS;
    }
#endif

    jmtx_result res;
    uint32_t* row_indices;
    float* row_values;
    const uint32_t n_row_elements = crs_get_row_entries(mtx, i, &row_indices, &row_values);
    //  Check if row has any elements
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

jmtx_result jmtx_matrix_crs_get_entry(const jmtx_matrix_crs* mtx, uint32_t i, uint32_t j, float* p_value)
{
//    CALL_FUNCTION(matrix_crs_get_element);
#ifndef JMTX_NO_VERIFY_PARAMS
    if (!mtx)
    {
//        REPORT_ERROR_MESSAGE("Matrix pointer was null");
//        LEAVE_FUNCTION();
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTX_TYPE_CRS)
    {
//        REPORT_ERROR_MESSAGE("Matrix was not compressed row sparse");
//        LEAVE_FUNCTION();
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (j >= mtx->base.cols)
    {
//        REPORT_ERROR_MESSAGE("Matrix has %u columns but column %u was requested", mtx->base.cols, j);
//        LEAVE_FUNCTION();
        return JMTX_RESULT_INDEX_OUT_OF_BOUNDS;
    }
    if (i >= mtx->base.rows)
    {
//        REPORT_ERROR_MESSAGE("Matrix has %u rows but row %u was requested", mtx->base.rows, i);
//        LEAVE_FUNCTION();
        return JMTX_RESULT_INDEX_OUT_OF_BOUNDS;
    }
#endif

    jmtx_result res = JMTX_RESULT_SUCCESS;
    uint32_t* row_indices;
    float* row_values;
    const uint32_t n_row_elements = crs_get_row_entries(mtx, i, &row_indices, &row_values);
    //  Check if row has any elements
    if (n_row_elements != 0)
    {
        //  Find first column entry less or equal to it
        const uint32_t possible = jmtx_internal_find_last_leq_value(n_row_elements, row_indices, j);
        if (row_indices[possible] == j)
        {
            *p_value = row_values[possible];
            return JMTX_RESULT_SUCCESS;
        }
    }
    *p_value = 0;
//    LEAVE_FUNCTION();
    return res;
}

jmtx_result jmtx_matrix_crs_get_row(const jmtx_matrix_crs* mtx, uint32_t row, uint32_t* n, uint32_t** p_indices, float** p_elements)
{
//    CALL_FUNCTION(matrix_crs_get_element);

#ifndef JMTX_NO_VERIFY_PARAMS
    if (!mtx)
    {
//        REPORT_ERROR_MESSAGE("Matrix pointer was null");
//        LEAVE_FUNCTION();
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTX_TYPE_CRS)
    {
//        REPORT_ERROR_MESSAGE("Matrix was not compressed row sparse");
//        LEAVE_FUNCTION();
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (row >= mtx->base.rows)
    {
//        REPORT_ERROR_MESSAGE("Matrix has %u rows but row %u was requested", mtx->base.rows, row);
//        LEAVE_FUNCTION();
        return JMTX_RESULT_INDEX_OUT_OF_BOUNDS;
    }
    if (!n)
    {
//        REPORT_ERROR_MESSAGE("Count pointer was null");
//        LEAVE_FUNCTION();
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!p_indices)
    {
//        REPORT_ERROR_MESSAGE("Index pointer was null");
//        LEAVE_FUNCTION();
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!p_elements)
    {
//        REPORT_ERROR_MESSAGE("Element pointer was null");
//        LEAVE_FUNCTION();
        return JMTX_RESULT_NULL_PARAM;
    }
#endif

    jmtx_result res = JMTX_RESULT_SUCCESS;
    uint32_t* row_indices;
    float* row_values;
    const uint32_t n_row_elements = crs_get_row_entries(mtx, row, &row_indices, &row_values);
    *p_indices = row_indices;
    *p_elements = row_values;
    *n = n_row_elements;
    //    LEAVE_FUNCTION();
    return res;
}

jmtx_result jmtx_matrix_crs_beef_check(const jmtx_matrix_crs* mtx, int* p_beef_status)
{
#ifndef JMTX_NO_VERIFY_PARAMS
    if (!mtx)
    {
//        REPORT_ERROR_MESSAGE("Matrix pointer was null");
//        LEAVE_FUNCTION();
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTX_TYPE_CRS)
    {
//        REPORT_ERROR_MESSAGE("Matrix was not compressed row sparse");
//        LEAVE_FUNCTION();
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (!p_beef_status)
    {
//        REPORT_ERROR_MESSAGE("Beef status pointer was null");
//        LEAVE_FUNCTION();
        return JMTX_RESULT_WRONG_TYPE;
    }
#endif

    const int beef_status = beef_check(mtx->values, mtx->n_entries) + beef_check((void*)mtx->indices, mtx->n_entries);
    *p_beef_status = (beef_status << 16) | 0x0000Beef;
//    LEAVE_FUNCTION();
    return 0;
}

jmtx_result jmtx_matrix_crs_apply_unary_fn(const jmtx_matrix_crs* mtx, int (*unary_fn)(uint32_t i, uint32_t j, float* p_value, void* param), void* param)
{
//    CALL_FUNCTION(matrix_crs_apply_unary_fn);
#ifndef JMTX_NO_VERIFY_PARAMS
    if (!mtx)
    {
//        REPORT_ERROR_MESSAGE("Matrix pointer was null");
//        LEAVE_FUNCTION();
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTX_TYPE_CRS)
    {
//        REPORT_ERROR_MESSAGE("Matrix was not compressed row sparse");
//        LEAVE_FUNCTION();
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (!unary_fn)
    {
//        REPORT_ERROR_MESSAGE("Unary function pointer was null");
//        LEAVE_FUNCTION();
        return JMTX_RESULT_NULL_PARAM;
    }
#endif

    for (uint32_t i = 0; i < mtx->base.rows; ++i)
    {
        float* p_elements;
        uint32_t* p_indices;
        const uint32_t n_in_row = crs_get_row_entries(mtx, i, &p_indices, &p_elements);
        for (uint32_t j = 0; j < n_in_row; ++j)
        {
            if ((unary_fn(i, p_indices[j], p_elements + j, param)))
            {
//                LEAVE_FUNCTION();
                return JMTX_RESULT_UNARY_RETURN;
            }
        }
    }
//    LEAVE_FUNCTION();
    return JMTX_RESULT_SUCCESS;
}

jmtx_result jmtx_matrix_crs_remove_zeros(jmtx_matrix_crs* mtx)
{
//    CALL_FUNCTION(matrix_crs_remove_zeros);
#ifndef JMTX_NO_VERIFY_PARAMS
    if (!mtx)
    {
//        REPORT_ERROR_MESSAGE("Matrix pointer was null");
//        LEAVE_FUNCTION();
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTX_TYPE_CRS)
    {
//        REPORT_ERROR_MESSAGE("Matrix was not compressed row sparse");
//        LEAVE_FUNCTION();
        return JMTX_RESULT_WRONG_TYPE;
    }
#endif
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


    //  Beef
    mtx->n_entries -= c;
#ifndef JMTX_NO_VERIFY_PARAMS
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

jmtx_result jmtx_matrix_crs_remove_bellow_magnitude(jmtx_matrix_crs* mtx, float v)
{
//    CALL_FUNCTION(matrix_crs_remove_zeros);
#ifndef JMTX_NO_VERIFY_PARAMS
    if (!mtx)
    {
//        REPORT_ERROR_MESSAGE("Matrix pointer was null");
//        LEAVE_FUNCTION();
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTX_TYPE_CRS)
    {
//        REPORT_ERROR_MESSAGE("Matrix was not compressed row sparse");
//        LEAVE_FUNCTION();
        return JMTX_RESULT_WRONG_TYPE;
    }
#endif
    //  Update offsets
    uint32_t p, c = 0, r = 0;
    while (mtx->end_of_row_offsets[r] == 0)
    {
        r += 1;
    }
    for (p = 0; p < mtx->n_entries; ++p)
    {
        assert(r < mtx->base.rows);
        if (is_bellow_magnitude(mtx->values[p], v))
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
    return 0;
}

jmtx_result jmtx_matrix_crs_entries_in_col(const jmtx_matrix_crs* mtx, uint32_t col, uint32_t* p_n)
{
#ifndef JMTX_NO_VERIFY_PARAMS
    if (!mtx)
    {
//        REPORT_ERROR_MESSAGE("Matrix pointer was null");
//        LEAVE_FUNCTION();
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTX_TYPE_CRS)
    {
//        REPORT_ERROR_MESSAGE("Matrix was not compressed row sparse");
//        LEAVE_FUNCTION();
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (col >= mtx->base.cols)
    {
//        REPORT_ERROR_MESSAGE("Matrix has %u columns but column %u was requested", mtx->base.cols, col);
//        LEAVE_FUNCTION();
        return JMTX_RESULT_INDEX_OUT_OF_BOUNDS;
    }
    if (!p_n)
    {
//        REPORT_ERROR_MESSAGE("Count pointer was null");
//        LEAVE_FUNCTION();
        return JMTX_RESULT_NULL_PARAM;
    }
#endif

    uint32_t element_count = 0;
    for (uint32_t row = 0; row < mtx->base.rows && (!row || (mtx->end_of_row_offsets[row - 1] != mtx->n_entries)); ++row)
    {
        uint32_t* row_indices;
        float* unused_row_values;
        const uint32_t n_row_elements = crs_get_row_entries(mtx, row, &row_indices, &unused_row_values);
        if (n_row_elements)
        {
            uint32_t current = jmtx_internal_find_last_leq_value(n_row_elements, row_indices, col);
            if (row_indices[current] == col)
            {
                element_count += 1;
            }
        }
    }
    *p_n = element_count;
    return JMTX_RESULT_SUCCESS;
}

jmtx_result jmtx_matrix_crs_get_col(
        const jmtx_matrix_crs* mtx, uint32_t col, uint32_t n, uint32_t* p_count, float* p_values, uint32_t* p_rows)
{
#ifndef JMTX_NO_VERIFY_PARAMS
    if (!mtx)
    {
//        REPORT_ERROR_MESSAGE("Matrix pointer was null");
//        LEAVE_FUNCTION();
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTX_TYPE_CRS)
    {
//        REPORT_ERROR_MESSAGE("Matrix was not compressed row sparse");
//        LEAVE_FUNCTION();
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (col >= mtx->base.cols)
    {
//        REPORT_ERROR_MESSAGE("Matrix has %u columns but column %u was requested", mtx->base.cols, col);
//        LEAVE_FUNCTION();
        return JMTX_RESULT_INDEX_OUT_OF_BOUNDS;
    }
    if (!p_values)
    {
//        REPORT_ERROR_MESSAGE("Element pointer was null");
//        LEAVE_FUNCTION();
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!p_rows)
    {
//        REPORT_ERROR_MESSAGE("Rows pointer was null");
//        LEAVE_FUNCTION();
        return JMTX_RESULT_NULL_PARAM;
    }
#endif
    uint32_t k = 0;
    for (uint32_t row = 0; k < n && row < mtx->base.rows && (!row || (mtx->end_of_row_offsets[row - 1] != mtx->n_entries)); ++row)
    {
        uint32_t* row_indices;
        float* unused_row_values;
        const uint32_t n_row_elements = crs_get_row_entries(mtx, row, &row_indices, &unused_row_values);
        if (n_row_elements)
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
    if (p_count)
    {
        *p_count = k;
    }
//    LEAVE_FUNCTION();
    return JMTX_RESULT_SUCCESS;
}

jmtx_result jmtx_matrix_crs_transpose(
        const jmtx_matrix_crs* mtx, jmtx_matrix_crs** p_out, const jmtx_allocator_callbacks* allocator_callbacks)
{
//    CALL_FUNCTION(matrix_crs_transpose);
#ifndef JMTX_NO_VERIFY_PARAMS
    if (!mtx)
    {
//        REPORT_ERROR_MESSAGE("Input matrix pointer was null");
//        LEAVE_FUNCTION();
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTX_TYPE_CRS)
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
    const uint32_t new_rows = mtx->base.cols;
    const uint32_t new_cols = mtx->base.rows;

    jmtx_matrix_crs* const out = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*out));
    if (!out)
    {
        return JMTX_RESULT_BAD_ALLOC;
    }

    uint32_t* const column_cum_counts = allocator_callbacks->alloc(allocator_callbacks->state, (new_rows) * sizeof(*column_cum_counts));
    if (!column_cum_counts)
    {
//        CALLOC_FAILED((new_rows + 1) * sizeof*column_cum_counts);
//        LEAVE_FUNCTION();
        allocator_callbacks->free(allocator_callbacks->state, out);
        return JMTX_RESULT_BAD_ALLOC;
    }
    uint32_t* const new_indices = allocator_callbacks->alloc(allocator_callbacks->state, (n_elements) * sizeof*new_indices);
    if (!new_indices)
    {
        allocator_callbacks->free(allocator_callbacks->state, out);
        allocator_callbacks->free(allocator_callbacks->state, column_cum_counts);
//        CALLOC_FAILED((n_elements + 1) * sizeof*new_indices);
//        LEAVE_FUNCTION();
        return JMTX_RESULT_BAD_ALLOC;
    }
    memset(new_indices, 0, (n_elements) * sizeof*new_indices);
    float* const new_values = allocator_callbacks->alloc(allocator_callbacks->state, (n_elements) * sizeof(*new_values));
    if (!new_values)
    {
        allocator_callbacks->free(allocator_callbacks->state, out);
        allocator_callbacks->free(allocator_callbacks->state, column_cum_counts);
        allocator_callbacks->free(allocator_callbacks->state, new_indices);
//        CALLOC_FAILED((n_elements + 1) * sizeof*new_elements);
//        LEAVE_FUNCTION();
        return JMTX_RESULT_BAD_ALLOC;
    }

    *column_cum_counts = 1;

    for (uint32_t j = 0, n, p = 0; j < mtx->base.cols; ++j)
    {
        //  This MUST NOT fail, since the parameters are all within the correct bounds, which is the only way it can fail
        jmtx_result res = jmtx_matrix_crs_entries_in_col(mtx, j, &n);
        assert(res == JMTX_RESULT_SUCCESS);

        res = jmtx_matrix_crs_get_col(mtx, j, n, NULL, new_values + p, new_indices + p);
        assert(res == JMTX_RESULT_SUCCESS);
        p += n;
        column_cum_counts[j] = (j != 0 ? column_cum_counts[j - 1] : 0) + n;
    }

    memcpy(out, mtx, sizeof*out);
    out->end_of_row_offsets = column_cum_counts;
    out->values = new_values;
    out->indices = new_indices;
    out->capacity = n_elements;
    out->base.rows = new_rows;
    out->base.cols = new_cols;
    out->base.allocator_callbacks = *allocator_callbacks;
    *p_out = out;
    //    LEAVE_FUNCTION();
    return JMTX_RESULT_SUCCESS;
}

jmtx_result jmtx_matrix_crs_copy(const jmtx_matrix_crs* mtx, jmtx_matrix_crs** p_out, const jmtx_allocator_callbacks* allocator_callbacks)
{
//    CALL_FUNCTION(matrix_crs_copy);
#ifndef JMTX_NO_VERIFY_PARAMS
    if (!mtx)
    {
//        REPORT_ERROR_MESSAGE("Input matrix pointer was null");
//        LEAVE_FUNCTION();
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTX_TYPE_CRS)
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
    jmtx_matrix_crs* const out = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*out));
    if (!out)
    {
        return JMTX_RESULT_BAD_ALLOC;
    }
    float* const elements = allocator_callbacks->alloc(allocator_callbacks->state, (mtx->n_entries) * sizeof (*elements));
    if (!elements)
    {
        allocator_callbacks->free(allocator_callbacks->state, out);
//        CALLOC_FAILED((1 + mtx->n_elements) * sizeof *elements);
//        LEAVE_FUNCTION();
        return JMTX_RESULT_BAD_ALLOC;
    }
    uint32_t* const indices = allocator_callbacks->alloc(allocator_callbacks->state, (mtx->n_entries) * sizeof *indices);
    if (!indices)
    {
        allocator_callbacks->free(allocator_callbacks->state, out);
        allocator_callbacks->free(allocator_callbacks->state, elements);
//        CALLOC_FAILED((1 + mtx->n_elements) * sizeof *indices);
//        LEAVE_FUNCTION();
        return JMTX_RESULT_BAD_ALLOC;
    }
    uint32_t* const cum_sum = allocator_callbacks->alloc(allocator_callbacks->state, (mtx->base.rows) * sizeof *cum_sum);
    if (!cum_sum)
    {
        allocator_callbacks->free(allocator_callbacks->state, out);
        allocator_callbacks->free(allocator_callbacks->state, indices);
        allocator_callbacks->free(allocator_callbacks->state, elements);
//        CALLOC_FAILED((1 + mtx->base.rows) * sizeof *cum_sum);
//        LEAVE_FUNCTION();
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
    out->base.allocator_callbacks = *allocator_callbacks;
    *p_out = out;
    return JMTX_RESULT_SUCCESS;
}

jmtx_result jmtx_matrix_crs_build_row(jmtx_matrix_crs* mtx, uint32_t row, uint32_t n, const uint32_t* indices, const float* values)
{
//    CALL_FUNCTION(matrix_crs_build_row);
//#ifndef JMTX_NO_VERIFY_PARAMS
//    if (!mtx)
//    {
////        REPORT_ERROR_MESSAGE("Matrix pointer was null");
////        LEAVE_FUNCTION();
//        return JMTX_RESULT_NULL_PARAM;
//    }
//    if (mtx->base.type != JMTX_TYPE_CRS)
//    {
////        REPORT_ERROR_MESSAGE("Matrix was not compressed row sparse");
////        LEAVE_FUNCTION();
//        return JMTX_RESULT_WRONG_TYPE;
//    }
//    if (mtx->base.rows <= row)
//    {
////        REPORT_ERROR_MESSAGE("Matrix has %u rows, but row %u was requested", mtx->base.rows, row);
////        LEAVE_FUNCTION();
//        return JMTX_RESULT_INDEX_OUT_OF_BOUNDS;
//    }
//    if (mtx->base.cols < n)
//    {
////        REPORT_ERROR_MESSAGE("Matrix has %u columns, but %u values were specified to be set in row %u", mtx->base.cols, n, row);
////        LEAVE_FUNCTION();
//        return JMTX_RESULT_INDEX_OUT_OF_BOUNDS;
//    }
//    if (!indices)
//    {
////        REPORT_ERROR_MESSAGE("Indices pointer was null");
////        LEAVE_FUNCTION();
//        return JMTX_RESULT_NULL_PARAM;
//    }
//    if (!values)
//    {
////        REPORT_ERROR_MESSAGE("Elements pointer was null");
////        LEAVE_FUNCTION();
//        return JMTX_RESULT_NULL_PARAM;
//    }
//#endif

    jmtx_result res = JMTX_RESULT_SUCCESS;
    const uint32_t required_capacity = (uint32_t)((int32_t)mtx->n_entries + (int32_t)n);
    if (mtx->capacity < required_capacity)
    {
        float* new_element_ptr = mtx->base.allocator_callbacks.realloc(mtx->base.allocator_callbacks.state, mtx->values, sizeof*mtx->values * (required_capacity + 1));
        if (!new_element_ptr)
        {
            res = JMTX_RESULT_BAD_ALLOC;
//            REALLOC_FAILED(sizeof*mtx->values * (required_capacity + 1));
            goto end;
        }
        mtx->values = new_element_ptr;
        uint32_t* new_indices_ptr = mtx->base.allocator_callbacks.realloc(mtx->base.allocator_callbacks.state, mtx->indices, sizeof*mtx->indices * (required_capacity + 1));
        if (!new_indices_ptr)
        {
            res = JMTX_RESULT_BAD_ALLOC;
//            REALLOC_FAILED(sizeof*mtx->indices * (required_capacity + 1));
            goto end;
        }
        mtx->indices = new_indices_ptr;
    }

    const uint32_t offset = row ? mtx->end_of_row_offsets[row - 1] : 0;
    memcpy(mtx->values + offset, values, sizeof*values * n);
    memcpy(mtx->indices + offset, indices, sizeof*indices * n);

    mtx->end_of_row_offsets[row] = n + offset;
    mtx->n_entries += n;
end:
//    LEAVE_FUNCTION();
    return res;
}

float jmtx_matrix_crs_vector_multiply_row_raw(const jmtx_matrix_crs* mtx, const float* x, uint32_t i)
{
    uint32_t* indices;
    float* values;
    const uint32_t n_row = crs_get_row_entries(mtx, i, &indices, &values);
    float v = 0;
    for (uint32_t j = 0; j < n_row; ++j)
    {
        v += values[j] * x[indices[j]];
    }
    return v;
}

jmtx_result jmtx_matrix_crs_vector_multiply_row(const jmtx_matrix_crs* mtx, const float* x, uint32_t i, float* p_r)
{
//    CALL_FUNCTION(matrix_crs_vector_multiply_row);
#ifndef JMTX_NO_VERIFY_PARAMS
    if (!mtx)
    {
//        REPORT_ERROR_MESSAGE("Matrix pointer was null");
//        LEAVE_FUNCTION();
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTX_TYPE_CRS)
    {
//        REPORT_ERROR_MESSAGE("Matrix was not compressed row sparse");
//        LEAVE_FUNCTION();
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (!x)
    {
//        REPORT_ERROR_MESSAGE("Vector x pointer was null");
//        LEAVE_FUNCTION();
        return JMTX_RESULT_NULL_PARAM;
    }
    if (i >= mtx->base.rows)
    {
//        REPORT_ERROR_MESSAGE("Matrix has %u rows but row %u was requested", mtx->base.rows, i);
//        LEAVE_FUNCTION();
        return JMTX_RESULT_INDEX_OUT_OF_BOUNDS;
    }
    if (!p_r)
    {
//        REPORT_ERROR_MESSAGE("Result pointer was null");
//        LEAVE_FUNCTION();
        return JMTX_RESULT_NULL_PARAM;
    }
#endif

    *p_r = jmtx_matrix_crs_vector_multiply_row_raw(mtx, x, i);
//    LEAVE_FUNCTION();
    return JMTX_RESULT_SUCCESS;
}

jmtx_result jmtx_matrix_crs_add_to_entry(jmtx_matrix_crs* mtx, uint32_t i, uint32_t j, float value)
{
//    CALL_FUNCTION(matrix_crs_set_element);
#ifndef JMTX_NO_VERIFY_PARAMS
    if (!mtx)
    {
//        REPORT_ERROR_MESSAGE("Matrix pointer was null");
//        LEAVE_FUNCTION();
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTX_TYPE_CRS)
    {
//        REPORT_ERROR_MESSAGE("Matrix was not compressed row sparse");
//        LEAVE_FUNCTION();
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (j >= mtx->base.cols)
    {
//        REPORT_ERROR_MESSAGE("Matrix has %u columns but column %u was requested", mtx->base.cols, j);
//        LEAVE_FUNCTION();
        return JMTX_RESULT_INDEX_OUT_OF_BOUNDS;
    }
    if (i >= mtx->base.rows)
    {
//        REPORT_ERROR_MESSAGE("Matrix has %u rows but row %u was requested", mtx->base.rows, i);
//        LEAVE_FUNCTION();
        return JMTX_RESULT_INDEX_OUT_OF_BOUNDS;
    }
#endif

    jmtx_result res;
    uint32_t* row_indices;
    float* row_values;
    const uint32_t n_row_elements = crs_get_row_entries(mtx, i, &row_indices, &row_values);
    //  Check if row has any elements
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

jmtx_result jmtx_matrix_crs_zero_all_entries(jmtx_matrix_crs* mtx)
{
#ifndef JMTX_NO_VERIFY_PARAMS
    if (!mtx)
    {
//        REPORT_ERROR_MESSAGE("Matrix pointer was null");
//        LEAVE_FUNCTION();
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTX_TYPE_CRS)
    {
//        REPORT_ERROR_MESSAGE("Matrix was not compressed row sparse");
//        LEAVE_FUNCTION();
        return JMTX_RESULT_WRONG_TYPE;
    }
#endif
    memset(mtx->values, 0, sizeof(*mtx->values) * mtx->n_entries);

    return JMTX_RESULT_SUCCESS;
}

jmtx_result jmtx_matrix_crs_set_all_entries(jmtx_matrix_crs* mtx, float x)
{
#ifndef JMTX_NO_VERIFY_PARAMS
    if (!mtx)
    {
//        REPORT_ERROR_MESSAGE("Matrix pointer was null");
//        LEAVE_FUNCTION();
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTX_TYPE_CRS)
    {
//        REPORT_ERROR_MESSAGE("Matrix was not compressed row sparse");
//        LEAVE_FUNCTION();
        return JMTX_RESULT_WRONG_TYPE;
    }
#endif
    for (float* ptr = mtx->values; ptr != mtx->values + mtx->n_entries; ++ptr)
    {
        *ptr = x;
    }

    return JMTX_RESULT_SUCCESS;

}

jmtx_result jmtx_matrix_crs_remove_row(jmtx_matrix_crs* mtx, uint32_t row)
{
#ifndef JMTX_NO_VERIFY_PARAMS
//    CALL_FUNCTION(jmtx_matrix_crs_remove_row);
    if (!mtx)
    {
//        REPORT_ERROR_MESSAGE("Matrix pointer was null");
//        LEAVE_FUNCTION();
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTX_TYPE_CRS)
    {
//        REPORT_ERROR_MESSAGE("Matrix was not compressed row sparse");
//        LEAVE_FUNCTION();
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (mtx->base.rows <= row)
    {
//        REPORT_ERROR_MESSAGE("Matrix has %u rows, but row %u was requested", mtx->base.rows, row);
//        LEAVE_FUNCTION();
        return JMTX_RESULT_INDEX_OUT_OF_BOUNDS;
    }
#endif
    jmtx_result res = JMTX_RESULT_SUCCESS;
    uint32_t* row_indices;
    float* row_values;
    const uint32_t removed_entry_count = crs_get_row_entries(mtx, row, &row_indices, &row_values);
    const uint32_t end_offset = mtx->end_of_row_offsets[row];
    //  Check if row is not empty
    if (removed_entry_count != 0)
    {
        //  Check if there are any entries following the ones being removed
        const uint32_t following_entries = mtx->n_entries - end_offset;
        if (following_entries != 0)
        {
            //  Move other elements to their new position
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
//    LEAVE_FUNCTION();
    return res;
}

jmtx_result jmtx_matrix_crs_remove_column(jmtx_matrix_crs* mtx, uint32_t col)
{
#ifndef JMTX_NO_VERIFY_PARAMS
//    CALL_FUNCTION(jmtx_matrix_crs_remove_row);
    if (!mtx)
    {
//        REPORT_ERROR_MESSAGE("Matrix pointer was null");
//        LEAVE_FUNCTION();
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTX_TYPE_CRS)
    {
//        REPORT_ERROR_MESSAGE("Matrix was not compressed row sparse");
//        LEAVE_FUNCTION();
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (mtx->base.cols <= col)
    {
//        REPORT_ERROR_MESSAGE("Matrix has %u rows, but row %u was requested", mtx->base.rows, row);
//        LEAVE_FUNCTION();
        return JMTX_RESULT_INDEX_OUT_OF_BOUNDS;
    }
#endif
    jmtx_result res = JMTX_RESULT_SUCCESS;
    //  i: tracks the current element to check
    //  j: how many elements to be removed were found
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
//    LEAVE_FUNCTION();
    return res;
}

jmtx_result jmtx_matrix_crs_clear(jmtx_matrix_crs* mtx)
{
#ifndef JMTX_NO_VERIFY_PARAMS
//    CALL_FUNCTION(jmtx_matrix_crs_remove_row);
    if (!mtx)
    {
//        REPORT_ERROR_MESSAGE("Matrix pointer was null");
//        LEAVE_FUNCTION();
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTX_TYPE_CRS)
    {
//        REPORT_ERROR_MESSAGE("Matrix was not compressed row sparse");
//        LEAVE_FUNCTION();
        return JMTX_RESULT_WRONG_TYPE;
    }
#endif
    mtx->n_entries = 0;
    memset(mtx->end_of_row_offsets, 0, sizeof(*mtx->end_of_row_offsets) * mtx->base.rows);
    return JMTX_RESULT_SUCCESS;
}

jmtx_result jmtx_matrix_crs_join_vertically(
        jmtx_matrix_crs** output, const jmtx_allocator_callbacks* allocators, unsigned int k,
        const jmtx_matrix_crs** matrix_list)
{
    const uint32_t col_count = matrix_list[0]->base.cols;
    uint32_t n_rows = matrix_list[0]->base.rows;

    uint32_t element_count = matrix_list[0]->n_entries;
    //  First one is already accounted for, so no need to start at 0
    for (unsigned i = 1; i < k; ++i)
    {
        const jmtx_matrix_crs* const e = matrix_list[i];
        if (e->base.cols != col_count)
        {
            return JMTX_RESULT_DIMS_MISMATCH;
        }
        n_rows += e->base.rows;
        element_count += e->n_entries;
    }

    jmtx_matrix_crs* out;
    jmtx_result res = jmtx_matrix_crs_new(&out, col_count, n_rows, element_count, allocators);
    if (res != JMTX_RESULT_SUCCESS)
    {
        return res;
    }

    uint32_t pos = 0;
    for (unsigned i = 0; i < k; ++i)
    {
        const jmtx_matrix_crs* const mtx = matrix_list[i];
        unsigned j;
        for (j = 0; j < mtx->base.rows; ++j)
        {
            uint32_t* indices;
            float* values;
            const uint32_t n = crs_get_row_entries(mtx, j, &indices, &values);
            //  All memory for this should be allocated in advance, so no need to check the return value
            const jmtx_result r1 = jmtx_matrix_crs_build_row(out, pos + j, n, indices, values);
            assert(r1 == JMTX_RESULT_SUCCESS); (void)r1;
        }
        pos += j;
    }
    assert(pos == n_rows);

    *output = out;
    return JMTX_RESULT_SUCCESS;
}
