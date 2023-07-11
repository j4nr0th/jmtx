//
// Created by jan on 15.6.2022.
//

#include <assert.h>
#include "sparse_column_compressed.h"
//
// Created by jan on 13.6.2022.
//



#define DEFAULT_RESERVED_ELEMENTS 64

#define FILL_VALUE 0xDEADBEEF

static void beef_it_up(jmtx_scalar_t* ptr, size_t elements)
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

static int beef_check(const jmtx_scalar_t* ptr, size_t elements)
{
    const uint32_t* const buffer = (const uint32_t*)ptr;
    int beef_count = 0;
    for (uint32_t i = 0; i < elements; ++i)
    {
        beef_count += (buffer[i] == 0xDEADBEEF);
    }
    return (beef_count << 16) | 0x0000BEEF;
}

jmtx_result matrix_ccs_new(
        jmtx_matrix_ccs* mtx, uint32_t columns, uint32_t rows, uint32_t reserved_elements,
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
    if (reserved_elements > columns * rows)
    {
//        //REPORT_ERROR_MESSAGE("Number of reserved elements (%u) exceeds product of columns (%u) by rows (%u)", reserved_elements, rows, columns);
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

    if (reserved_elements == 0)
    {
        reserved_elements = DEFAULT_RESERVED_ELEMENTS;
        reserved_elements = reserved_elements < columns * rows ? reserved_elements : columns * rows;
    }
    jmtx_result res = JMTX_RESULT_SUCCESS;
    uint32_t* elements_per_col = NULL;
    uint32_t* indices = NULL;
    jmtx_scalar_t* p_elements = allocator_callbacks->alloc(allocator_callbacks->state, (1 + reserved_elements) * sizeof(*p_elements));
    if (!p_elements)
    {
        res = JMTX_RESULT_BAD_ALLOC;
//        CALLOC_FAILED((1 + reserved_elements) * sizeof*p_elements);
        goto fail1;
    }

    indices = allocator_callbacks->alloc(allocator_callbacks->state, (1 + reserved_elements) * sizeof*indices);
    if (!indices)
    {
        res = JMTX_RESULT_BAD_ALLOC;
//        CALLOC_FAILED((1 + reserved_elements) * sizeof*indices);
        goto fail2;
    }

    elements_per_col = allocator_callbacks->alloc(allocator_callbacks->state, (columns + 1) * sizeof*elements_per_col);
    if (!elements_per_col)
    {
        res = JMTX_RESULT_BAD_ALLOC;
//        CALLOC_FAILED((columns + 1) * sizeof*elements_per_col);
        goto fail3;
    }
    beef_it_up(p_elements + 1, reserved_elements);
    static_assert(sizeof(jmtx_scalar_t) == sizeof(uint32_t), "element and index sizes must be the same");
    beef_it_up((jmtx_scalar_t*)indices + 1, reserved_elements);
    for (uint32_t i = 0; i < columns + 1; ++i)
    {
        elements_per_col[i] = 1;
    }
    memset(mtx, 0, sizeof*mtx);
    mtx->base.cols = columns;
    mtx->base.type = JMTX_TYPE_CCS;
    mtx->base.rows = rows;
    mtx->base.allocator_callbacks = *allocator_callbacks;
    mtx->indices = indices;
    mtx->elements = p_elements;
    mtx->capacity = reserved_elements;
    mtx->n_elements = 0;
    mtx->elements_before = elements_per_col;
//    //LEAVE_FUNCTION();
    return res;
    fail3: allocator_callbacks->free(allocator_callbacks->state, indices);
    fail2: allocator_callbacks->free(allocator_callbacks->state, elements_per_col);
    fail1: allocator_callbacks->free(allocator_callbacks->state, p_elements);
//    //LEAVE_FUNCTION();
    return res;
}

jmtx_result matrix_ccs_destroy(jmtx_matrix_ccs* mtx)
{
    if (!mtx)
    {
//        //REPORT_ERROR_MESSAGE("Matrix pointer was null");
//        //LEAVE_FUNCTION();
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTX_TYPE_CCS)
    {
//        //REPORT_ERROR_MESSAGE("Matrix was not compressed column sparse");
//        //LEAVE_FUNCTION();
        return JMTX_RESULT_WRONG_TYPE;
    }

    jmtx_result res = JMTX_RESULT_SUCCESS;
    free(mtx->indices);
    free(mtx->elements_before);
    free(mtx->elements);
    //LEAVE_FUNCTION();
    return res;
}

jmtx_result matrix_ccs_shrink(jmtx_matrix_ccs* mtx)
{
    if (!mtx)
    {
        //REPORT_ERROR_MESSAGE("Matrix pointer was null");
        //LEAVE_FUNCTION();
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTX_TYPE_CCS)
    {
        //REPORT_ERROR_MESSAGE("Matrix was not compressed column sparse");
        //LEAVE_FUNCTION();
        return JMTX_RESULT_WRONG_TYPE;
    }
    jmtx_result res = JMTX_RESULT_SUCCESS;

    if (mtx->n_elements == mtx->capacity)
    {
        //LEAVE_FUNCTION();
        return JMTX_RESULT_SUCCESS;
    }
    mtx->capacity = mtx->n_elements;
    jmtx_scalar_t* element_new_ptr = mtx->base.allocator_callbacks.realloc(mtx->base.allocator_callbacks.state, mtx->elements, sizeof*mtx->elements * (mtx->n_elements + 1));
    if (!element_new_ptr)
    {
        res = JMTX_RESULT_BAD_ALLOC;
        //REALLOC_FAILED(sizeof*mtx->elements * (mtx->n_elements + 1));
        goto end;
    }
    mtx->elements = element_new_ptr;
    uint32_t* new_indices_ptr = mtx->base.allocator_callbacks.realloc(mtx->base.allocator_callbacks.state, mtx->indices, sizeof*mtx->indices * (mtx->n_elements + 1));
    if (!new_indices_ptr)
    {
        res = JMTX_RESULT_BAD_ALLOC;
        //REALLOC_FAILED(sizeof*mtx->indices * (mtx->n_elements + 1));
        goto end;
    }
    mtx->indices = new_indices_ptr;

    end:
    //LEAVE_FUNCTION();
    return res;
}

jmtx_result matrix_ccs_set_col(jmtx_matrix_ccs* mtx, uint32_t col, uint32_t n, const uint32_t* indices, const jmtx_scalar_t* elements)
{
    if (!mtx)
    {
        //REPORT_ERROR_MESSAGE("Matrix pointer was null");
        //LEAVE_FUNCTION();
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTX_TYPE_CCS)
    {
        //REPORT_ERROR_MESSAGE("Matrix was not compressed column sparse");
        //LEAVE_FUNCTION();
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (mtx->base.cols <= col)
    {
        //REPORT_ERROR_MESSAGE("Matrix has %u columns, but column %u was requested", mtx->base.cols, col);
        //LEAVE_FUNCTION();
        return JMTX_RESULT_INDEX_OUT_OF_BOUNDS;
    }
    if (mtx->base.rows < n)
    {
        //REPORT_ERROR_MESSAGE("Matrix has %u rows, but %u elements were specified to be set in column %u", mtx->base.rows, n, col);
        //LEAVE_FUNCTION();
        return JMTX_RESULT_INDEX_OUT_OF_BOUNDS;
    }
    if (!indices)
    {
        //REPORT_ERROR_MESSAGE("Indices pointer was null");
        //LEAVE_FUNCTION();
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!elements)
    {
        //REPORT_ERROR_MESSAGE("Elements pointer was null");
        //LEAVE_FUNCTION();
        return JMTX_RESULT_NULL_PARAM;
    }

    jmtx_result res = 0;
    const int32_t new_elements = (int32_t)n - (int32_t)(mtx->elements_before[col + 1] - mtx->elements_before[col]);
    const uint32_t required_capacity = (uint32_t)((int32_t)mtx->n_elements + new_elements);
    if (mtx->capacity < required_capacity)
    {
        jmtx_scalar_t* new_element_ptr = mtx->base.allocator_callbacks.realloc(mtx->base.allocator_callbacks.state, mtx->elements, sizeof*mtx->elements * (required_capacity + 1));
        if (!new_element_ptr)
        {
            res = JMTX_RESULT_BAD_ALLOC;
            //REALLOC_FAILED(sizeof*mtx->elements * (required_capacity + 1));
            goto end;
        }
        mtx->elements = new_element_ptr;
        uint32_t* new_indices_ptr = mtx->base.allocator_callbacks.realloc(mtx->base.allocator_callbacks.state, mtx->indices, sizeof*mtx->indices * (required_capacity + 1));
        if (!new_indices_ptr)
        {
            res = JMTX_RESULT_BAD_ALLOC;
            //REALLOC_FAILED(sizeof*mtx->indices * (required_capacity + 1));
            goto end;
        }
        mtx->indices = new_indices_ptr;
    }

    const uint32_t elements_after = mtx->elements_before[mtx->base.cols] - mtx->elements_before[col + 1];
    if (elements_after)
    {
        memmove(mtx->elements + mtx->elements_before[col + 1] + new_elements, mtx->elements + mtx->elements_before[col + 1],
                sizeof*mtx->elements * (elements_after));
        memmove(mtx->indices + mtx->elements_before[col + 1] + new_elements, mtx->indices + mtx->elements_before[col + 1],
                sizeof*mtx->indices * (elements_after));
    }
    memcpy(mtx->elements + mtx->elements_before[col], elements, sizeof*elements * n);
    memcpy(mtx->indices + mtx->elements_before[col], indices, sizeof*indices * n);

    for (uint32_t i = col; i < mtx->base.cols; ++i)
    {
        mtx->elements_before[i + 1] += new_elements;
    }
    mtx->n_elements += new_elements;
end:
    //LEAVE_FUNCTION();
    return res;
}

jmtx_result matrix_ccs_vector_multiply(const jmtx_matrix_ccs* mtx, const jmtx_scalar_t* restrict x, jmtx_scalar_t* restrict y)
{

    if (!mtx)
    {
        //REPORT_ERROR_MESSAGE("Matrix pointer was null");
        //LEAVE_FUNCTION();
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTX_TYPE_CCS)
    {
        //REPORT_ERROR_MESSAGE("Matrix was not compressed column sparse");
        //LEAVE_FUNCTION();
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (!x)
    {
        //REPORT_ERROR_MESSAGE("Vector x was null");
        //LEAVE_FUNCTION();
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!y)
    {
        //REPORT_ERROR_MESSAGE("Vector y was null");
        //LEAVE_FUNCTION();
        return JMTX_RESULT_NULL_PARAM;
    }
    jmtx_result res = 0;

    for (uint32_t i = 0; i < mtx->base.cols; ++i)
    {
        const uint32_t* indices = mtx->indices + mtx->elements_before[i];
        const jmtx_scalar_t* col_ptr = mtx->elements + mtx->elements_before[i];
        jmtx_scalar_t v = 0;
        const uint32_t n_elements = mtx->elements_before[i + 1] - mtx->elements_before[i];
        for (uint32_t j = 0; j < n_elements; ++j)
        {
            const uint32_t k = indices[j];
            v += col_ptr[j] * x[k];
        }
        y[i] = v;
    }
    //LEAVE_FUNCTION();
    return res;
}

jmtx_result matrix_ccs_set_element(jmtx_matrix_ccs* mtx, uint32_t i, uint32_t j, jmtx_scalar_t x)
{

    if (!mtx)
    {
        //REPORT_ERROR_MESSAGE("Matrix pointer was null");
        //LEAVE_FUNCTION();
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTX_TYPE_CCS)
    {
        //REPORT_ERROR_MESSAGE("Matrix was not compressed column sparse");
        //LEAVE_FUNCTION();
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (j >= mtx->base.cols)
    {
        //REPORT_ERROR_MESSAGE("Matrix has %u columns but column %u was requested", mtx->base.cols, j);
        //LEAVE_FUNCTION();
        return JMTX_RESULT_INDEX_OUT_OF_BOUNDS;
    }
    if (i >= mtx->base.rows)
    {
        //REPORT_ERROR_MESSAGE("Matrix has %u rows but row %u was requested", mtx->base.rows, i);
        //LEAVE_FUNCTION();
        return JMTX_RESULT_INDEX_OUT_OF_BOUNDS;
    }

    jmtx_result res = JMTX_RESULT_SUCCESS;

    const uint32_t n_col_elements = mtx->elements_before[j + 1] - mtx->elements_before[j];
    const uint32_t* const col_indices = mtx->indices + mtx->elements_before[j];
    uint32_t current = 0;
    for (
            uint32_t size = n_col_elements, step = n_col_elements / 2;
            step != 0;
            step = size / 2
            )
    {
        if (col_indices[current + step] < j)
        {
            current += step;
            size -= step;
        }
        else if (col_indices[current + step] == j)
        {
            current += step;
            break;
        }
        else
        {
            size = step;
        }
    }

    if (n_col_elements && col_indices[current] == j)
    {
        *(mtx->elements + mtx->elements_before[j] + current) = x;
    }
    else
    {
        if (mtx->capacity == mtx->n_elements)
        {
            const uint32_t new_capacity = mtx->capacity + DEFAULT_RESERVED_ELEMENTS;
            jmtx_scalar_t* const new_element_ptr = mtx->base.allocator_callbacks.realloc(mtx->base.allocator_callbacks.state, mtx->elements, sizeof*mtx->elements * (new_capacity + 1));
            if (!new_element_ptr)
            {
                res = JMTX_RESULT_BAD_ALLOC;
                //REALLOC_FAILED(sizeof*mtx->elements * (new_capacity + 1));
                goto end;
            }
            mtx->elements = new_element_ptr;
            beef_it_up(new_element_ptr + 1 + mtx->n_elements, new_capacity - mtx->capacity);
            uint32_t* const new_index_ptr = mtx->base.allocator_callbacks.realloc(mtx->base.allocator_callbacks.state, mtx->indices, sizeof*mtx->indices * (new_capacity + 1));
            if (!new_index_ptr)
            {
                res = JMTX_RESULT_BAD_ALLOC;
                //REALLOC_FAILED(sizeof*mtx->indices * (new_capacity + 1));
                goto end;
            }
            mtx->indices = new_index_ptr;
            beef_it_up((jmtx_scalar_t *)(new_index_ptr + 1 + mtx->n_elements), new_capacity - mtx->capacity);
            mtx->capacity = new_capacity;
        }
        if (mtx->elements_before[mtx->base.cols] - mtx->elements_before[j + 1] != 0)
        {
            current += (n_col_elements != 0);
            memmove(mtx->elements + mtx->elements_before[j] + current + 1, mtx->elements + mtx->elements_before[j] + current, (mtx->elements_before[mtx->base.cols] - mtx->elements_before[j + 1] + n_col_elements - current) * sizeof(*mtx->elements));
            memmove(mtx->indices + mtx->elements_before[j] + current + 1, mtx->indices + mtx->elements_before[j] + current, (mtx->elements_before[mtx->base.cols] - mtx->elements_before[j + 1] + n_col_elements - current) * sizeof(*mtx->indices));
        }
        *(mtx->elements + mtx->elements_before[j] + current) = x;
        *(mtx->indices + mtx->elements_before[j] + current) = j;
        mtx->n_elements += 1;
        for (uint32_t k = j; k < mtx->base.cols; ++k)
        {
            mtx->elements_before[k + 1] += 1;
        }
    }
end:
    //LEAVE_FUNCTION();
    return res;
}

jmtx_result matrix_ccs_get_element(const jmtx_matrix_ccs* mtx, uint32_t i, uint32_t j, jmtx_scalar_t* x)
{

    if (!mtx)
    {
        //REPORT_ERROR_MESSAGE("Matrix pointer was null");
        //LEAVE_FUNCTION();
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTX_TYPE_CCS)
    {
        //REPORT_ERROR_MESSAGE("Matrix was not compressed column sparse");
        //LEAVE_FUNCTION();
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (j >= mtx->base.cols)
    {
        //REPORT_ERROR_MESSAGE("Matrix has %u columns but column %u was requested", mtx->base.cols, j);
        //LEAVE_FUNCTION();
        return JMTX_RESULT_INDEX_OUT_OF_BOUNDS;
    }
    if (i >= mtx->base.rows)
    {
        //REPORT_ERROR_MESSAGE("Matrix has %u rows but row %u was requested", mtx->base.rows, i);
        //LEAVE_FUNCTION();
        return JMTX_RESULT_INDEX_OUT_OF_BOUNDS;
    }

    jmtx_result res = JMTX_RESULT_SUCCESS;
    const uint32_t n_col_elements = mtx->elements_before[j + 1] - mtx->elements_before[j];
    const uint32_t* const col_indices = mtx->indices + mtx->elements_before[j];
    uint32_t current = 0;
    for (
            uint32_t size = n_col_elements, step = n_col_elements / 2;
            step != 0;
            step = size / 2
            )
    {
        if (col_indices[current + step] < i)
        {
            current += step;
            size -= step;
        }
        else if (col_indices[current + step] == i)
        {
            current += step;
            break;
        }
        else
        {
            size = step;
        }
    }

    if (n_col_elements && col_indices[current] == i)
    {
        *x = *(mtx->elements + mtx->elements_before[j] + current);
    }
    else
    {
        *x = (jmtx_scalar_t)0.0;
    }

    //LEAVE_FUNCTION();
    return res;
}

jmtx_result matrix_ccs_get_col(const jmtx_matrix_ccs* mtx, uint32_t col, uint32_t* n, uint32_t** p_indices, jmtx_scalar_t** p_elements)
{

    if (!mtx)
    {
        //REPORT_ERROR_MESSAGE("Matrix pointer was null");
        //LEAVE_FUNCTION();
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTX_TYPE_CCS)
    {
        //REPORT_ERROR_MESSAGE("Matrix was not compressed column sparse");
        //LEAVE_FUNCTION();
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (col >= mtx->base.cols)
    {
        //REPORT_ERROR_MESSAGE("Matrix has %u columns but column %u was requested", mtx->base.cols, col);
        //LEAVE_FUNCTION();
        return JMTX_RESULT_INDEX_OUT_OF_BOUNDS;
    }
    if (!n)
    {
        //REPORT_ERROR_MESSAGE("Count pointer was null");
        //LEAVE_FUNCTION();
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!p_indices)
    {
        //REPORT_ERROR_MESSAGE("Index pointer was null");
        //LEAVE_FUNCTION();
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!p_elements)
    {
        //REPORT_ERROR_MESSAGE("Element pointer was null");
        //LEAVE_FUNCTION();
        return JMTX_RESULT_NULL_PARAM;
    }

    jmtx_result res = 0;
    const uint32_t n_col = mtx->elements_before[col + 1] - mtx->elements_before[col];
    *n = n_col;
    if (n_col)
    {
        *p_indices = mtx->indices + mtx->elements_before[col];
        *p_elements = mtx->elements + mtx->elements_before[col];
    }

    //LEAVE_FUNCTION();
    return res;
}

jmtx_result matrix_ccs_beef_check(const jmtx_matrix_ccs* mtx, int* p_beef_status)
{

    if (!mtx)
    {
        //REPORT_ERROR_MESSAGE("Matrix pointer was null");
        //LEAVE_FUNCTION();
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTX_TYPE_CCS)
    {
        //REPORT_ERROR_MESSAGE("Matrix was not compressed column sparse");
        //LEAVE_FUNCTION();
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (!p_beef_status)
    {
        //REPORT_ERROR_MESSAGE("Beef status pointer was null");
        //LEAVE_FUNCTION();
        return JMTX_RESULT_WRONG_TYPE;
    }

    const int beef_status = beef_check(mtx->elements + 1, mtx->n_elements);
    *p_beef_status = beef_status;
    //LEAVE_FUNCTION();
    return 0;
}

jmtx_result matrix_ccs_apply_unary_fn(const jmtx_matrix_ccs* mtx, int (*unary_fn)(uint32_t i, uint32_t j, jmtx_scalar_t* p_element, void* param), void* param)
{

    if (!mtx)
    {
        //REPORT_ERROR_MESSAGE("Matrix pointer was null");
        //LEAVE_FUNCTION();
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTX_TYPE_CCS)
    {
        //REPORT_ERROR_MESSAGE("Matrix was not compressed column sparse");
        //LEAVE_FUNCTION();
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (!unary_fn)
    {
        //REPORT_ERROR_MESSAGE("Unary function pointer was null");
        //LEAVE_FUNCTION();
        return JMTX_RESULT_NULL_PARAM;
    }

    for (uint32_t j = 0; j < mtx->base.cols; ++j)
    {
        const uint32_t n_in_col = mtx->elements_before[j + 1] - mtx->elements_before[j];
        jmtx_scalar_t* const p_elements = mtx->elements + mtx->elements_before[j];
        const uint32_t* const p_indices = mtx->indices + mtx->elements_before[j];
        for (uint32_t i = 0; i < n_in_col; ++i)
        {
            int res;
            if ((res = unary_fn(p_indices[i], j, p_elements + j, param)))
            {
                //LEAVE_FUNCTION();
                return res;
            }
        }
    }
    //LEAVE_FUNCTION();
    return JMTX_RESULT_SUCCESS;
}

jmtx_result matrix_ccs_remove_zeros(jmtx_matrix_ccs* mtx)
{

    if (!mtx)
    {
        //REPORT_ERROR_MESSAGE("Matrix pointer was null");
        //LEAVE_FUNCTION();
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTX_TYPE_CCS)
    {
        //REPORT_ERROR_MESSAGE("Matrix was not compressed column sparse");
        //LEAVE_FUNCTION();
        return JMTX_RESULT_WRONG_TYPE;
    }

    uint32_t zero_count, k, l;
    for (zero_count = 0, k = 0; k < mtx->n_elements; ++k)
    {
        zero_count += (mtx->elements[k + 1] == (jmtx_scalar_t)0.0);
    }
    if (!zero_count)
    {
        //LEAVE_FUNCTION();
        return JMTX_RESULT_SUCCESS;
    }

    uint32_t* const zero_indices = mtx->base.allocator_callbacks.alloc(mtx->base.allocator_callbacks.state, zero_count * sizeof*zero_indices);
    if (!zero_indices)
    {
//        CALLOC_FAILED(zero_count * sizeof*zero_indices);
        //LEAVE_FUNCTION();
        return JMTX_RESULT_BAD_ALLOC;
    }

    for (k = 0, l = 0; l < zero_count && k < mtx->n_elements; ++k)
    {
        if (mtx->elements[k + 1] == (jmtx_scalar_t)0.0)
        {
            zero_indices[l++] = k;
        }
    }

    const uint32_t original_zero_count = zero_count;
    while (zero_count)
    {
        uint32_t consecutive_zeros;
        for (consecutive_zeros = 1; consecutive_zeros < zero_count; ++consecutive_zeros)
        {
            if (zero_indices[zero_count - consecutive_zeros] != zero_indices[zero_count - consecutive_zeros - 1])
                break;
        }

        memmove(mtx->elements + 1 + zero_indices[zero_count - consecutive_zeros],
                mtx->elements + 1 + zero_indices[zero_count - 1] + 1,
                (mtx->n_elements - zero_indices[zero_count - 1] - 1) * sizeof*mtx->elements);
        memmove(mtx->indices + 1 + zero_indices[zero_count - consecutive_zeros],
                mtx->indices + 1 + zero_indices[zero_count - 1] + 1,
                (mtx->n_elements - zero_indices[zero_count - 1] - 1) * sizeof*mtx->indices);

        mtx->n_elements -= consecutive_zeros;
        zero_count -= consecutive_zeros;
    }

    for (k = 0, l = 0, zero_count = 0; k < mtx->base.cols + 1; ++k)
    {
        while (l < original_zero_count && zero_indices[l] + 1 < mtx->elements_before[k])
        {
            zero_count += 1;
            l += 1;
        }
        mtx->elements_before[k] -= zero_count;
    }

    //  Beef
    beef_it_up(mtx->elements + 1 + mtx->n_elements, original_zero_count);
    static_assert(sizeof(jmtx_scalar_t) == sizeof(uint32_t), "Size of index and scalar must be the same for beef");
    beef_it_up((jmtx_scalar_t*)(mtx->indices + 1 + mtx->n_elements), original_zero_count);


    mtx->base.allocator_callbacks.free(mtx->base.allocator_callbacks.state, zero_indices);
    //LEAVE_FUNCTION();
    return 0;
}

jmtx_result matrix_ccs_remove_bellow(jmtx_matrix_ccs* mtx, jmtx_scalar_t v)
{

    if (!mtx)
    {
        //REPORT_ERROR_MESSAGE("Matrix pointer was null");
        //LEAVE_FUNCTION();
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTX_TYPE_CCS)
    {
        //REPORT_ERROR_MESSAGE("Matrix was not compressed column sparse");
        //LEAVE_FUNCTION();
        return JMTX_RESULT_WRONG_TYPE;
    }

    //  AND-ing with 0x7FFFFFFF should remove the sign bit, allowing to use fast comparison of float's absolute values
    const uint32_t int_abs_v = ((*(uint32_t*)&v) & 0x7FFFFFFF);
    const jmtx_scalar_t abs_v = *(jmtx_scalar_t*)&int_abs_v;
    static_assert(sizeof(jmtx_scalar_t) == sizeof(uint32_t), "Size of scalar and uint32_t must be the same for this to work");
    uint32_t zero_count, k, l;
    for (zero_count = 0, k = 0; k < mtx->n_elements; ++k)
    {
        const uint32_t element_abs = ((uint32_t*)mtx->elements)[k + 1] & 0x7FFFFFFF;
        zero_count += (*(jmtx_scalar_t*)&element_abs < abs_v);
    }
    if (!zero_count)
    {
        //LEAVE_FUNCTION();
        return JMTX_RESULT_SUCCESS;
    }

    uint32_t* const zero_indices = mtx->base.allocator_callbacks.alloc(mtx->base.allocator_callbacks.state, zero_count * sizeof*zero_indices);
    if (!zero_indices)
    {
//        CALLOC_FAILED(zero_count * sizeof*zero_indices);
        //LEAVE_FUNCTION();
        return JMTX_RESULT_BAD_ALLOC;
    }

    for (k = 0, l = 0; l < zero_count && k < mtx->n_elements; ++k)
    {
        const uint32_t element_abs = ((uint32_t*)mtx->elements)[k + 1] & 0x7FFFFFFF;
        if (*(jmtx_scalar_t*)&element_abs < abs_v)
        {
            zero_indices[l++] = k;
        }
    }

    const uint32_t original_zero_count = zero_count;
    while (zero_count)
    {
        uint32_t consecutive_zeros;
        for (consecutive_zeros = 1; consecutive_zeros < zero_count; ++consecutive_zeros)
        {
            if (zero_indices[zero_count - consecutive_zeros] != zero_indices[zero_count - consecutive_zeros - 1])
                break;
        }

        memmove(mtx->elements + 1 + zero_indices[zero_count - consecutive_zeros],
                mtx->elements + 1 + zero_indices[zero_count - 1] + 1,
                (mtx->n_elements - zero_indices[zero_count - 1] - 1) * sizeof*mtx->elements);
        memmove(mtx->indices + 1 + zero_indices[zero_count - consecutive_zeros],
                mtx->indices + 1 + zero_indices[zero_count - 1] + 1,
                (mtx->n_elements - zero_indices[zero_count - 1] - 1) * sizeof*mtx->indices);

        mtx->n_elements -= consecutive_zeros;
        zero_count -= consecutive_zeros;
    }

    for (k = 0, l = 0, zero_count = 0; k < mtx->base.cols + 1; ++k)
    {
        while (l < original_zero_count && zero_indices[l] + 1 < mtx->elements_before[k])
        {
            zero_count += 1;
            l += 1;
        }
        mtx->elements_before[k] -= zero_count;
    }

    //  Beef
    beef_it_up(mtx->elements + 1 + mtx->n_elements, original_zero_count);
    static_assert(sizeof(jmtx_scalar_t) == sizeof(uint32_t), "Size of index and scalar must be the same for beef");
    beef_it_up((jmtx_scalar_t*)(mtx->indices + 1 + mtx->n_elements), original_zero_count);


    mtx->base.allocator_callbacks.free(mtx->base.allocator_callbacks.state, zero_indices);
    //LEAVE_FUNCTION();
    return 0;
}

jmtx_result matrix_ccs_elements_in_row(const jmtx_matrix_ccs* mtx, uint32_t row, uint32_t* p_n)
{

    if (!mtx)
    {
        //REPORT_ERROR_MESSAGE("Matrix pointer was null");
        //LEAVE_FUNCTION();
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTX_TYPE_CCS)
    {
        //REPORT_ERROR_MESSAGE("Matrix was not compressed column sparse");
        //LEAVE_FUNCTION();
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (row >= mtx->base.rows)
    {
        //REPORT_ERROR_MESSAGE("Matrix has %u rows but row %u was requested", mtx->base.rows, row);
        //LEAVE_FUNCTION();
        return JMTX_RESULT_INDEX_OUT_OF_BOUNDS;
    }
    if (!p_n)
    {
        //REPORT_ERROR_MESSAGE("Count pointer was null");
        //LEAVE_FUNCTION();
        return JMTX_RESULT_NULL_PARAM;
    }

    uint32_t element_count = 0;
    for (uint32_t col = 0; col < mtx->base.cols && mtx->elements_before[col] != mtx->elements_before[mtx->base.cols]; ++col)
    {
        uint32_t current = 0;
        const uint32_t n_col_elements = mtx->elements_before[col + 1] - mtx->elements_before[col];
        const uint32_t* col_indices = mtx->indices + mtx->elements_before[col];
        for (
                uint32_t size = n_col_elements, step = n_col_elements / 2;
                step != 0;
                step = size / 2
                )
        {
            if (col_indices[current + step] < row)
            {
                current += step;
                size -= step;
            }
            else if (col_indices[current + step] == row)
            {
                current += step;
                break;
            }
            else
            {
                size = step;
            }
        }

        if (n_col_elements && col_indices[current] == row)
        {
            element_count += 1;
        }
    }
    *p_n = element_count;
    //LEAVE_FUNCTION();
    return JMTX_RESULT_SUCCESS;
}

jmtx_result matrix_ccs_get_row(const jmtx_matrix_ccs* mtx, uint32_t row, uint32_t n, jmtx_scalar_t* p_elements, uint32_t* p_cols)
{

    if (!mtx)
    {
        //REPORT_ERROR_MESSAGE("Matrix pointer was null");
        //LEAVE_FUNCTION();
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTX_TYPE_CCS)
    {
        //REPORT_ERROR_MESSAGE("Matrix was not compressed column sparse");
        //LEAVE_FUNCTION();
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (row >= mtx->base.rows)
    {
        //REPORT_ERROR_MESSAGE("Matrix has %u rows but row %u was requested", mtx->base.rows, row);
        //LEAVE_FUNCTION();
        return JMTX_RESULT_INDEX_OUT_OF_BOUNDS;
    }
    if (!p_elements)
    {
        //REPORT_ERROR_MESSAGE("Element pointer was null");
        //LEAVE_FUNCTION();
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!p_cols)
    {
        //REPORT_ERROR_MESSAGE("Columns pointer was null");
        //LEAVE_FUNCTION();
        return JMTX_RESULT_NULL_PARAM;
    }

    uint32_t k = 0;
    for (uint32_t col = 0; col < mtx->base.cols && mtx->elements_before[col] != mtx->elements_before[mtx->base.cols] && k != n; ++col)
    {
        uint32_t current = 0;
        const uint32_t n_col_elements = mtx->elements_before[col + 1] - mtx->elements_before[col];
        const uint32_t* col_indices = mtx->indices + mtx->elements_before[col];
        for (
                uint32_t size = n_col_elements, step = n_col_elements / 2;
                step != 0;
                step = size / 2
                )
        {
            if (col_indices[current + step] < row)
            {
                current += step;
                size -= step;
            }
            else if (col_indices[current + step] == row)
            {
                current += step;
                break;
            }
            else
            {
                size = step;
            }
        }

        if (n_col_elements && col_indices[current] == row)
        {
            p_elements[k] = mtx->elements[mtx->elements_before[col] + current];
            p_cols[k] = col;
            k += 1;
        }
    }
    //LEAVE_FUNCTION();
    return JMTX_RESULT_SUCCESS;
}

jmtx_result matrix_ccs_transpose(const jmtx_matrix_ccs* restrict mtx, jmtx_matrix_ccs* restrict out)
{

    if (!mtx)
    {
        //REPORT_ERROR_MESSAGE("Input matrix pointer was null");
        //LEAVE_FUNCTION();
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTX_TYPE_CCS)
    {
        //REPORT_ERROR_MESSAGE("Input matrix was not compressed column sparse");
        //LEAVE_FUNCTION();
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (!out)
    {
        //REPORT_ERROR_MESSAGE("Output matrix pointer was null");
        //LEAVE_FUNCTION();
        return JMTX_RESULT_NULL_PARAM;
    }

    const uint32_t n_elements = mtx->n_elements;
//    const uint32_t new_rows = mtx->base.cols;
    const uint32_t new_cols = mtx->base.rows;

    uint32_t* const row_cum_counts = mtx->base.allocator_callbacks.alloc(mtx->base.allocator_callbacks.state, (new_cols + 1) * sizeof(*row_cum_counts));
    if (!row_cum_counts)
    {
//        CALLOC_FAILED((new_rows + 1) * sizeof*row_cum_counts);
        //LEAVE_FUNCTION();
        return JMTX_RESULT_BAD_ALLOC;
    }
    uint32_t* const new_indices = mtx->base.allocator_callbacks.alloc(mtx->base.allocator_callbacks.state, (n_elements + 1) * sizeof(*new_indices));
    if (!new_indices)
    {
        mtx->base.allocator_callbacks.free(mtx->base.allocator_callbacks.state, row_cum_counts);
//        CALLOC_FAILED((n_elements + 1) * sizeof*new_indices);
        //LEAVE_FUNCTION();
        return JMTX_RESULT_BAD_ALLOC;
    }
    jmtx_scalar_t* const new_elements = mtx->base.allocator_callbacks.alloc(mtx->base.allocator_callbacks.state, (n_elements + 1) * sizeof*new_elements);
    if (!new_elements)
    {
        mtx->base.allocator_callbacks.free(mtx->base.allocator_callbacks.state, row_cum_counts);
        mtx->base.allocator_callbacks.free(mtx->base.allocator_callbacks.state, new_indices);
//        CALLOC_FAILED((n_elements + 1) * sizeof*new_elements);
        //LEAVE_FUNCTION();
        return JMTX_RESULT_BAD_ALLOC;
    }

    *row_cum_counts = 1;

    for (uint32_t j = 0, n, p = 1; j < mtx->base.rows; ++j)
    {
        //  This MUST NOT fail, since the parameters are all within the bounds
        jmtx_result res = matrix_ccs_elements_in_row(mtx, j, &n);
        assert(res == JMTX_RESULT_SUCCESS);

        res = matrix_ccs_get_row(mtx, j, n, new_elements + p, new_indices + p);
        assert(res == JMTX_RESULT_SUCCESS);
        p += n;
        row_cum_counts[j + 1] = row_cum_counts[j] + n;
    }

    memcpy(out, mtx, sizeof*out);
    out->elements_before = row_cum_counts;
    out->elements = new_elements;
    out->indices = new_indices;
    out->capacity = n_elements;
    out->base = mtx->base;

    //LEAVE_FUNCTION();
    return 0;
}

jmtx_result matrix_ccs_copy(const jmtx_matrix_ccs* mtx, jmtx_matrix_ccs* out)
{
    if (!mtx)
    {
        //REPORT_ERROR_MESSAGE("Input matrix pointer was null");
        //LEAVE_FUNCTION();
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTX_TYPE_CCS)
    {
        //REPORT_ERROR_MESSAGE("Input matrix was not compressed column sparse");
        //LEAVE_FUNCTION();
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (!out)
    {
        //REPORT_ERROR_MESSAGE("Output matrix pointer was null");
        //LEAVE_FUNCTION();
        return JMTX_RESULT_NULL_PARAM;
    }
    
    jmtx_scalar_t* const elements = mtx->base.allocator_callbacks.alloc(mtx->base.allocator_callbacks.state, (1 + mtx->n_elements) * sizeof *elements);
    if (!elements)
    {
//        CALLOC_FAILED((1 + mtx->n_elements) * sizeof *elements);
        //LEAVE_FUNCTION();
        return JMTX_RESULT_BAD_ALLOC;
    }
    *elements = 0.0f;
    uint32_t* const indices = mtx->base.allocator_callbacks.alloc(mtx->base.allocator_callbacks.state, (1 + mtx->n_elements) * sizeof *indices);
    if (!indices)
    {
        mtx->base.allocator_callbacks.free(mtx->base.allocator_callbacks.state, elements);
//        CALLOC_FAILED((1 + mtx->n_elements) * sizeof *indices);
        //LEAVE_FUNCTION();
        return JMTX_RESULT_BAD_ALLOC;
    }
    *indices = 0;
    uint32_t* const cum_sum = mtx->base.allocator_callbacks.alloc(mtx->base.allocator_callbacks.state, (1 + mtx->base.cols) * sizeof *cum_sum);
    if (!cum_sum)
    {
        mtx->base.allocator_callbacks.free(mtx->base.allocator_callbacks.state, indices);
        mtx->base.allocator_callbacks.free(mtx->base.allocator_callbacks.state, elements);
//        CALLOC_FAILED((1 + mtx->base.rows) * sizeof *cum_sum);
        //LEAVE_FUNCTION();
        return JMTX_RESULT_BAD_ALLOC;
    }

    memcpy(elements + 1, mtx->elements + 1, sizeof* elements * mtx->n_elements);
    memcpy(indices + 1, mtx->indices + 1, sizeof* indices * mtx->n_elements);
    memcpy(cum_sum, mtx->elements_before, sizeof* cum_sum * (mtx->base.cols + 1));
    memcpy(out, mtx, sizeof *out);
    out->elements = elements;
    out->indices = indices;
    out->elements_before = cum_sum;
    //LEAVE_FUNCTION();
    return 0;
}
