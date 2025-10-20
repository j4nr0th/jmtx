#include "ccs_conversion.h"

/**
 * Creates a new CCS matrix with single precision from a CCS matrix with JMTX_SCALAR_T precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result JMTX_NAME_CONVERSION(matrix_ccs)(JMTX_NAME_OUT(matrix_ccs) * *p_mtx, const JMTX_NAME_IN(matrix_ccs) * in,
                                             const jmtx_allocator_callbacks *allocator_callbacks)
{
    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }

    JMTX_INDEX_T *offsets = NULL;
    JMTX_INDEX_T *indices = NULL;

    JMTX_NAME_OUT(matrix_ccs) *mtx = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*mtx));
    if (!mtx)
    {
        return JMTX_RESULT_BAD_ALLOC;
    }

    JMTX_OUTPUT_TYPE *values =
        allocator_callbacks->alloc(allocator_callbacks->state, (in->n_entries) * sizeof(*values));
    if (!values)
    {
        allocator_callbacks->free(allocator_callbacks->state, mtx);
        return JMTX_RESULT_BAD_ALLOC;
    }

    indices = allocator_callbacks->alloc(allocator_callbacks->state, (in->n_entries) * sizeof(*indices));
    if (!indices)
    {
        allocator_callbacks->free(allocator_callbacks->state, indices);
        allocator_callbacks->free(allocator_callbacks->state, mtx);
        return JMTX_RESULT_BAD_ALLOC;
    }

    offsets = allocator_callbacks->alloc(allocator_callbacks->state, (in->base.cols) * sizeof(*offsets));
    if (!offsets)
    {
        allocator_callbacks->free(allocator_callbacks->state, offsets);
        allocator_callbacks->free(allocator_callbacks->state, indices);
        allocator_callbacks->free(allocator_callbacks->state, mtx);
        return JMTX_RESULT_BAD_ALLOC;
    }
    memcpy(offsets, in->end_of_column_offsets, (in->base.cols) * sizeof(*offsets));
    memcpy(indices, in->indices, (in->n_entries) * sizeof(*indices));

    for (JMTX_FAST_INT_T i = 0; i < in->n_entries; ++i)
    {
        values[i] = (JMTX_OUTPUT_TYPE)in->values[i];
    }

    mtx->base.cols = in->base.cols;
    mtx->base.type = JMTX_OUTPUT_ENUM(CCS);
    mtx->base.rows = in->base.rows;
    mtx->base.allocator_callbacks = *allocator_callbacks;
    mtx->indices = indices;
    mtx->values = values;
    mtx->capacity = in->n_entries;
    mtx->n_entries = in->n_entries;
    mtx->end_of_column_offsets = offsets;
    *p_mtx = mtx;

    return JMTX_RESULT_SUCCESS;
}

/**
 * Creates a new CCS matrix with single precision from a CCS matrix with JMTX_SCALAR_T precision. Requires minimum
 * memory allocation by reusing the memory of the initial matrix.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert (will be invalid if function succeeds)
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
JMTX_NAME_OUT(matrix_ccs) * JMTX_NAME_CONVERSION(matrix_ccs_inplace)(JMTX_NAME_IN(matrix_ccs) * in)
{
    const JMTX_FAST_INT_T count = in->n_entries;
    if (sizeof(JMTX_INPUT_TYPE) < sizeof(JMTX_OUTPUT_TYPE))
    {
        JMTX_OUTPUT_TYPE *const values = (JMTX_OUTPUT_TYPE *)in->values;
        // Front-to-back copy
        for (JMTX_FAST_INT_T i = 0; i < count; ++i)
        {
            values[i] = (JMTX_OUTPUT_TYPE)in->values[i];
        }

        // Reallocation
        JMTX_OUTPUT_TYPE *const new_ptr =
            in->base.allocator_callbacks.realloc(in->base.allocator_callbacks.state, values, sizeof(*new_ptr) * count);
        if (new_ptr)
        {
            // Failure is not bad
            in->values = (JMTX_INPUT_TYPE *)new_ptr;
        }
    }
    else
    {
        // First, try to reallocate
        JMTX_OUTPUT_TYPE *const new_ptr = in->base.allocator_callbacks.realloc(in->base.allocator_callbacks.state,
                                                                               in->values, sizeof(*new_ptr) * count);
        if (!new_ptr)
        {
            // Failure is bad here
            return NULL;
        }
        in->values = (JMTX_INPUT_TYPE *)new_ptr;
        JMTX_OUTPUT_TYPE *const values = (JMTX_OUTPUT_TYPE *)in->values;

        // Copy back-to-front
        for (JMTX_FAST_INT_T i = 0; i < count; ++i)
        {
            values[count - 1 - i] = (JMTX_OUTPUT_TYPE)in->values[count - 1 - i];
        }
    }

    in->capacity = count;
    in->base.type = JMTX_OUTPUT_ENUM(CCS);

    return (JMTX_NAME_OUT(matrix_ccs) *)in;
}
