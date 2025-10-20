#include "brm_conversion.h"

static JMTX_FAST_INT_T JMTX_NAME_CONVERSION(brm_row_offset)(const JMTX_NAME_IN(matrix_brm)
                                                                mtx[JMTX_ARRAY_ATTRIB(const static 1)],
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

/**
 * Creates a new BRM matrix with single precision from a BRM matrix with JMTX_SCALAR_T precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result JMTX_NAME_CONVERSION(matrix_brm)(JMTX_NAME_OUT(matrix_brm) * *p_mtx, const JMTX_NAME_IN(matrix_brm) * in,
                                             const jmtx_allocator_callbacks *allocator_callbacks)
{
    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }

    JMTX_NAME_OUT(matrix_brm) *mtx = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*mtx));
    if (!mtx)
    {
        return JMTX_RESULT_BAD_ALLOC;
    }

    const JMTX_FAST_INT_T count = JMTX_NAME_CONVERSION(brm_row_offset)(in, in->base.rows);
    JMTX_OUTPUT_TYPE *values = allocator_callbacks->alloc(allocator_callbacks->state, count * sizeof(*values));
    if (!values)
    {
        allocator_callbacks->free(allocator_callbacks->state, mtx);
        return JMTX_RESULT_BAD_ALLOC;
    }

    for (JMTX_FAST_INT_T i = 0; i < count; ++i)
    {
        values[i] = (JMTX_OUTPUT_TYPE)in->values[i];
    }

    mtx->base.cols = in->base.cols;
    mtx->base.type = JMTX_OUTPUT_ENUM(BRM);
    mtx->base.rows = in->base.rows;
    mtx->base.allocator_callbacks = *allocator_callbacks;
    mtx->values = values;
    mtx->lower_bandwidth = in->lower_bandwidth;
    mtx->upper_bandwidth = in->upper_bandwidth;

    *p_mtx = mtx;

    return JMTX_RESULT_SUCCESS;
}
/**
 * Creates a new BRM matrix with single precision from a BRM matrix with JMTX_SCALAR_T precision. Requires no memory
 * allocation by reusing the memory of the initial matrix. Can not fail if the input matrix is valid.
 * @param in matrix which to convert (will be invalid if function succeeds)
 * @return converted matrix
 */
JMTX_NAME_OUT(matrix_brm) * JMTX_NAME_CONVERSION(matrix_brm_inplace)(JMTX_NAME_IN(matrix_brm) * in)
{
    const JMTX_FAST_INT_T count = JMTX_NAME_CONVERSION(brm_row_offset)(in, in->base.rows);

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

    in->base.type = JMTX_OUTPUT_ENUM(BRM);
    return (JMTX_NAME_OUT(matrix_brm) *)in;
}
