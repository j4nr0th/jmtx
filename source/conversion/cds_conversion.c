#include "cds_conversion.h"

/**
 * Creates a new CDS matrix with single precision from a CDS matrix with JMTX_INPUT_TYPE precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result JMTX_NAME_CONVERSION(matrix_cds)(JMTX_NAME_OUT(matrix_cds) * *p_mtx, const JMTX_NAME_IN(matrix_cds) * in,
                                             const jmtx_allocator_callbacks *allocator_callbacks)
{
    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }
    JMTX_NAME_OUT(matrix_cds) * mtx;
    jmtx_result res =
        JMTX_NAME_OUT(matrix_cds_new)(&mtx, in->base.rows, in->base.cols, 0, (int32_t[]){0}, allocator_callbacks);
    if (res != JMTX_RESULT_SUCCESS)
    {
        return res;
    }

    for (JMTX_FAST_INT_T i = 0; i < in->sub_diagonals.count; ++i)
    {
        const int32_t dia_idx = -(int32_t)in->sub_diagonals.indices[i];
        const JMTX_INPUT_TYPE *const dia_ptr = in->sub_diagonals.diagonals[i];
        JMTX_INDEX_T len;
        JMTX_OUTPUT_TYPE *const ptr = JMTX_NAME_OUT(matrix_cds_allocate_diagonal)(mtx, dia_idx, &len);
        if (ptr == NULL)
        {
            JMTX_NAME_OUT(matrix_cds_destroy)(mtx);
            return JMTX_RESULT_BAD_ALLOC;
        }
        for (JMTX_FAST_INT_T j = 0; j < len; ++j)
        {
            ptr[j] = (JMTX_OUTPUT_TYPE)dia_ptr[j];
        }
    }
    if (in->main_diagonal)
    {
        JMTX_INDEX_T len;
        JMTX_OUTPUT_TYPE *const ptr = JMTX_NAME_OUT(matrix_cds_allocate_diagonal)(mtx, 0, &len);
        if (ptr == NULL)
        {
            JMTX_NAME_OUT(matrix_cds_destroy)(mtx);
            return JMTX_RESULT_BAD_ALLOC;
        }
        for (JMTX_FAST_INT_T j = 0; j < len; ++j)
        {
            ptr[j] = (JMTX_OUTPUT_TYPE)in->main_diagonal[j];
        }
    }
    for (JMTX_FAST_INT_T i = 0; i < in->super_diagonals.count; ++i)
    {
        const int32_t dia_idx = (int32_t)in->super_diagonals.indices[i];
        const JMTX_INPUT_TYPE *const dia_ptr = in->super_diagonals.diagonals[i];
        JMTX_INDEX_T len;
        JMTX_OUTPUT_TYPE *const ptr = JMTX_NAME_OUT(matrix_cds_allocate_diagonal)(mtx, dia_idx, &len);
        if (ptr == NULL)
        {
            JMTX_NAME_OUT(matrix_cds_destroy)(mtx);
            return JMTX_RESULT_BAD_ALLOC;
        }
        for (JMTX_FAST_INT_T j = 0; j < len; ++j)
        {
            ptr[j] = (JMTX_OUTPUT_TYPE)dia_ptr[j];
        }
    }

    mtx->base.cols = in->base.cols;
    mtx->base.type = JMTX_OUTPUT_ENUM(CDS);
    mtx->base.rows = in->base.rows;
    mtx->base.allocator_callbacks = *allocator_callbacks;

    *p_mtx = mtx;

    return JMTX_RESULT_SUCCESS;
}

/**
 * Creates a new CDS matrix with single precision from a CDS matrix with JMTX_INPUT_TYPE precision. Requires no memory
 * allocation by reusing the memory of the initial matrix. Can not fail if the input matrix is valid.
 * @param in matrix which to convert (will be invalid if function succeeds)
 * @return converted matrix
 */
JMTX_NAME_OUT(matrix_cds) * JMTX_NAME_CONVERSION(matrix_cds_inplace)(JMTX_NAME_IN(matrix_cds) * in)
{
    if (sizeof(JMTX_INPUT_TYPE) < sizeof(JMTX_OUTPUT_TYPE))
    {

        //  Convert diagonals
        if (in->main_diagonal)
        {
            JMTX_OUTPUT_TYPE *const values = (JMTX_OUTPUT_TYPE *)in->main_diagonal;
            const JMTX_INDEX_T len = JMTX_NAME_IN(matrix_cds_entries_in_dia)(in, 0);
            for (JMTX_INDEX_T i = 0; i < len; ++i)
            {
                values[i] = (JMTX_OUTPUT_TYPE)in->main_diagonal[i];
            }
        }

        for (JMTX_FAST_INT_T i = 0; i < in->sub_diagonals.count; ++i)
        {
            JMTX_INPUT_TYPE *ptr = in->sub_diagonals.diagonals[i];
            JMTX_OUTPUT_TYPE *const values = (JMTX_OUTPUT_TYPE *)ptr;
            const JMTX_INDEX_T len = JMTX_NAME_IN(matrix_cds_entries_in_dia)(in, -((int32_t)i));
            for (JMTX_INDEX_T j = 0; j < len; ++j)
            {
                values[j] = (JMTX_OUTPUT_TYPE)ptr[j];
            }
        }

        for (JMTX_FAST_INT_T i = 0; i < in->super_diagonals.count; ++i)
        {
            JMTX_INPUT_TYPE *ptr = in->super_diagonals.diagonals[i];
            JMTX_OUTPUT_TYPE *const values = (JMTX_OUTPUT_TYPE *)ptr;
            const JMTX_INDEX_T len = JMTX_NAME_IN(matrix_cds_entries_in_dia)(in, +((int32_t)i));
            for (JMTX_INDEX_T j = 0; j < len; ++j)
            {
                values[j] = (JMTX_OUTPUT_TYPE)ptr[j];
            }
        }

        //  Shrink the diagonals
        if (in->main_diagonal)
        {
            const JMTX_INDEX_T len = JMTX_NAME_IN(matrix_cds_entries_in_dia)(in, 0);
            JMTX_OUTPUT_TYPE *const new_ptr = in->base.allocator_callbacks.realloc(
                in->base.allocator_callbacks.state, in->main_diagonal, sizeof(*new_ptr) * len);
            if (new_ptr)
            {
                in->main_diagonal = (JMTX_INPUT_TYPE *)new_ptr;
            }
        }

        for (JMTX_FAST_INT_T i = 0; i < in->sub_diagonals.count; ++i)
        {
            const JMTX_INDEX_T len = JMTX_NAME_IN(matrix_cds_entries_in_dia)(in, -((int32_t)i));
            JMTX_OUTPUT_TYPE *const new_ptr = in->base.allocator_callbacks.realloc(
                in->base.allocator_callbacks.state, in->sub_diagonals.diagonals[i], sizeof(*new_ptr) * len);
            if (new_ptr)
            {
                in->sub_diagonals.diagonals[i] = (JMTX_INPUT_TYPE *)new_ptr;
            }
        }

        for (JMTX_FAST_INT_T i = 0; i < in->super_diagonals.count; ++i)
        {
            const JMTX_INDEX_T len = JMTX_NAME_IN(matrix_cds_entries_in_dia)(in, +((int32_t)i));
            JMTX_OUTPUT_TYPE *const new_ptr = in->base.allocator_callbacks.realloc(
                in->base.allocator_callbacks.state, in->super_diagonals.diagonals[i], sizeof(*new_ptr) * len);
            if (new_ptr)
            {
                in->super_diagonals.diagonals[i] = (JMTX_INPUT_TYPE *)new_ptr;
            }
        }
    }
    else
    {
        // Expand the diagonals
        // Main diagonal
        if (in->main_diagonal)
        {
            const JMTX_INDEX_T len = JMTX_NAME_IN(matrix_cds_entries_in_dia)(in, 0);
            JMTX_OUTPUT_TYPE *const new_ptr = in->base.allocator_callbacks.realloc(
                in->base.allocator_callbacks.state, in->main_diagonal, sizeof(*new_ptr) * len);
            if (!new_ptr)
            {
                return NULL;
            }
            in->main_diagonal = (JMTX_INPUT_TYPE *)new_ptr;
        }
        // sub-diagonals
        for (JMTX_FAST_INT_T i = 0; i < in->sub_diagonals.count; ++i)
        {
            const JMTX_INDEX_T len = JMTX_NAME_IN(matrix_cds_entries_in_dia)(in, -((int32_t)i));
            JMTX_OUTPUT_TYPE *const new_ptr = in->base.allocator_callbacks.realloc(
                in->base.allocator_callbacks.state, in->sub_diagonals.diagonals[i], sizeof(*new_ptr) * len);
            if (!new_ptr)
            {
                return NULL;
            }
            in->sub_diagonals.diagonals[i] = (JMTX_INPUT_TYPE *)new_ptr;
        }
        // super-diagonals
        for (JMTX_FAST_INT_T i = 0; i < in->super_diagonals.count; ++i)
        {
            const JMTX_INDEX_T len = JMTX_NAME_IN(matrix_cds_entries_in_dia)(in, +((int32_t)i));
            JMTX_OUTPUT_TYPE *const new_ptr = in->base.allocator_callbacks.realloc(
                in->base.allocator_callbacks.state, in->super_diagonals.diagonals[i], sizeof(*new_ptr) * len);
            if (!new_ptr)
            {
                return NULL;
            }
            in->super_diagonals.diagonals[i] = (JMTX_INPUT_TYPE *)new_ptr;
        }

        // Convert the diagonals going back-to-front
        // Main diagonal
        if (in->main_diagonal)
        {
            JMTX_OUTPUT_TYPE *const values = (JMTX_OUTPUT_TYPE *)in->main_diagonal;
            const JMTX_INDEX_T len = JMTX_NAME_IN(matrix_cds_entries_in_dia)(in, 0);
            for (JMTX_INDEX_T i = 0; i < len; ++i)
            {
                values[len - 1 - i] = (JMTX_OUTPUT_TYPE)in->main_diagonal[len - 1 - i];
            }
        }
        // Sub-diagonals
        for (JMTX_FAST_INT_T i = 0; i < in->sub_diagonals.count; ++i)
        {
            JMTX_INPUT_TYPE *ptr = in->sub_diagonals.diagonals[i];
            JMTX_OUTPUT_TYPE *const values = (JMTX_OUTPUT_TYPE *)ptr;
            const JMTX_INDEX_T len = JMTX_NAME_IN(matrix_cds_entries_in_dia)(in, -((int32_t)i));
            for (JMTX_INDEX_T j = 0; j < len; ++j)
            {
                values[len - 1 - j] = (JMTX_OUTPUT_TYPE)ptr[len - 1 - j];
            }
        }
        // Super-diagonals
        for (JMTX_FAST_INT_T i = 0; i < in->super_diagonals.count; ++i)
        {
            JMTX_INPUT_TYPE *ptr = in->super_diagonals.diagonals[i];
            JMTX_OUTPUT_TYPE *const values = (JMTX_OUTPUT_TYPE *)ptr;
            const JMTX_INDEX_T len = JMTX_NAME_IN(matrix_cds_entries_in_dia)(in, +((int32_t)i));
            for (JMTX_INDEX_T j = 0; j < len; ++j)
            {
                values[len - 1 - j] = (JMTX_OUTPUT_TYPE)ptr[len - 1 - j];
            }
        }
    }

    in->base.type = JMTX_OUTPUT_ENUM(CDS);

    return (JMTX_NAME_OUT(matrix_cds) *)in;
}
