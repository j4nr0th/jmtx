//
// Created by jan on 1.12.2023.
//

#include <assert.h>
#include <complex.h>

#include "../double/matrices/sparse_diagonal_compressed.h"
#include "../float/matrices/sparse_diagonal_compressed.h"
#include "../matrix_base.h"
#ifndef _MSC_BUILD
#    include "../cdouble/matrices/sparse_diagonal_compressed.h"
#    include "../cfloat/matrices/sparse_diagonal_compressed.h"
#endif
#include "cds_conversion.h"

/***********************************************************************************************************************
 *                                                                                                                     *
 *                                          FLOAT <-> DOUBLE                                                           *
 *                                                                                                                     *
 **********************************************************************************************************************/

/**
 * Creates a new CDS matrix with single precision from a CDS matrix with JMTX_SCALAR_T precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtx_matrix_cds_from_double(jmtx_matrix_cds **p_mtx, const JMTX_NAME_TYPED(matrix_cds) * in,
                                        const jmtx_allocator_callbacks *allocator_callbacks)
{
    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }
    jmtx_matrix_cds *mtx;
    jmtx_result res = jmtx_matrix_cds_new(&mtx, in->base.rows, in->base.cols, 0, (int32_t[]){0}, allocator_callbacks);
    if (res != JMTX_RESULT_SUCCESS)
    {
        return res;
    }

    for (JMTX_FAST_INT_T i = 0; i < in->sub_diagonals.count; ++i)
    {
        const int32_t dia_idx = -(int32_t)in->sub_diagonals.indices[i];
        const JMTX_SCALAR_T *const dia_ptr = in->sub_diagonals.diagonals[i];
        JMTX_INDEX_T len;
        float *const ptr = jmtx_matrix_cds_allocate_diagonal(mtx, dia_idx, &len);
        for (JMTX_FAST_INT_T j = 0; j < len; ++j)
        {
            ptr[j] = (float)dia_ptr[j];
        }
    }
    if (in->main_diagonal)
    {
        JMTX_INDEX_T len;
        float *const ptr = jmtx_matrix_cds_allocate_diagonal(mtx, 0, &len);
        for (JMTX_FAST_INT_T j = 0; j < len; ++j)
        {
            ptr[j] = (float)in->main_diagonal[j];
        }
    }
    for (JMTX_FAST_INT_T i = 0; i < in->super_diagonals.count; ++i)
    {
        const int32_t dia_idx = (int32_t)in->super_diagonals.indices[i];
        const JMTX_SCALAR_T *const dia_ptr = in->super_diagonals.diagonals[i];
        JMTX_INDEX_T len;
        float *const ptr = jmtx_matrix_cds_allocate_diagonal(mtx, dia_idx, &len);
        for (JMTX_FAST_INT_T j = 0; j < len; ++j)
        {
            ptr[j] = (float)dia_ptr[j];
        }
    }

    mtx->base.cols = in->base.cols;
    mtx->base.type = JMTX_TYPE_CDS;
    mtx->base.rows = in->base.rows;
    mtx->base.allocator_callbacks = *allocator_callbacks;

    *p_mtx = mtx;

    return JMTX_RESULT_SUCCESS;
}

/**
 * Creates a new CDS matrix with single precision from a CDS matrix with JMTX_SCALAR_T precision. Requires no memory
 * allocation by reusing the memory of the initial matrix. Can not fail if the input matrix is valid.
 * @param in matrix which to convert (will be invalid if function succeeds)
 * @return converted matrix
 */
jmtx_matrix_cds *jmtx_matrix_cds_from_double_inplace(JMTX_NAME_TYPED(matrix_cds) * in)
{
    //  Convert diagonals
    if (in->main_diagonal)
    {
        float *const values = (float *)in->main_diagonal;
        const JMTX_INDEX_T len = JMTX_NAME_TYPED(matrix_cds_entries_in_dia(in, 0);
        for (JMTX_INDEX_T i = 0; i < len; ++i)
        {
            values[i] = (float)in->main_diagonal[i];
        }
    }

    for (JMTX_FAST_INT_T i = 0; i < in->sub_diagonals.count; ++i)
    {
        JMTX_SCALAR_T *ptr = in->sub_diagonals.diagonals[i];
        float *const values = (float *)ptr;
        const JMTX_INDEX_T len = JMTX_NAME_TYPED(matrix_cds_entries_in_dia(in, -((int32_t)i));
        for (JMTX_INDEX_T j = 0; j < len; ++j)
        {
            values[j] = (float)ptr[j];
        }
    }

    for (JMTX_FAST_INT_T i = 0; i < in->super_diagonals.count; ++i)
    {
        JMTX_SCALAR_T *ptr = in->super_diagonals.diagonals[i];
        float *const values = (float *)ptr;
        const JMTX_INDEX_T len = JMTX_NAME_TYPED(matrix_cds_entries_in_dia(in, +((int32_t)i));
        for (JMTX_INDEX_T j = 0; j < len; ++j)
        {
            values[j] = (float)ptr[j];
        }
    }

    //  Shrink the diagonals
    if (in->main_diagonal)
    {
        const JMTX_INDEX_T len = JMTX_NAME_TYPED(matrix_cds_entries_in_dia(in, 0);
        float *const new_ptr = in->base.allocator_callbacks.realloc(in->base.allocator_callbacks.state,
                                                                    in->main_diagonal, sizeof(*new_ptr) * len);
        if (new_ptr)
        {
            in->main_diagonal = (JMTX_SCALAR_T *)new_ptr;
        }
    }

    for (JMTX_FAST_INT_T i = 0; i < in->sub_diagonals.count; ++i)
    {
        const JMTX_INDEX_T len = JMTX_NAME_TYPED(matrix_cds_entries_in_dia(in, -((int32_t)i));
        float *const new_ptr = in->base.allocator_callbacks.realloc(
            in->base.allocator_callbacks.state, in->sub_diagonals.diagonals[i], sizeof(*new_ptr) * len);
        if (new_ptr)
        {
            in->sub_diagonals.diagonals[i] = (JMTX_SCALAR_T *)new_ptr;
        }
    }

    for (JMTX_FAST_INT_T i = 0; i < in->super_diagonals.count; ++i)
    {
        const JMTX_INDEX_T len = JMTX_NAME_TYPED(matrix_cds_entries_in_dia(in, +((int32_t)i));
        float *const new_ptr = in->base.allocator_callbacks.realloc(
            in->base.allocator_callbacks.state, in->super_diagonals.diagonals[i], sizeof(*new_ptr) * len);
        if (new_ptr)
        {
            in->super_diagonals.diagonals[i] = (JMTX_SCALAR_T *)new_ptr;
        }
    }

    in->base.type = JMTX_TYPE_CDS;

    return (jmtx_matrix_cds *)in;
}

/**
 * Creates a new CDS matrix with single precision from a CDS matrix with JMTX_SCALAR_T precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxs_matrix_cds_from_double(jmtx_matrix_cds **p_mtx, const JMTX_NAME_TYPED(matrix_cds) * in,
                                         const jmtx_allocator_callbacks *allocator_callbacks)
{
    if (!p_mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!in)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (allocator_callbacks &&
        (!allocator_callbacks->free || !allocator_callbacks->alloc || !allocator_callbacks->realloc))
    {
        return JMTX_RESULT_BAD_PARAM;
    }
    if (in->base.type != JMTXD_TYPE_CDS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }

    return jmtx_matrix_cds_from_double(p_mtx, in, allocator_callbacks);
}

/**
 * Creates a new CDS matrix with JMTX_SCALAR_T precision from a CDS matrix with single precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result JMTX_NAME_TYPED(matrix_cds_from_float(JMTX_NAME_TYPED(matrix_cds) **p_mtx, const jmtx_matrix_cds *in,
                                        const jmtx_allocator_callbacks *allocator_callbacks)
{
    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }
    JMTX_NAME_TYPED(matrix_cds) * mtx;
    jmtx_result res = JMTX_NAME_TYPED(matrix_cds_new(&mtx, in->base.rows, in->base.cols, 0, (int32_t[]){0}, allocator_callbacks);
    if (res != JMTX_RESULT_SUCCESS)
    {
        return res;
    }

    for (JMTX_FAST_INT_T i = 0; i < in->sub_diagonals.count; ++i)
    {
        const int32_t dia_idx = -(int32_t)in->sub_diagonals.indices[i];
        const float *const dia_ptr = in->sub_diagonals.diagonals[i];
        JMTX_INDEX_T len;
        JMTX_SCALAR_T *const ptr = JMTX_NAME_TYPED(matrix_cds_allocate_diagonal(mtx, dia_idx, &len);
        for (JMTX_FAST_INT_T j = 0; j < len; ++j)
        {
            ptr[j] = (JMTX_SCALAR_T)dia_ptr[j];
        }
    }
    if (in->main_diagonal)
    {
        JMTX_INDEX_T len;
        JMTX_SCALAR_T *const ptr = JMTX_NAME_TYPED(matrix_cds_allocate_diagonal(mtx, 0, &len);
        for (JMTX_FAST_INT_T j = 0; j < len; ++j)
        {
            ptr[j] = (JMTX_SCALAR_T)in->main_diagonal[j];
        }
    }
    for (JMTX_FAST_INT_T i = 0; i < in->super_diagonals.count; ++i)
    {
        const int32_t dia_idx = (int32_t)in->super_diagonals.indices[i];
        const float *const dia_ptr = in->super_diagonals.diagonals[i];
        JMTX_INDEX_T len;
        JMTX_SCALAR_T *const ptr = JMTX_NAME_TYPED(matrix_cds_allocate_diagonal(mtx, dia_idx, &len);
        for (JMTX_FAST_INT_T j = 0; j < len; ++j)
        {
            ptr[j] = (JMTX_SCALAR_T)dia_ptr[j];
        }
    }

    mtx->base.cols = in->base.cols;
    mtx->base.type = JMTXD_TYPE_CDS;
    mtx->base.rows = in->base.rows;
    mtx->base.allocator_callbacks = *allocator_callbacks;

    *p_mtx = mtx;

    return JMTX_RESULT_SUCCESS;
}

/**
 * Creates a new CDS matrix with JMTX_SCALAR_T precision from a CDS matrix with single precision. Only one memory reallocation
 * may be needed. Can not fail if the input matrix is valid.
 * @param in matrix which to convert
 * @return converted matrix, or NULL in case of allocation failure
 */
JMTX_NAME_TYPED(matrix_cds) *JMTX_NAME_TYPED(matrix_cds_from_float_inplace(jmtx_matrix_cds *in)
{
    //  Expand the diagonals
    if (in->main_diagonal)
    {
        const JMTX_INDEX_T len = jmtx_matrix_cds_entries_in_dia(in, 0);
        JMTX_SCALAR_T *const new_ptr = in->base.allocator_callbacks.realloc(in->base.allocator_callbacks.state,
                                                                            in->main_diagonal, sizeof(*new_ptr) * len);
        if (!new_ptr)
        {
            return NULL;
        }
        in->main_diagonal = (float *)new_ptr;
    }

    for (JMTX_FAST_INT_T i = 0; i < in->sub_diagonals.count; ++i)
    {
        const JMTX_INDEX_T len = jmtx_matrix_cds_entries_in_dia(in, -((int32_t)i));
        JMTX_SCALAR_T *const new_ptr = in->base.allocator_callbacks.realloc(
            in->base.allocator_callbacks.state, in->sub_diagonals.diagonals[i], sizeof(*new_ptr) * len);
        if (!new_ptr)
        {
            return NULL;
        }
        in->sub_diagonals.diagonals[i] = (float *)new_ptr;
    }

    for (JMTX_FAST_INT_T i = 0; i < in->super_diagonals.count; ++i)
    {
        const JMTX_INDEX_T len = jmtx_matrix_cds_entries_in_dia(in, +((int32_t)i));
        JMTX_SCALAR_T *const new_ptr = in->base.allocator_callbacks.realloc(
            in->base.allocator_callbacks.state, in->super_diagonals.diagonals[i], sizeof(*new_ptr) * len);
        if (!new_ptr)
        {
            return NULL;
        }
        in->super_diagonals.diagonals[i] = (float *)new_ptr;
    }

    //  Convert diagonals
    if (in->main_diagonal)
    {
        JMTX_SCALAR_T *const values = (JMTX_SCALAR_T *)in->main_diagonal;
        const JMTX_INDEX_T len = jmtx_matrix_cds_entries_in_dia(in, 0);
        for (JMTX_INDEX_T i = 0; i < len; ++i)
        {
            values[len - 1 - i] = (JMTX_SCALAR_T)in->main_diagonal[len - 1 - i];
        }
    }

    for (JMTX_FAST_INT_T i = 0; i < in->sub_diagonals.count; ++i)
    {
        float *ptr = in->sub_diagonals.diagonals[i];
        JMTX_SCALAR_T *const values = (JMTX_SCALAR_T *)ptr;
        const JMTX_INDEX_T len = jmtx_matrix_cds_entries_in_dia(in, -((int32_t)i));
        for (JMTX_INDEX_T j = 0; j < len; ++j)
        {
            values[len - 1 - j] = (JMTX_SCALAR_T)ptr[len - 1 - j];
        }
    }

    for (JMTX_FAST_INT_T i = 0; i < in->super_diagonals.count; ++i)
    {
        float *ptr = in->super_diagonals.diagonals[i];
        JMTX_SCALAR_T *const values = (JMTX_SCALAR_T *)ptr;
        const JMTX_INDEX_T len = jmtx_matrix_cds_entries_in_dia(in, +((int32_t)i));
        for (JMTX_INDEX_T j = 0; j < len; ++j)
        {
            values[len - 1 - j] = (JMTX_SCALAR_T)ptr[len - 1 - j];
        }
    }

    in->base.type = JMTXD_TYPE_CDS;

    return (JMTX_NAME_TYPED(matrix_cds) *)in;
}

/**
 * Creates a new CDS matrix with JMTX_SCALAR_T precision from a CDS matrix with single precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result .,.,.,.,..matrix_cds_from_float(JMTX_NAME_TYPED(matrix_cds) **p_mtx, const jmtx_matrix_cds *in,
                                         const jmtx_allocator_callbacks *allocator_callbacks)
{
    if (!p_mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!in)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (allocator_callbacks &&
        (!allocator_callbacks->free || !allocator_callbacks->alloc || !allocator_callbacks->realloc))
    {
        return JMTX_RESULT_BAD_PARAM;
    }
    if (in->base.type != JMTX_TYPE_CDS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }

    return JMTX_NAME_TYPED(matrix_cds_from_float(p_mtx, in, allocator_callbacks);
}

#ifndef _MSC_BUILD
/***********************************************************************************************************************
 *                                                                                                                     *
 *                                          FLOAT <-> COMPLEX FLOAT                                                    *
 *                                                                                                                     *
 **********************************************************************************************************************/

/**
 * Creates a new CDS matrix with single precision from the real part of a complex CDS matrix with single precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtx_matrix_cds_from_cfloat_real(jmtx_matrix_cds **p_mtx, const jmtxc_matrix_cds *in,
                                             const jmtx_allocator_callbacks *allocator_callbacks)
{
    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }
    jmtx_matrix_cds *mtx;
    jmtx_result res = jmtx_matrix_cds_new(&mtx, in->base.rows, in->base.cols, 0, (int32_t[]){0}, allocator_callbacks);
    if (res != JMTX_RESULT_SUCCESS)
    {
        return res;
    }

    for (JMTX_FAST_INT_T i = 0; i < in->sub_diagonals.count; ++i)
    {
        const int32_t dia_idx = -(int32_t)in->sub_diagonals.indices[i];
        const _Complex float *const dia_ptr = in->sub_diagonals.diagonals[i];
        JMTX_INDEX_T len;
        float *const ptr = jmtx_matrix_cds_allocate_diagonal(mtx, dia_idx, &len);
        for (JMTX_FAST_INT_T j = 0; j < len; ++j)
        {
            ptr[j] = crealf(dia_ptr[j]);
        }
    }
    if (in->main_diagonal)
    {
        JMTX_INDEX_T len;
        float *const ptr = jmtx_matrix_cds_allocate_diagonal(mtx, 0, &len);
        for (JMTX_FAST_INT_T j = 0; j < len; ++j)
        {
            ptr[j] = crealf(in->main_diagonal[j]);
        }
    }
    for (JMTX_FAST_INT_T i = 0; i < in->super_diagonals.count; ++i)
    {
        const int32_t dia_idx = (int32_t)in->super_diagonals.indices[i];
        const _Complex float *const dia_ptr = in->super_diagonals.diagonals[i];
        JMTX_INDEX_T len;
        float *const ptr = jmtx_matrix_cds_allocate_diagonal(mtx, dia_idx, &len);
        for (JMTX_FAST_INT_T j = 0; j < len; ++j)
        {
            ptr[j] = crealf(dia_ptr[j]);
        }
    }

    mtx->base.cols = in->base.cols;
    mtx->base.type = JMTX_TYPE_CDS;
    mtx->base.rows = in->base.rows;
    mtx->base.allocator_callbacks = *allocator_callbacks;

    *p_mtx = mtx;

    return JMTX_RESULT_SUCCESS;
}

/**
 * Creates a new CDS matrix with single precision from the imaginary part of a complex CDS matrix with single precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtx_matrix_cds_from_cfloat_imag(jmtx_matrix_cds **p_mtx, const jmtxc_matrix_cds *in,
                                             const jmtx_allocator_callbacks *allocator_callbacks)
{
    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }
    jmtx_matrix_cds *mtx;
    jmtx_result res = jmtx_matrix_cds_new(&mtx, in->base.rows, in->base.cols, 0, (int32_t[]){0}, allocator_callbacks);
    if (res != JMTX_RESULT_SUCCESS)
    {
        return res;
    }

    for (JMTX_FAST_INT_T i = 0; i < in->sub_diagonals.count; ++i)
    {
        const int32_t dia_idx = -(int32_t)in->sub_diagonals.indices[i];
        const _Complex float *const dia_ptr = in->sub_diagonals.diagonals[i];
        JMTX_INDEX_T len;
        float *const ptr = jmtx_matrix_cds_allocate_diagonal(mtx, dia_idx, &len);
        for (JMTX_FAST_INT_T j = 0; j < len; ++j)
        {
            ptr[j] = cimagf(dia_ptr[j]);
        }
    }
    if (in->main_diagonal)
    {
        JMTX_INDEX_T len;
        float *const ptr = jmtx_matrix_cds_allocate_diagonal(mtx, 0, &len);
        for (JMTX_FAST_INT_T j = 0; j < len; ++j)
        {
            ptr[j] = cimagf(in->main_diagonal[j]);
        }
    }
    for (JMTX_FAST_INT_T i = 0; i < in->super_diagonals.count; ++i)
    {
        const int32_t dia_idx = (int32_t)in->super_diagonals.indices[i];
        const _Complex float *const dia_ptr = in->super_diagonals.diagonals[i];
        JMTX_INDEX_T len;
        float *const ptr = jmtx_matrix_cds_allocate_diagonal(mtx, dia_idx, &len);
        for (JMTX_FAST_INT_T j = 0; j < len; ++j)
        {
            ptr[j] = cimagf(dia_ptr[j]);
        }
    }

    mtx->base.cols = in->base.cols;
    mtx->base.type = JMTX_TYPE_CDS;
    mtx->base.rows = in->base.rows;
    mtx->base.allocator_callbacks = *allocator_callbacks;

    *p_mtx = mtx;

    return JMTX_RESULT_SUCCESS;
}

/**
 * Creates a new CDS matrix with single precision from the real part of a complex CDS matrix with single precision.
 * Requires no memory allocation by reusing the memory of the initial matrix. Can not fail if the input matrix is valid.
 * @param in matrix which to convert (will be invalid if function succeeds)
 * @return converted matrix
 */
jmtx_matrix_cds *jmtx_matrix_cds_from_cfloat_real_inplace(jmtxc_matrix_cds *in)
{
    //  Convert diagonals
    if (in->main_diagonal)
    {
        float *const values = (float *)in->main_diagonal;
        const JMTX_INDEX_T len = jmtxc_matrix_cds_entries_in_dia(in, 0);
        for (JMTX_INDEX_T i = 0; i < len; ++i)
        {
            values[i] = crealf(in->main_diagonal[i]);
        }
    }

    for (JMTX_FAST_INT_T i = 0; i < in->sub_diagonals.count; ++i)
    {
        _Complex float *ptr = in->sub_diagonals.diagonals[i];
        float *const values = (float *)ptr;
        const JMTX_INDEX_T len = jmtxc_matrix_cds_entries_in_dia(in, -((int32_t)i));
        for (JMTX_INDEX_T j = 0; j < len; ++j)
        {
            values[j] = crealf(ptr[j]);
        }
    }

    for (JMTX_FAST_INT_T i = 0; i < in->super_diagonals.count; ++i)
    {
        _Complex float *ptr = in->super_diagonals.diagonals[i];
        float *const values = (float *)ptr;
        const JMTX_INDEX_T len = jmtxc_matrix_cds_entries_in_dia(in, +((int32_t)i));
        for (JMTX_INDEX_T j = 0; j < len; ++j)
        {
            values[j] = crealf(ptr[j]);
        }
    }

    //  Shrink the diagonals
    if (in->main_diagonal)
    {
        const JMTX_INDEX_T len = jmtxc_matrix_cds_entries_in_dia(in, 0);
        float *const new_ptr = in->base.allocator_callbacks.realloc(in->base.allocator_callbacks.state,
                                                                    in->main_diagonal, sizeof(*new_ptr) * len);
        if (new_ptr)
        {
            in->main_diagonal = (_Complex float *)new_ptr;
        }
    }

    for (JMTX_FAST_INT_T i = 0; i < in->sub_diagonals.count; ++i)
    {
        const JMTX_INDEX_T len = jmtxc_matrix_cds_entries_in_dia(in, -((int32_t)i));
        float *const new_ptr = in->base.allocator_callbacks.realloc(
            in->base.allocator_callbacks.state, in->sub_diagonals.diagonals[i], sizeof(*new_ptr) * len);
        if (new_ptr)
        {
            in->sub_diagonals.diagonals[i] = (_Complex float *)new_ptr;
        }
    }

    for (JMTX_FAST_INT_T i = 0; i < in->super_diagonals.count; ++i)
    {
        const JMTX_INDEX_T len = jmtxc_matrix_cds_entries_in_dia(in, +((int32_t)i));
        float *const new_ptr = in->base.allocator_callbacks.realloc(
            in->base.allocator_callbacks.state, in->super_diagonals.diagonals[i], sizeof(*new_ptr) * len);
        if (new_ptr)
        {
            in->super_diagonals.diagonals[i] = (_Complex float *)new_ptr;
        }
    }

    in->base.type = JMTX_TYPE_CDS;

    return (jmtx_matrix_cds *)in;
}

/**
 * Creates a new CDS matrix with single precision from the imaginary part of a complex CDS matrix with single precision.
 * Requires no memory allocation by reusing the memory of the initial matrix. Can not fail if the input matrix is valid.
 * @param in matrix which to convert (will be invalid if function succeeds)
 * @return converted matrix
 */
jmtx_matrix_cds *jmtx_matrix_cds_from_cfloat_imag_inplace(jmtxc_matrix_cds *in)
{
    //  Convert diagonals
    if (in->main_diagonal)
    {
        float *const values = (float *)in->main_diagonal;
        const JMTX_INDEX_T len = jmtxc_matrix_cds_entries_in_dia(in, 0);
        for (JMTX_INDEX_T i = 0; i < len; ++i)
        {
            values[i] = cimagf(in->main_diagonal[i]);
        }
    }

    for (JMTX_FAST_INT_T i = 0; i < in->sub_diagonals.count; ++i)
    {
        _Complex float *ptr = in->sub_diagonals.diagonals[i];
        float *const values = (float *)ptr;
        const JMTX_INDEX_T len = jmtxc_matrix_cds_entries_in_dia(in, -((int32_t)i));
        for (JMTX_INDEX_T j = 0; j < len; ++j)
        {
            values[j] = cimagf(ptr[j]);
        }
    }

    for (JMTX_FAST_INT_T i = 0; i < in->super_diagonals.count; ++i)
    {
        _Complex float *ptr = in->super_diagonals.diagonals[i];
        float *const values = (float *)ptr;
        const JMTX_INDEX_T len = jmtxc_matrix_cds_entries_in_dia(in, +((int32_t)i));
        for (JMTX_INDEX_T j = 0; j < len; ++j)
        {
            values[j] = cimagf(ptr[j]);
        }
    }

    //  Shrink the diagonals
    if (in->main_diagonal)
    {
        const JMTX_INDEX_T len = jmtxc_matrix_cds_entries_in_dia(in, 0);
        float *const new_ptr = in->base.allocator_callbacks.realloc(in->base.allocator_callbacks.state,
                                                                    in->main_diagonal, sizeof(*new_ptr) * len);
        if (new_ptr)
        {
            in->main_diagonal = (_Complex float *)new_ptr;
        }
    }

    for (JMTX_FAST_INT_T i = 0; i < in->sub_diagonals.count; ++i)
    {
        const JMTX_INDEX_T len = jmtxc_matrix_cds_entries_in_dia(in, -((int32_t)i));
        float *const new_ptr = in->base.allocator_callbacks.realloc(
            in->base.allocator_callbacks.state, in->sub_diagonals.diagonals[i], sizeof(*new_ptr) * len);
        if (new_ptr)
        {
            in->sub_diagonals.diagonals[i] = (_Complex float *)new_ptr;
        }
    }

    for (JMTX_FAST_INT_T i = 0; i < in->super_diagonals.count; ++i)
    {
        const JMTX_INDEX_T len = jmtxc_matrix_cds_entries_in_dia(in, +((int32_t)i));
        float *const new_ptr = in->base.allocator_callbacks.realloc(
            in->base.allocator_callbacks.state, in->super_diagonals.diagonals[i], sizeof(*new_ptr) * len);
        if (new_ptr)
        {
            in->super_diagonals.diagonals[i] = (_Complex float *)new_ptr;
        }
    }

    in->base.type = JMTX_TYPE_CDS;

    return (jmtx_matrix_cds *)in;
}

/**
 * Creates a new CDS matrix with single precision from the real part of a complex CDS matrix with single precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxs_matrix_cds_from_cfloat_real(jmtx_matrix_cds **p_mtx, const jmtxc_matrix_cds *in,
                                              const jmtx_allocator_callbacks *allocator_callbacks)
{
    if (!p_mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!in)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (allocator_callbacks &&
        (!allocator_callbacks->free || !allocator_callbacks->alloc || !allocator_callbacks->realloc))
    {
        return JMTX_RESULT_BAD_PARAM;
    }
    if (in->base.type != JMTXC_TYPE_CDS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }

    return jmtx_matrix_cds_from_cfloat_real(p_mtx, in, allocator_callbacks);
}

/**
 * Creates a new CDS matrix with single precision from the imaginary part of a complex CDS matrix with single precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxs_matrix_cds_from_cfloat_imag(jmtx_matrix_cds **p_mtx, const jmtxc_matrix_cds *in,
                                              const jmtx_allocator_callbacks *allocator_callbacks)
{
    if (!p_mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!in)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (allocator_callbacks &&
        (!allocator_callbacks->free || !allocator_callbacks->alloc || !allocator_callbacks->realloc))
    {
        return JMTX_RESULT_BAD_PARAM;
    }
    if (in->base.type != JMTXC_TYPE_CDS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }

    return jmtx_matrix_cds_from_cfloat_imag(p_mtx, in, allocator_callbacks);
}

/**
 * Creates a new complex CDS matrix with single precision from a CDS matrix with single precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in_real matrix to use as the real component of the new matrix (if NULL real part is zero)
 * @param in_imag matrix to use as the real component of the new matrix (if NULL imaginary part is zero)
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxc_matrix_cds_from_float(jmtxc_matrix_cds **p_mtx, const jmtx_matrix_cds *in_real,
                                        const jmtx_matrix_cds *in_imag,
                                        const jmtx_allocator_callbacks *allocator_callbacks)
{
    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }
    JMTX_INDEX_T r, c;
    assert(in_real || in_imag);
    if (in_real)
    {
        r = in_real->base.rows;
        c = in_real->base.cols;
    }
    else
    {
        r = in_imag->base.rows;
        c = in_imag->base.cols;
    }

    jmtxc_matrix_cds *mtx;
    jmtx_result res = jmtxc_matrix_cds_new(&mtx, r, c, 0, (int32_t[]){0}, allocator_callbacks);
    if (res != JMTX_RESULT_SUCCESS)
    {
        return res;
    }

    if (in_real)
    {
        for (JMTX_FAST_INT_T i = 0; i < in_real->sub_diagonals.count; ++i)
        {
            const int32_t dia_idx = -(int32_t)in_real->sub_diagonals.indices[i];
            const float *const dia_ptr = in_real->sub_diagonals.diagonals[i];
            JMTX_INDEX_T len;
            _Complex float *const ptr = jmtxc_matrix_cds_allocate_zero_diagonal(mtx, dia_idx, &len);
            for (JMTX_FAST_INT_T j = 0; j < len; ++j)
            {
                ptr[j] += (_Complex float)dia_ptr[j];
            }
        }
        if (in_real->main_diagonal)
        {
            JMTX_INDEX_T len;
            _Complex float *const ptr = jmtxc_matrix_cds_allocate_zero_diagonal(mtx, 0, &len);
            for (JMTX_FAST_INT_T j = 0; j < len; ++j)
            {
                ptr[j] += (_Complex float)in_real->main_diagonal[j];
            }
        }
        for (JMTX_FAST_INT_T i = 0; i < in_real->super_diagonals.count; ++i)
        {
            const int32_t dia_idx = (int32_t)in_real->super_diagonals.indices[i];
            const float *const dia_ptr = in_real->super_diagonals.diagonals[i];
            JMTX_INDEX_T len;
            _Complex float *const ptr = jmtxc_matrix_cds_allocate_zero_diagonal(mtx, dia_idx, &len);
            for (JMTX_FAST_INT_T j = 0; j < len; ++j)
            {
                ptr[j] += (_Complex float)dia_ptr[j];
            }
        }
    }
    if (in_imag)
    {
        for (JMTX_FAST_INT_T i = 0; i < in_imag->sub_diagonals.count; ++i)
        {
            const int32_t dia_idx = -(int32_t)in_imag->sub_diagonals.indices[i];
            const float *const dia_ptr = in_imag->sub_diagonals.diagonals[i];
            JMTX_INDEX_T len;
            _Complex float *const ptr = jmtxc_matrix_cds_allocate_zero_diagonal(mtx, dia_idx, &len);
            for (JMTX_FAST_INT_T j = 0; j < len; ++j)
            {
                ptr[j] += (_Complex float)dia_ptr[j] * _Complex_I;
            }
        }
        if (in_imag->main_diagonal)
        {
            JMTX_INDEX_T len;
            _Complex float *const ptr = jmtxc_matrix_cds_allocate_zero_diagonal(mtx, 0, &len);
            for (JMTX_FAST_INT_T j = 0; j < len; ++j)
            {
                ptr[j] += (_Complex float)in_imag->main_diagonal[j] * _Complex_I;
            }
        }
        for (JMTX_FAST_INT_T i = 0; i < in_imag->super_diagonals.count; ++i)
        {
            const int32_t dia_idx = (int32_t)in_imag->super_diagonals.indices[i];
            const float *const dia_ptr = in_imag->super_diagonals.diagonals[i];
            JMTX_INDEX_T len;
            _Complex float *const ptr = jmtxc_matrix_cds_allocate_zero_diagonal(mtx, dia_idx, &len);
            for (JMTX_FAST_INT_T j = 0; j < len; ++j)
            {
                ptr[j] += (_Complex float)dia_ptr[j] * _Complex_I;
            }
        }
    }

    mtx->base.cols = c;
    mtx->base.type = JMTXC_TYPE_CDS;
    mtx->base.rows = r;
    mtx->base.allocator_callbacks = *allocator_callbacks;

    *p_mtx = mtx;

    return JMTX_RESULT_SUCCESS;
}

/**
 * Creates a new complex CDS matrix with single precision from a CDS matrix with single precision as its real part.
 * Only one memory reallocation.
 * may be needed. Can not fail if the input matrix is valid.
 * @param in matrix which to convert
 * @return converted matrix, or NULL in case of allocation failure
 */
jmtxc_matrix_cds *jmtxc_matrix_cds_from_float_real_inplace(jmtx_matrix_cds *in)
{
    //  Expand the diagonals
    if (in->main_diagonal)
    {
        const JMTX_INDEX_T len = jmtx_matrix_cds_entries_in_dia(in, 0);
        _Complex float *const new_ptr = in->base.allocator_callbacks.realloc(in->base.allocator_callbacks.state,
                                                                             in->main_diagonal, sizeof(*new_ptr) * len);
        if (new_ptr)
        {
            in->main_diagonal = (float *)new_ptr;
        }
    }

    for (JMTX_FAST_INT_T i = 0; i < in->sub_diagonals.count; ++i)
    {
        const JMTX_INDEX_T len = jmtx_matrix_cds_entries_in_dia(in, -((int32_t)i));
        _Complex float *const new_ptr = in->base.allocator_callbacks.realloc(
            in->base.allocator_callbacks.state, in->sub_diagonals.diagonals[i], sizeof(*new_ptr) * len);
        if (new_ptr)
        {
            in->sub_diagonals.diagonals[i] = (float *)new_ptr;
        }
    }

    for (JMTX_FAST_INT_T i = 0; i < in->super_diagonals.count; ++i)
    {
        const JMTX_INDEX_T len = jmtx_matrix_cds_entries_in_dia(in, +((int32_t)i));
        _Complex float *const new_ptr = in->base.allocator_callbacks.realloc(
            in->base.allocator_callbacks.state, in->super_diagonals.diagonals[i], sizeof(*new_ptr) * len);
        if (new_ptr)
        {
            in->super_diagonals.diagonals[i] = (float *)new_ptr;
        }
    }

    //  Convert diagonals
    if (in->main_diagonal)
    {
        _Complex float *const values = (_Complex float *)in->main_diagonal;
        const JMTX_INDEX_T len = jmtx_matrix_cds_entries_in_dia(in, 0);
        for (JMTX_INDEX_T i = 0; i < len; ++i)
        {
            values[i] = in->main_diagonal[i];
        }
    }

    for (JMTX_FAST_INT_T i = 0; i < in->sub_diagonals.count; ++i)
    {
        float *ptr = in->sub_diagonals.diagonals[i];
        _Complex float *const values = (_Complex float *)ptr;
        const JMTX_INDEX_T len = jmtx_matrix_cds_entries_in_dia(in, -((int32_t)i));
        for (JMTX_INDEX_T j = 0; j < len; ++j)
        {
            values[j] = ptr[j];
        }
    }

    for (JMTX_FAST_INT_T i = 0; i < in->super_diagonals.count; ++i)
    {
        float *ptr = in->super_diagonals.diagonals[i];
        _Complex float *const values = (_Complex float *)ptr;
        const JMTX_INDEX_T len = jmtx_matrix_cds_entries_in_dia(in, +((int32_t)i));
        for (JMTX_INDEX_T j = 0; j < len; ++j)
        {
            values[j] = (ptr[j]);
        }
    }

    in->base.type = JMTXC_TYPE_CDS;

    return (jmtxc_matrix_cds *)in;
}

/**
 * Creates a new CDS matrix with single precision from a CDS matrix with single precision as its imaginary part.
 * Only one memory reallocation.
 * may be needed. Can not fail if the input matrix is valid.
 * @param in matrix which to convert
 * @return converted matrix, or NULL in case of allocation failure
 */
jmtxc_matrix_cds *jmtxc_matrix_cds_from_float_imag_inplace(jmtx_matrix_cds *in)
{
    //  Expand the diagonals
    if (in->main_diagonal)
    {
        const JMTX_INDEX_T len = jmtx_matrix_cds_entries_in_dia(in, 0);
        _Complex float *const new_ptr = in->base.allocator_callbacks.realloc(in->base.allocator_callbacks.state,
                                                                             in->main_diagonal, sizeof(*new_ptr) * len);
        if (new_ptr)
        {
            in->main_diagonal = (float *)new_ptr;
        }
    }

    for (JMTX_FAST_INT_T i = 0; i < in->sub_diagonals.count; ++i)
    {
        const JMTX_INDEX_T len = jmtx_matrix_cds_entries_in_dia(in, -((int32_t)i));
        _Complex float *const new_ptr = in->base.allocator_callbacks.realloc(
            in->base.allocator_callbacks.state, in->sub_diagonals.diagonals[i], sizeof(*new_ptr) * len);
        if (new_ptr)
        {
            in->sub_diagonals.diagonals[i] = (float *)new_ptr;
        }
    }

    for (JMTX_FAST_INT_T i = 0; i < in->super_diagonals.count; ++i)
    {
        const JMTX_INDEX_T len = jmtx_matrix_cds_entries_in_dia(in, +((int32_t)i));
        _Complex float *const new_ptr = in->base.allocator_callbacks.realloc(
            in->base.allocator_callbacks.state, in->super_diagonals.diagonals[i], sizeof(*new_ptr) * len);
        if (new_ptr)
        {
            in->super_diagonals.diagonals[i] = (float *)new_ptr;
        }
    }

    //  Convert diagonals
    if (in->main_diagonal)
    {
        _Complex float *const values = (_Complex float *)in->main_diagonal;
        const JMTX_INDEX_T len = jmtx_matrix_cds_entries_in_dia(in, 0);
        for (JMTX_INDEX_T i = 0; i < len; ++i)
        {
            values[i] = in->main_diagonal[i] * _Complex_I;
        }
    }

    for (JMTX_FAST_INT_T i = 0; i < in->sub_diagonals.count; ++i)
    {
        float *ptr = in->sub_diagonals.diagonals[i];
        _Complex float *const values = (_Complex float *)ptr;
        const JMTX_INDEX_T len = jmtx_matrix_cds_entries_in_dia(in, -((int32_t)i));
        for (JMTX_INDEX_T j = 0; j < len; ++j)
        {
            values[j] = ptr[j] * _Complex_I;
        }
    }

    for (JMTX_FAST_INT_T i = 0; i < in->super_diagonals.count; ++i)
    {
        float *ptr = in->super_diagonals.diagonals[i];
        _Complex float *const values = (_Complex float *)ptr;
        const JMTX_INDEX_T len = jmtx_matrix_cds_entries_in_dia(in, +((int32_t)i));
        for (JMTX_INDEX_T j = 0; j < len; ++j)
        {
            values[j] = (ptr[j]) * _Complex_I;
        }
    }

    in->base.type = JMTXC_TYPE_CDS;

    return (jmtxc_matrix_cds *)in;
}

/**
 * Creates a new complex CDS matrix with single precision from a CDS matrix with single precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in_real matrix to use as the real component of the new matrix (if NULL real part is zero)
 * @param in_imag matrix to use as the real component of the new matrix (if NULL imaginary part is zero)
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxcs_matrix_cds_from_float(jmtxc_matrix_cds **p_mtx, const jmtx_matrix_cds *in_real,
                                         const jmtx_matrix_cds *in_imag,
                                         const jmtx_allocator_callbacks *allocator_callbacks)
{
    if (!p_mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!in_real && !in_imag)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (allocator_callbacks &&
        (!allocator_callbacks->free || !allocator_callbacks->alloc || !allocator_callbacks->realloc))
    {
        return JMTX_RESULT_BAD_PARAM;
    }
    if ((in_real && in_real->base.type != JMTX_TYPE_CDS) || (in_imag && in_imag->base.type != JMTX_TYPE_CDS))
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if ((in_real && in_imag) && (in_real->base.rows != in_imag->base.rows || in_real->base.cols != in_imag->base.cols))
    {
        return JMTX_RESULT_BAD_MATRIX;
    }

    return jmtxc_matrix_cds_from_float(p_mtx, in_real, in_imag, allocator_callbacks);
}

/***********************************************************************************************************************
 *                                                                                                                     *
 *                                          DOUBLE <-> COMPLEX DOUBLE                                                  *
 *                                                                                                                     *
 **********************************************************************************************************************/

/**
 * Creates a new CDS matrix with JMTX_SCALAR_T precision from the real part of a complex CDS matrix with JMTX_SCALAR_T precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result JMTX_NAME_TYPED(matrix_cds_from_cdouble_real(JMTX_NAME_TYPED(matrix_cds) **p_mtx, const jmtxz_matrix_cds *in,
                                               const jmtx_allocator_callbacks *allocator_callbacks)
{
    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }
    JMTX_NAME_TYPED(matrix_cds) * mtx;
    jmtx_result res = JMTX_NAME_TYPED(matrix_cds_new(&mtx, in->base.rows, in->base.cols, 0, (int32_t[]){0}, allocator_callbacks);
    if (res != JMTX_RESULT_SUCCESS)
    {
        return res;
    }

    for (JMTX_FAST_INT_T i = 0; i < in->sub_diagonals.count; ++i)
    {
        const int32_t dia_idx = -(int32_t)in->sub_diagonals.indices[i];
        const _Complex JMTX_SCALAR_T *const dia_ptr = in->sub_diagonals.diagonals[i];
        JMTX_INDEX_T len;
        JMTX_SCALAR_T *const ptr = JMTX_NAME_TYPED(matrix_cds_allocate_diagonal(mtx, dia_idx, &len);
        for (JMTX_FAST_INT_T j = 0; j < len; ++j)
        {
            ptr[j] = creal(dia_ptr[j]);
        }
    }
    if (in->main_diagonal)
    {
        JMTX_INDEX_T len;
        JMTX_SCALAR_T *const ptr = JMTX_NAME_TYPED(matrix_cds_allocate_diagonal(mtx, 0, &len);
        for (JMTX_FAST_INT_T j = 0; j < len; ++j)
        {
            ptr[j] = creal(in->main_diagonal[j]);
        }
    }
    for (JMTX_FAST_INT_T i = 0; i < in->super_diagonals.count; ++i)
    {
        const int32_t dia_idx = (int32_t)in->super_diagonals.indices[i];
        const _Complex JMTX_SCALAR_T *const dia_ptr = in->super_diagonals.diagonals[i];
        JMTX_INDEX_T len;
        JMTX_SCALAR_T *const ptr = JMTX_NAME_TYPED(matrix_cds_allocate_diagonal(mtx, dia_idx, &len);
        for (JMTX_FAST_INT_T j = 0; j < len; ++j)
        {
            ptr[j] = creal(dia_ptr[j]);
        }
    }

    mtx->base.cols = in->base.cols;
    mtx->base.type = JMTXD_TYPE_CDS;
    mtx->base.rows = in->base.rows;
    mtx->base.allocator_callbacks = *allocator_callbacks;

    *p_mtx = mtx;

    return JMTX_RESULT_SUCCESS;
}
/**
 * Creates a new CDS matrix with JMTX_SCALAR_T precision from the imaginary part of a complex CDS matrix with JMTX_SCALAR_T precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result JMTX_NAME_TYPED(matrix_cds_from_cdouble_imag(JMTX_NAME_TYPED(matrix_cds) **p_mtx, const jmtxz_matrix_cds *in,
                                               const jmtx_allocator_callbacks *allocator_callbacks)
{
    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }
    JMTX_NAME_TYPED(matrix_cds) * mtx;
    jmtx_result res = JMTX_NAME_TYPED(matrix_cds_new(&mtx, in->base.rows, in->base.cols, 0, (int32_t[]){0}, allocator_callbacks);
    if (res != JMTX_RESULT_SUCCESS)
    {
        return res;
    }

    for (JMTX_FAST_INT_T i = 0; i < in->sub_diagonals.count; ++i)
    {
        const int32_t dia_idx = -(int32_t)in->sub_diagonals.indices[i];
        const _Complex JMTX_SCALAR_T *const dia_ptr = in->sub_diagonals.diagonals[i];
        JMTX_INDEX_T len;
        JMTX_SCALAR_T *const ptr = JMTX_NAME_TYPED(matrix_cds_allocate_diagonal(mtx, dia_idx, &len);
        for (JMTX_FAST_INT_T j = 0; j < len; ++j)
        {
            ptr[j] = cimag(dia_ptr[j]);
        }
    }
    if (in->main_diagonal)
    {
        JMTX_INDEX_T len;
        JMTX_SCALAR_T *const ptr = JMTX_NAME_TYPED(matrix_cds_allocate_diagonal(mtx, 0, &len);
        for (JMTX_FAST_INT_T j = 0; j < len; ++j)
        {
            ptr[j] = cimag(in->main_diagonal[j]);
        }
    }
    for (JMTX_FAST_INT_T i = 0; i < in->super_diagonals.count; ++i)
    {
        const int32_t dia_idx = (int32_t)in->super_diagonals.indices[i];
        const _Complex JMTX_SCALAR_T *const dia_ptr = in->super_diagonals.diagonals[i];
        JMTX_INDEX_T len;
        JMTX_SCALAR_T *const ptr = JMTX_NAME_TYPED(matrix_cds_allocate_diagonal(mtx, dia_idx, &len);
        for (JMTX_FAST_INT_T j = 0; j < len; ++j)
        {
            ptr[j] = cimag(dia_ptr[j]);
        }
    }

    mtx->base.cols = in->base.cols;
    mtx->base.type = JMTXD_TYPE_CDS;
    mtx->base.rows = in->base.rows;
    mtx->base.allocator_callbacks = *allocator_callbacks;

    *p_mtx = mtx;

    return JMTX_RESULT_SUCCESS;
}

/**
 * Creates a new CDS matrix with JMTX_SCALAR_T precision from the real part of a complex CDS matrix with JMTX_SCALAR_T precision.
 * Requires no memory allocation by reusing the memory of the initial matrix. Can not fail if the input matrix is valid.
 * @param in matrix which to convert (will be invalid if function succeeds)
 * @return converted matrix
 */
JMTX_NAME_TYPED(matrix_cds) *JMTX_NAME_TYPED(matrix_cds_from_cdouble_real_inplace(jmtxz_matrix_cds *in)
{
    //  Convert diagonals
    if (in->main_diagonal)
    {
        JMTX_SCALAR_T *const values = (JMTX_SCALAR_T *)in->main_diagonal;
        const JMTX_INDEX_T len = jmtxz_matrix_cds_entries_in_dia(in, 0);
        for (JMTX_INDEX_T i = 0; i < len; ++i)
        {
            values[i] = creal(in->main_diagonal[i]);
        }
    }

    for (JMTX_FAST_INT_T i = 0; i < in->sub_diagonals.count; ++i)
    {
        _Complex JMTX_SCALAR_T *ptr = in->sub_diagonals.diagonals[i];
        JMTX_SCALAR_T *const values = (JMTX_SCALAR_T *)ptr;
        const JMTX_INDEX_T len = jmtxz_matrix_cds_entries_in_dia(in, -((int32_t)i));
        for (JMTX_INDEX_T j = 0; j < len; ++j)
        {
            values[j] = creal(ptr[j]);
        }
    }

    for (JMTX_FAST_INT_T i = 0; i < in->super_diagonals.count; ++i)
    {
        _Complex JMTX_SCALAR_T *ptr = in->super_diagonals.diagonals[i];
        JMTX_SCALAR_T *const values = (JMTX_SCALAR_T *)ptr;
        const JMTX_INDEX_T len = jmtxz_matrix_cds_entries_in_dia(in, +((int32_t)i));
        for (JMTX_INDEX_T j = 0; j < len; ++j)
        {
            values[j] = creal(ptr[j]);
        }
    }

    //  Shrink the diagonals
    if (in->main_diagonal)
    {
        const JMTX_INDEX_T len = jmtxz_matrix_cds_entries_in_dia(in, 0);
        JMTX_SCALAR_T *const new_ptr = in->base.allocator_callbacks.realloc(in->base.allocator_callbacks.state,
                                                                            in->main_diagonal, sizeof(*new_ptr) * len);
        if (new_ptr)
        {
            in->main_diagonal = (_Complex JMTX_SCALAR_T *)new_ptr;
        }
    }

    for (JMTX_FAST_INT_T i = 0; i < in->sub_diagonals.count; ++i)
    {
        const JMTX_INDEX_T len = jmtxz_matrix_cds_entries_in_dia(in, -((int32_t)i));
        JMTX_SCALAR_T *const new_ptr = in->base.allocator_callbacks.realloc(
            in->base.allocator_callbacks.state, in->sub_diagonals.diagonals[i], sizeof(*new_ptr) * len);
        if (new_ptr)
        {
            in->sub_diagonals.diagonals[i] = (_Complex JMTX_SCALAR_T *)new_ptr;
        }
    }

    for (JMTX_FAST_INT_T i = 0; i < in->super_diagonals.count; ++i)
    {
        const JMTX_INDEX_T len = jmtxz_matrix_cds_entries_in_dia(in, +((int32_t)i));
        JMTX_SCALAR_T *const new_ptr = in->base.allocator_callbacks.realloc(
            in->base.allocator_callbacks.state, in->super_diagonals.diagonals[i], sizeof(*new_ptr) * len);
        if (new_ptr)
        {
            in->super_diagonals.diagonals[i] = (_Complex JMTX_SCALAR_T *)new_ptr;
        }
    }

    in->base.type = JMTXD_TYPE_CDS;

    return (JMTX_NAME_TYPED(matrix_cds) *)in;
}

/**
 * Creates a new CDS matrix with JMTX_SCALAR_T precision from the imaginary part of a complex CDS matrix with JMTX_SCALAR_T precision.
 * Requires no memory allocation by reusing the memory of the initial matrix. Can not fail if the input matrix is valid.
 * @param in matrix which to convert (will be invalid if function succeeds)
 * @return converted matrix
 */
JMTX_NAME_TYPED(matrix_cds) *JMTX_NAME_TYPED(matrix_cds_from_cdouble_imag_inplace(jmtxz_matrix_cds *in)
{
    //  Convert diagonals
    if (in->main_diagonal)
    {
        JMTX_SCALAR_T *const values = (JMTX_SCALAR_T *)in->main_diagonal;
        const JMTX_INDEX_T len = jmtxz_matrix_cds_entries_in_dia(in, 0);
        for (JMTX_INDEX_T i = 0; i < len; ++i)
        {
            values[i] = cimag(in->main_diagonal[i]);
        }
    }

    for (JMTX_FAST_INT_T i = 0; i < in->sub_diagonals.count; ++i)
    {
        _Complex JMTX_SCALAR_T *ptr = in->sub_diagonals.diagonals[i];
        JMTX_SCALAR_T *const values = (JMTX_SCALAR_T *)ptr;
        const JMTX_INDEX_T len = jmtxz_matrix_cds_entries_in_dia(in, -((int32_t)i));
        for (JMTX_INDEX_T j = 0; j < len; ++j)
        {
            values[j] = cimag(ptr[j]);
        }
    }

    for (JMTX_FAST_INT_T i = 0; i < in->super_diagonals.count; ++i)
    {
        _Complex JMTX_SCALAR_T *ptr = in->super_diagonals.diagonals[i];
        JMTX_SCALAR_T *const values = (JMTX_SCALAR_T *)ptr;
        const JMTX_INDEX_T len = jmtxz_matrix_cds_entries_in_dia(in, +((int32_t)i));
        for (JMTX_INDEX_T j = 0; j < len; ++j)
        {
            values[j] = cimag(ptr[j]);
        }
    }

    //  Shrink the diagonals
    if (in->main_diagonal)
    {
        const JMTX_INDEX_T len = jmtxz_matrix_cds_entries_in_dia(in, 0);
        JMTX_SCALAR_T *const new_ptr = in->base.allocator_callbacks.realloc(in->base.allocator_callbacks.state,
                                                                            in->main_diagonal, sizeof(*new_ptr) * len);
        if (new_ptr)
        {
            in->main_diagonal = (_Complex JMTX_SCALAR_T *)new_ptr;
        }
    }

    for (JMTX_FAST_INT_T i = 0; i < in->sub_diagonals.count; ++i)
    {
        const JMTX_INDEX_T len = jmtxz_matrix_cds_entries_in_dia(in, -((int32_t)i));
        JMTX_SCALAR_T *const new_ptr = in->base.allocator_callbacks.realloc(
            in->base.allocator_callbacks.state, in->sub_diagonals.diagonals[i], sizeof(*new_ptr) * len);
        if (new_ptr)
        {
            in->sub_diagonals.diagonals[i] = (_Complex JMTX_SCALAR_T *)new_ptr;
        }
    }

    for (JMTX_FAST_INT_T i = 0; i < in->super_diagonals.count; ++i)
    {
        const JMTX_INDEX_T len = jmtxz_matrix_cds_entries_in_dia(in, +((int32_t)i));
        JMTX_SCALAR_T *const new_ptr = in->base.allocator_callbacks.realloc(
            in->base.allocator_callbacks.state, in->super_diagonals.diagonals[i], sizeof(*new_ptr) * len);
        if (new_ptr)
        {
            in->super_diagonals.diagonals[i] = (_Complex JMTX_SCALAR_T *)new_ptr;
        }
    }

    in->base.type = JMTXD_TYPE_CDS;

    return (JMTX_NAME_TYPED(matrix_cds) *)in;
}

/**
 * Creates a new CDS matrix with JMTX_SCALAR_T precision from the real part of a complex CDS matrix with JMTX_SCALAR_T precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result .,.,.,.,..matrix_cds_from_cdouble_real(JMTX_NAME_TYPED(matrix_cds) **p_mtx, const jmtxz_matrix_cds *in,
                                                const jmtx_allocator_callbacks *allocator_callbacks)
{
    if (!p_mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!in)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (allocator_callbacks &&
        (!allocator_callbacks->free || !allocator_callbacks->alloc || !allocator_callbacks->realloc))
    {
        return JMTX_RESULT_BAD_PARAM;
    }
    if (in->base.type != JMTXZ_TYPE_CDS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }

    return JMTX_NAME_TYPED(matrix_cds_from_cdouble_real(p_mtx, in, allocator_callbacks);
}

/**
 * Creates a new CDS matrix with JMTX_SCALAR_T precision from the imaginary part of a complex CDS matrix with JMTX_SCALAR_T precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result .,.,.,.,..matrix_cds_from_cdouble_imag(JMTX_NAME_TYPED(matrix_cds) **p_mtx, const jmtxz_matrix_cds *in,
                                                const jmtx_allocator_callbacks *allocator_callbacks)
{
    if (!p_mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!in)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (allocator_callbacks &&
        (!allocator_callbacks->free || !allocator_callbacks->alloc || !allocator_callbacks->realloc))
    {
        return JMTX_RESULT_BAD_PARAM;
    }
    if (in->base.type != JMTXZ_TYPE_CDS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }

    return JMTX_NAME_TYPED(matrix_cds_from_cdouble_imag(p_mtx, in, allocator_callbacks);
}

/**
 * Creates a new complex CDS matrix with JMTX_SCALAR_T precision from a CDS matrix with JMTX_SCALAR_T precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in_real matrix to use as the real component of the new matrix (if NULL real part is zero)
 * @param in_imag matrix to use as the real component of the new matrix (if NULL imaginary part is zero)
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxz_matrix_cds_from_double(jmtxz_matrix_cds **p_mtx, const JMTX_NAME_TYPED(matrix_cds) *in_real,
                                         const JMTX_NAME_TYPED(matrix_cds) *in_imag,
                                         const jmtx_allocator_callbacks *allocator_callbacks)
{
    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }
    JMTX_INDEX_T r, c;
    assert(in_real || in_imag);
    if (in_real)
    {
        r = in_real->base.rows;
        c = in_real->base.cols;
    }
    else
    {
        r = in_imag->base.rows;
        c = in_imag->base.cols;
    }

    jmtxz_matrix_cds *mtx;
    jmtx_result res = jmtxz_matrix_cds_new(&mtx, r, c, 0, (int32_t[]){0}, allocator_callbacks);
    if (res != JMTX_RESULT_SUCCESS)
    {
        return res;
    }

    if (in_real)
    {
        for (JMTX_FAST_INT_T i = 0; i < in_real->sub_diagonals.count; ++i)
        {
            const int32_t dia_idx = -(int32_t)in_real->sub_diagonals.indices[i];
            const JMTX_SCALAR_T *const dia_ptr = in_real->sub_diagonals.diagonals[i];
            JMTX_INDEX_T len;
            _Complex JMTX_SCALAR_T *const ptr = jmtxz_matrix_cds_allocate_zero_diagonal(mtx, dia_idx, &len);
            for (JMTX_FAST_INT_T j = 0; j < len; ++j)
            {
                ptr[j] += (_Complex double)dia_ptr[j];
            }
        }
        if (in_real->main_diagonal)
        {
            JMTX_INDEX_T len;
            _Complex JMTX_SCALAR_T *const ptr = jmtxz_matrix_cds_allocate_zero_diagonal(mtx, 0, &len);
            for (JMTX_FAST_INT_T j = 0; j < len; ++j)
            {
                ptr[j] += (_Complex double)in_real->main_diagonal[j];
            }
        }
        for (JMTX_FAST_INT_T i = 0; i < in_real->super_diagonals.count; ++i)
        {
            const int32_t dia_idx = (int32_t)in_real->super_diagonals.indices[i];
            const JMTX_SCALAR_T *const dia_ptr = in_real->super_diagonals.diagonals[i];
            JMTX_INDEX_T len;
            _Complex JMTX_SCALAR_T *const ptr = jmtxz_matrix_cds_allocate_zero_diagonal(mtx, dia_idx, &len);
            for (JMTX_FAST_INT_T j = 0; j < len; ++j)
            {
                ptr[j] += (_Complex double)dia_ptr[j];
            }
        }
    }
    if (in_imag)
    {
        for (JMTX_FAST_INT_T i = 0; i < in_imag->sub_diagonals.count; ++i)
        {
            const int32_t dia_idx = -(int32_t)in_imag->sub_diagonals.indices[i];
            const JMTX_SCALAR_T *const dia_ptr = in_imag->sub_diagonals.diagonals[i];
            JMTX_INDEX_T len;
            _Complex JMTX_SCALAR_T *const ptr = jmtxz_matrix_cds_allocate_zero_diagonal(mtx, dia_idx, &len);
            for (JMTX_FAST_INT_T j = 0; j < len; ++j)
            {
                ptr[j] += (_Complex double)dia_ptr[j] * _Complex_I;
            }
        }
        if (in_imag->main_diagonal)
        {
            JMTX_INDEX_T len;
            _Complex JMTX_SCALAR_T *const ptr = jmtxz_matrix_cds_allocate_zero_diagonal(mtx, 0, &len);
            for (JMTX_FAST_INT_T j = 0; j < len; ++j)
            {
                ptr[j] += (_Complex double)in_imag->main_diagonal[j] * _Complex_I;
            }
        }
        for (JMTX_FAST_INT_T i = 0; i < in_imag->super_diagonals.count; ++i)
        {
            const int32_t dia_idx = (int32_t)in_imag->super_diagonals.indices[i];
            const JMTX_SCALAR_T *const dia_ptr = in_imag->super_diagonals.diagonals[i];
            JMTX_INDEX_T len;
            _Complex JMTX_SCALAR_T *const ptr = jmtxz_matrix_cds_allocate_zero_diagonal(mtx, dia_idx, &len);
            for (JMTX_FAST_INT_T j = 0; j < len; ++j)
            {
                ptr[j] += (_Complex double)dia_ptr[j] * _Complex_I;
            }
        }
    }

    mtx->base.cols = c;
    mtx->base.type = JMTXZ_TYPE_CDS;
    mtx->base.rows = r;
    mtx->base.allocator_callbacks = *allocator_callbacks;

    *p_mtx = mtx;

    return JMTX_RESULT_SUCCESS;
}

/**
 * Creates a new complex CDS matrix with JMTX_SCALAR_T precision from a CDS matrix with JMTX_SCALAR_T precision as its real part.
 * Only one memory reallocation.
 * may be needed. Can not fail if the input matrix is valid.
 * @param in matrix which to convert
 * @return converted matrix, or NULL in case of allocation failure
 */
jmtxz_matrix_cds *jmtxz_matrix_cds_from_double_real_inplace(JMTX_NAME_TYPED(matrix_cds) *in)
{
    //  Expand the diagonals
    if (in->main_diagonal)
    {
        const JMTX_INDEX_T len = JMTX_NAME_TYPED(matrix_cds_entries_in_dia(in, 0);
        _Complex JMTX_SCALAR_T *const new_ptr = in->base.allocator_callbacks.realloc(
            in->base.allocator_callbacks.state, in->main_diagonal, sizeof(*new_ptr) * len);
        if (new_ptr)
        {
            in->main_diagonal = (JMTX_SCALAR_T *)new_ptr;
        }
    }

    for (JMTX_FAST_INT_T i = 0; i < in->sub_diagonals.count; ++i)
    {
        const JMTX_INDEX_T len = JMTX_NAME_TYPED(matrix_cds_entries_in_dia(in, -((int32_t)i));
        _Complex JMTX_SCALAR_T *const new_ptr = in->base.allocator_callbacks.realloc(
            in->base.allocator_callbacks.state, in->sub_diagonals.diagonals[i], sizeof(*new_ptr) * len);
        if (new_ptr)
        {
            in->sub_diagonals.diagonals[i] = (JMTX_SCALAR_T *)new_ptr;
        }
    }

    for (JMTX_FAST_INT_T i = 0; i < in->super_diagonals.count; ++i)
    {
        const JMTX_INDEX_T len = JMTX_NAME_TYPED(matrix_cds_entries_in_dia(in, +((int32_t)i));
        _Complex JMTX_SCALAR_T *const new_ptr = in->base.allocator_callbacks.realloc(
            in->base.allocator_callbacks.state, in->super_diagonals.diagonals[i], sizeof(*new_ptr) * len);
        if (new_ptr)
        {
            in->super_diagonals.diagonals[i] = (JMTX_SCALAR_T *)new_ptr;
        }
    }

    //  Convert diagonals
    if (in->main_diagonal)
    {
        _Complex JMTX_SCALAR_T *const values = (_Complex JMTX_SCALAR_T *)in->main_diagonal;
        const JMTX_INDEX_T len = JMTX_NAME_TYPED(matrix_cds_entries_in_dia(in, 0);
        for (JMTX_INDEX_T i = 0; i < len; ++i)
        {
            values[i] = in->main_diagonal[i];
        }
    }

    for (JMTX_FAST_INT_T i = 0; i < in->sub_diagonals.count; ++i)
    {
        JMTX_SCALAR_T *ptr = in->sub_diagonals.diagonals[i];
        _Complex JMTX_SCALAR_T *const values = (_Complex JMTX_SCALAR_T *)ptr;
        const JMTX_INDEX_T len = JMTX_NAME_TYPED(matrix_cds_entries_in_dia(in, -((int32_t)i));
        for (JMTX_INDEX_T j = 0; j < len; ++j)
        {
            values[j] = ptr[j];
        }
    }

    for (JMTX_FAST_INT_T i = 0; i < in->super_diagonals.count; ++i)
    {
        JMTX_SCALAR_T *ptr = in->super_diagonals.diagonals[i];
        _Complex JMTX_SCALAR_T *const values = (_Complex JMTX_SCALAR_T *)ptr;
        const JMTX_INDEX_T len = JMTX_NAME_TYPED(matrix_cds_entries_in_dia(in, +((int32_t)i));
        for (JMTX_INDEX_T j = 0; j < len; ++j)
        {
            values[j] = (ptr[j]);
        }
    }

    in->base.type = JMTXZ_TYPE_CDS;

    return (jmtxz_matrix_cds *)in;
}

/**
 * Creates a new CDS matrix with JMTX_SCALAR_T precision from a CDS matrix with JMTX_SCALAR_T precision as its imaginary part.
 * Only one memory reallocation.
 * may be needed. Can not fail if the input matrix is valid.
 * @param in matrix which to convert
 * @return converted matrix, or NULL in case of allocation failure
 */
jmtxz_matrix_cds *jmtxz_matrix_cds_from_double_imag_inplace(JMTX_NAME_TYPED(matrix_cds) *in)
{
    //  Expand the diagonals
    if (in->main_diagonal)
    {
        const JMTX_INDEX_T len = JMTX_NAME_TYPED(matrix_cds_entries_in_dia(in, 0);
        _Complex JMTX_SCALAR_T *const new_ptr = in->base.allocator_callbacks.realloc(
            in->base.allocator_callbacks.state, in->main_diagonal, sizeof(*new_ptr) * len);
        if (new_ptr)
        {
            in->main_diagonal = (JMTX_SCALAR_T *)new_ptr;
        }
    }

    for (JMTX_FAST_INT_T i = 0; i < in->sub_diagonals.count; ++i)
    {
        const JMTX_INDEX_T len = JMTX_NAME_TYPED(matrix_cds_entries_in_dia(in, -((int32_t)i));
        _Complex JMTX_SCALAR_T *const new_ptr = in->base.allocator_callbacks.realloc(
            in->base.allocator_callbacks.state, in->sub_diagonals.diagonals[i], sizeof(*new_ptr) * len);
        if (new_ptr)
        {
            in->sub_diagonals.diagonals[i] = (JMTX_SCALAR_T *)new_ptr;
        }
    }

    for (JMTX_FAST_INT_T i = 0; i < in->super_diagonals.count; ++i)
    {
        const JMTX_INDEX_T len = JMTX_NAME_TYPED(matrix_cds_entries_in_dia(in, +((int32_t)i));
        _Complex JMTX_SCALAR_T *const new_ptr = in->base.allocator_callbacks.realloc(
            in->base.allocator_callbacks.state, in->super_diagonals.diagonals[i], sizeof(*new_ptr) * len);
        if (new_ptr)
        {
            in->super_diagonals.diagonals[i] = (JMTX_SCALAR_T *)new_ptr;
        }
    }

    //  Convert diagonals
    if (in->main_diagonal)
    {
        _Complex JMTX_SCALAR_T *const values = (_Complex JMTX_SCALAR_T *)in->main_diagonal;
        const JMTX_INDEX_T len = JMTX_NAME_TYPED(matrix_cds_entries_in_dia(in, 0);
        for (JMTX_INDEX_T i = 0; i < len; ++i)
        {
            values[i] = in->main_diagonal[i] * _Complex_I;
        }
    }

    for (JMTX_FAST_INT_T i = 0; i < in->sub_diagonals.count; ++i)
    {
        JMTX_SCALAR_T *ptr = in->sub_diagonals.diagonals[i];
        _Complex JMTX_SCALAR_T *const values = (_Complex JMTX_SCALAR_T *)ptr;
        const JMTX_INDEX_T len = JMTX_NAME_TYPED(matrix_cds_entries_in_dia(in, -((int32_t)i));
        for (JMTX_INDEX_T j = 0; j < len; ++j)
        {
            values[j] = ptr[j] * _Complex_I;
        }
    }

    for (JMTX_FAST_INT_T i = 0; i < in->super_diagonals.count; ++i)
    {
        JMTX_SCALAR_T *ptr = in->super_diagonals.diagonals[i];
        _Complex JMTX_SCALAR_T *const values = (_Complex JMTX_SCALAR_T *)ptr;
        const JMTX_INDEX_T len = JMTX_NAME_TYPED(matrix_cds_entries_in_dia(in, +((int32_t)i));
        for (JMTX_INDEX_T j = 0; j < len; ++j)
        {
            values[j] = (ptr[j]) * _Complex_I;
        }
    }

    in->base.type = JMTXC_TYPE_CDS;

    return (jmtxz_matrix_cds *)in;
}

/**
 * Creates a new complex CDS matrix with JMTX_SCALAR_T precision from a CDS matrix with JMTX_SCALAR_T precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in_real matrix to use as the real component of the new matrix (if NULL real part is zero)
 * @param in_imag matrix to use as the real component of the new matrix (if NULL imaginary part is zero)
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxzs_matrix_cds_from_double(jmtxz_matrix_cds **p_mtx, const JMTX_NAME_TYPED(matrix_cds) *in_real,
                                          const JMTX_NAME_TYPED(matrix_cds) *in_imag,
                                          const jmtx_allocator_callbacks *allocator_callbacks)
{
    if (!p_mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!in_real && !in_imag)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (allocator_callbacks &&
        (!allocator_callbacks->free || !allocator_callbacks->alloc || !allocator_callbacks->realloc))
    {
        return JMTX_RESULT_BAD_PARAM;
    }
    if ((in_real && in_real->base.type != JMTXD_TYPE_CDS) || (in_imag && in_imag->base.type != JMTXD_TYPE_CDS))
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if ((in_real && in_imag) && (in_real->base.rows != in_imag->base.rows || in_real->base.cols != in_imag->base.cols))
    {
        return JMTX_RESULT_BAD_MATRIX;
    }

    return jmtxz_matrix_cds_from_double(p_mtx, in_real, in_imag, allocator_callbacks);
}

/***********************************************************************************************************************
 *                                                                                                                     *
 *                                 COMPLEX FLOAT <->  COMPLEX DOUBLE                                                   *
 *                                                                                                                     *
 **********************************************************************************************************************/
/**
 * Creates a new complex CDS matrix with single precision from a complex CDS matrix with JMTX_SCALAR_T precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxc_matrix_cds_from_cdouble(jmtxc_matrix_cds **p_mtx, const jmtxz_matrix_cds *in,
                                          const jmtx_allocator_callbacks *allocator_callbacks)
{
    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }
    jmtxc_matrix_cds *mtx;
    jmtx_result res = jmtxc_matrix_cds_new(&mtx, in->base.rows, in->base.cols, 0, (int32_t[]){0}, allocator_callbacks);
    if (res != JMTX_RESULT_SUCCESS)
    {
        return res;
    }

    for (JMTX_FAST_INT_T i = 0; i < in->sub_diagonals.count; ++i)
    {
        const int32_t dia_idx = -(int32_t)in->sub_diagonals.indices[i];
        const _Complex JMTX_SCALAR_T *const dia_ptr = in->sub_diagonals.diagonals[i];
        JMTX_INDEX_T len;
        _Complex float *const ptr = jmtxc_matrix_cds_allocate_diagonal(mtx, dia_idx, &len);
        for (JMTX_FAST_INT_T j = 0; j < len; ++j)
        {
            ptr[j] = (_Complex float)dia_ptr[j];
        }
    }
    if (in->main_diagonal)
    {
        JMTX_INDEX_T len;
        _Complex float *const ptr = jmtxc_matrix_cds_allocate_diagonal(mtx, 0, &len);
        for (JMTX_FAST_INT_T j = 0; j < len; ++j)
        {
            ptr[j] = (_Complex float)in->main_diagonal[j];
        }
    }
    for (JMTX_FAST_INT_T i = 0; i < in->super_diagonals.count; ++i)
    {
        const int32_t dia_idx = (int32_t)in->super_diagonals.indices[i];
        const _Complex JMTX_SCALAR_T *const dia_ptr = in->super_diagonals.diagonals[i];
        JMTX_INDEX_T len;
        _Complex float *const ptr = jmtxc_matrix_cds_allocate_diagonal(mtx, dia_idx, &len);
        for (JMTX_FAST_INT_T j = 0; j < len; ++j)
        {
            ptr[j] = (_Complex float)dia_ptr[j];
        }
    }

    mtx->base.cols = in->base.cols;
    mtx->base.type = JMTXC_TYPE_CDS;
    mtx->base.rows = in->base.rows;
    mtx->base.allocator_callbacks = *allocator_callbacks;

    *p_mtx = mtx;

    return JMTX_RESULT_SUCCESS;
}
/**
 * Creates a new complex CDS matrix with single precision from a complex CDS matrix with JMTX_SCALAR_T precision. Requires no
 * memory allocation by reusing the memory of the initial matrix. Can not fail if the input matrix is valid.
 * @param in matrix which to convert (will be invalid if function succeeds)
 * @return converted matrix
 */
jmtxc_matrix_cds *jmtxc_matrix_cds_from_cdouble_inplace(jmtxz_matrix_cds *in)
{
    //  Convert diagonals
    if (in->main_diagonal)
    {
        _Complex float *const values = (_Complex float *)in->main_diagonal;
        const JMTX_INDEX_T len = jmtxz_matrix_cds_entries_in_dia(in, 0);
        for (JMTX_INDEX_T i = 0; i < len; ++i)
        {
            values[i] = (_Complex float)in->main_diagonal[i];
        }
    }

    for (JMTX_FAST_INT_T i = 0; i < in->sub_diagonals.count; ++i)
    {
        _Complex JMTX_SCALAR_T *ptr = in->sub_diagonals.diagonals[i];
        _Complex float *const values = (_Complex float *)ptr;
        const JMTX_INDEX_T len = jmtxz_matrix_cds_entries_in_dia(in, -((int32_t)i));
        for (JMTX_INDEX_T j = 0; j < len; ++j)
        {
            values[j] = (_Complex float)ptr[j];
        }
    }

    for (JMTX_FAST_INT_T i = 0; i < in->super_diagonals.count; ++i)
    {
        _Complex JMTX_SCALAR_T *ptr = in->super_diagonals.diagonals[i];
        _Complex float *const values = (_Complex float *)ptr;
        const JMTX_INDEX_T len = jmtxz_matrix_cds_entries_in_dia(in, +((int32_t)i));
        for (JMTX_INDEX_T j = 0; j < len; ++j)
        {
            values[j] = (_Complex float)ptr[j];
        }
    }

    //  Shrink the diagonals
    if (in->main_diagonal)
    {
        const JMTX_INDEX_T len = jmtxz_matrix_cds_entries_in_dia(in, 0);
        _Complex float *const new_ptr = in->base.allocator_callbacks.realloc(in->base.allocator_callbacks.state,
                                                                             in->main_diagonal, sizeof(*new_ptr) * len);
        if (new_ptr)
        {
            in->main_diagonal = (_Complex JMTX_SCALAR_T *)new_ptr;
        }
    }

    for (JMTX_FAST_INT_T i = 0; i < in->sub_diagonals.count; ++i)
    {
        const JMTX_INDEX_T len = jmtxz_matrix_cds_entries_in_dia(in, -((int32_t)i));
        _Complex float *const new_ptr = in->base.allocator_callbacks.realloc(
            in->base.allocator_callbacks.state, in->sub_diagonals.diagonals[i], sizeof(*new_ptr) * len);
        if (new_ptr)
        {
            in->sub_diagonals.diagonals[i] = (_Complex JMTX_SCALAR_T *)new_ptr;
        }
    }

    for (JMTX_FAST_INT_T i = 0; i < in->super_diagonals.count; ++i)
    {
        const JMTX_INDEX_T len = jmtxz_matrix_cds_entries_in_dia(in, +((int32_t)i));
        _Complex float *const new_ptr = in->base.allocator_callbacks.realloc(
            in->base.allocator_callbacks.state, in->super_diagonals.diagonals[i], sizeof(*new_ptr) * len);
        if (new_ptr)
        {
            in->super_diagonals.diagonals[i] = (_Complex JMTX_SCALAR_T *)new_ptr;
        }
    }

    in->base.type = JMTXC_TYPE_CDS;

    return (jmtxc_matrix_cds *)in;
}

/**
 * Creates a new complex CDS matrix with single precision from a complex CDS matrix with JMTX_SCALAR_T precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxcs_matrix_cds_from_cdouble(jmtxc_matrix_cds **p_mtx, const jmtxz_matrix_cds *in,
                                           const jmtx_allocator_callbacks *allocator_callbacks)
{
    if (!p_mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!in)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (allocator_callbacks &&
        (!allocator_callbacks->free || !allocator_callbacks->alloc || !allocator_callbacks->realloc))
    {
        return JMTX_RESULT_BAD_PARAM;
    }
    if (in->base.type != JMTXZ_TYPE_CDS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }

    return jmtxc_matrix_cds_from_cdouble(p_mtx, in, allocator_callbacks);
}

/**
 * Creates a new complex CDS matrix with JMTX_SCALAR_T precision from a complex CDS matrix with single precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxz_matrix_cds_from_cfloat(jmtxz_matrix_cds **p_mtx, const jmtxc_matrix_cds *in,
                                         const jmtx_allocator_callbacks *allocator_callbacks)
{
    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }
    jmtxz_matrix_cds *mtx;
    jmtx_result res = jmtxz_matrix_cds_new(&mtx, in->base.rows, in->base.cols, 0, (int32_t[]){0}, allocator_callbacks);
    if (res != JMTX_RESULT_SUCCESS)
    {
        return res;
    }

    for (JMTX_FAST_INT_T i = 0; i < in->sub_diagonals.count; ++i)
    {
        const int32_t dia_idx = -(int32_t)in->sub_diagonals.indices[i];
        const _Complex float *const dia_ptr = in->sub_diagonals.diagonals[i];
        JMTX_INDEX_T len;
        _Complex JMTX_SCALAR_T *const ptr = jmtxz_matrix_cds_allocate_diagonal(mtx, dia_idx, &len);
        for (JMTX_FAST_INT_T j = 0; j < len; ++j)
        {
            ptr[j] = (_Complex double)dia_ptr[j];
        }
    }
    if (in->main_diagonal)
    {
        JMTX_INDEX_T len;
        _Complex JMTX_SCALAR_T *const ptr = jmtxz_matrix_cds_allocate_diagonal(mtx, 0, &len);
        for (JMTX_FAST_INT_T j = 0; j < len; ++j)
        {
            ptr[j] = (_Complex double)in->main_diagonal[j];
        }
    }
    for (JMTX_FAST_INT_T i = 0; i < in->super_diagonals.count; ++i)
    {
        const int32_t dia_idx = (int32_t)in->super_diagonals.indices[i];
        const _Complex float *const dia_ptr = in->super_diagonals.diagonals[i];
        JMTX_INDEX_T len;
        _Complex JMTX_SCALAR_T *const ptr = jmtxz_matrix_cds_allocate_diagonal(mtx, dia_idx, &len);
        for (JMTX_FAST_INT_T j = 0; j < len; ++j)
        {
            ptr[j] = (_Complex double)dia_ptr[j];
        }
    }

    mtx->base.cols = in->base.cols;
    mtx->base.type = JMTXZ_TYPE_CDS;
    mtx->base.rows = in->base.rows;
    mtx->base.allocator_callbacks = *allocator_callbacks;

    *p_mtx = mtx;

    return JMTX_RESULT_SUCCESS;
}

/**
 * Creates a new complex CDS matrix with JMTX_SCALAR_T precision from a complex CDS matrix with single precision. Only one
 * memory reallocation may be needed. Can not fail if the input matrix is valid.
 * @param in matrix which to convert
 * @return converted matrix, or NULL in case of allocation failure
 */
jmtxz_matrix_cds *jmtxz_matrix_cds_from_cfloat_inplace(jmtxc_matrix_cds *in)
{
    //  Expand the diagonals
    if (in->main_diagonal)
    {
        const JMTX_INDEX_T len = jmtxc_matrix_cds_entries_in_dia(in, 0);
        _Complex JMTX_SCALAR_T *const new_ptr = in->base.allocator_callbacks.realloc(
            in->base.allocator_callbacks.state, in->main_diagonal, sizeof(*new_ptr) * len);
        if (!new_ptr)
        {
            return NULL;
        }
        in->main_diagonal = (_Complex float *)new_ptr;
    }

    for (JMTX_FAST_INT_T i = 0; i < in->sub_diagonals.count; ++i)
    {
        const JMTX_INDEX_T len = jmtxc_matrix_cds_entries_in_dia(in, -((int32_t)i));
        _Complex JMTX_SCALAR_T *const new_ptr = in->base.allocator_callbacks.realloc(
            in->base.allocator_callbacks.state, in->sub_diagonals.diagonals[i], sizeof(*new_ptr) * len);
        if (!new_ptr)
        {
            return NULL;
        }
        in->sub_diagonals.diagonals[i] = (_Complex float *)new_ptr;
    }

    for (JMTX_FAST_INT_T i = 0; i < in->super_diagonals.count; ++i)
    {
        const JMTX_INDEX_T len = jmtxc_matrix_cds_entries_in_dia(in, +((int32_t)i));
        _Complex JMTX_SCALAR_T *const new_ptr = in->base.allocator_callbacks.realloc(
            in->base.allocator_callbacks.state, in->super_diagonals.diagonals[i], sizeof(*new_ptr) * len);
        if (!new_ptr)
        {
            return NULL;
        }
        in->super_diagonals.diagonals[i] = (_Complex float *)new_ptr;
    }

    //  Convert diagonals
    if (in->main_diagonal)
    {
        _Complex JMTX_SCALAR_T *const values = (_Complex JMTX_SCALAR_T *)in->main_diagonal;
        const JMTX_INDEX_T len = jmtxc_matrix_cds_entries_in_dia(in, 0);
        for (JMTX_INDEX_T i = 0; i < len; ++i)
        {
            values[len - 1 - i] = (_Complex double)in->main_diagonal[len - 1 - i];
        }
    }

    for (JMTX_FAST_INT_T i = 0; i < in->sub_diagonals.count; ++i)
    {
        _Complex float *ptr = in->sub_diagonals.diagonals[i];
        _Complex JMTX_SCALAR_T *const values = (_Complex JMTX_SCALAR_T *)ptr;
        const JMTX_INDEX_T len = jmtxc_matrix_cds_entries_in_dia(in, -((int32_t)i));
        for (JMTX_INDEX_T j = 0; j < len; ++j)
        {
            values[len - 1 - j] = (_Complex double)ptr[len - 1 - j];
        }
    }

    for (JMTX_FAST_INT_T i = 0; i < in->super_diagonals.count; ++i)
    {
        _Complex float *ptr = in->super_diagonals.diagonals[i];
        _Complex JMTX_SCALAR_T *const values = (_Complex JMTX_SCALAR_T *)ptr;
        const JMTX_INDEX_T len = jmtxc_matrix_cds_entries_in_dia(in, +((int32_t)i));
        for (JMTX_INDEX_T j = 0; j < len; ++j)
        {
            values[len - 1 - j] = (_Complex double)ptr[len - 1 - j];
        }
    }

    in->base.type = JMTXZ_TYPE_CDS;

    return (jmtxz_matrix_cds *)in;
}

/**
 * Creates a new complex CDS matrix with JMTX_SCALAR_T precision from a complex CDS matrix with single precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxzs_matrix_cds_from_cfloat(jmtxz_matrix_cds **p_mtx, const jmtxc_matrix_cds *in,
                                          const jmtx_allocator_callbacks *allocator_callbacks)
{
    if (!p_mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!in)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (allocator_callbacks &&
        (!allocator_callbacks->free || !allocator_callbacks->alloc || !allocator_callbacks->realloc))
    {
        return JMTX_RESULT_BAD_PARAM;
    }
    if (in->base.type != JMTXC_TYPE_CDS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }

    return jmtxz_matrix_cds_from_cfloat(p_mtx, in, allocator_callbacks);
}
#endif //!_MSC_BUILD
