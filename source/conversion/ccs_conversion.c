//
// Created by jan on 1.12.2023.
//

#include <assert.h>
#include <complex.h>

#include "../double/matrices/sparse_column_compressed.h"
#include "../float/matrices/sparse_column_compressed.h"
#include "../matrix_base.h"
#ifndef _MSC_BUILD
#    include "../cdouble/matrices/sparse_column_compressed.h"
#    include "../cfloat/matrices/sparse_column_compressed.h"
#endif
#include "ccs_conversion.h"

/***********************************************************************************************************************
 *                                                                                                                     *
 *                                          FLOAT <-> DOUBLE                                                           *
 *                                                                                                                     *
 **********************************************************************************************************************/

/**
 * Creates a new CCS matrix with single precision from a CCS matrix with JMTX_SCALAR_T precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxf_matrix_ccs_from_double(jmtxf_matrix_ccs **p_mtx, const JMTX_NAME_TYPED(matrix_ccs) * in,
                                         const jmtx_allocator_callbacks *allocator_callbacks)
{
    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }

    JMTX_INDEX_T *offsets = NULL;
    JMTX_INDEX_T *indices = NULL;

    jmtxf_matrix_ccs *mtx = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*mtx));
    if (!mtx)
    {
        return JMTX_RESULT_BAD_ALLOC;
    }

    float *values = allocator_callbacks->alloc(allocator_callbacks->state, (in->n_entries) * sizeof(*values));
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
        values[i] = (float)in->values[i];
    }

    mtx->base.cols = in->base.cols;
    mtx->base.type = JMTX_TYPE_CCS;
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
jmtxf_matrix_ccs *jmtxf_matrix_ccs_from_double_inplace(JMTX_NAME_TYPED(matrix_ccs) * in)
{
    float *const values = (float *)in->values;

    for (JMTX_FAST_INT_T i = 0; i < in->n_entries; ++i)
    {
        values[i] = (float)in->values[i];
    }

    in->capacity = (in->capacity * sizeof(*in->values)) / sizeof(*values);
    in->base.type = JMTX_TYPE_CCS;

    return (jmtxf_matrix_ccs *)in;
}

/**
 * Creates a new CCS matrix with single precision from a CCS matrix with JMTX_SCALAR_T precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxs_matrix_ccs_from_double(jmtxf_matrix_ccs **p_mtx, const JMTX_NAME_TYPED(matrix_ccs) * in,
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
    if (in->base.type != JMTXD_TYPE_CCS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }

    return jmtxf_matrix_ccs_from_double(p_mtx, in, allocator_callbacks);
}

/**
 * Creates a new CCS matrix with JMTX_SCALAR_T precision from a CCS matrix with single precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result JMTX_NAME_TYPED(matrix_ccs_from_float(JMTX_NAME_TYPED(matrix_ccs) **p_mtx, const jmtxf_matrix_ccs *in,
                                        const jmtx_allocator_callbacks *allocator_callbacks)
{
    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }

    JMTX_INDEX_T *offsets = NULL;
    JMTX_INDEX_T *indices = NULL;

    JMTX_NAME_TYPED(matrix_ccs) *mtx = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*mtx));
    if (!mtx)
    {
        return JMTX_RESULT_BAD_ALLOC;
    }

    JMTX_SCALAR_T *values = allocator_callbacks->alloc(allocator_callbacks->state, (in->n_entries) * sizeof(*values));
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
        values[i] = (JMTX_SCALAR_T)in->values[i];
    }

    mtx->base.cols = in->base.cols;
    mtx->base.type = JMTXD_TYPE_CCS;
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
 * Creates a new CCS matrix with JMTX_SCALAR_T precision from a CCS matrix with single precision. Only one memory reallocation
 * may be needed. Can not fail if the input matrix is valid.
 * @param in matrix which to convert
 * @return converted matrix, or NULL in case of allocation failure
 */
JMTX_NAME_TYPED(matrix_ccs) *JMTX_NAME_TYPED(matrix_ccs_from_float_inplace(jmtxf_matrix_ccs *in)
{
    JMTX_SCALAR_T *vals = NULL;
    JMTX_INDEX_T capacity_for_doubles = in->capacity * sizeof(*in->values) / sizeof(*vals);
    if (capacity_for_doubles < in->n_entries)
    {
        capacity_for_doubles = in->n_entries;
        vals = in->base.allocator_callbacks.realloc(in->base.allocator_callbacks.state, in->values,
                                                    sizeof(*vals) * capacity_for_doubles);
        if (!vals)
        {
            return NULL;
        }
        in->values = (float *)vals;
    }

    vals = (JMTX_SCALAR_T *)in->values;
    in->capacity = capacity_for_doubles;

    for (JMTX_INDEX_T i = 0; i < in->n_entries; ++i)
    {
        vals[in->n_entries - 1 - i] = (JMTX_SCALAR_T)in->values[in->n_entries - 1 - i];
    }
    in->base.type = JMTXD_TYPE_CCS;

    return (JMTX_NAME_TYPED(matrix_ccs) *)in;
}

/**
 * Creates a new CCS matrix with JMTX_SCALAR_T precision from a CCS matrix with single precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result .,.,.,.,..matrix_ccs_from_float(JMTX_NAME_TYPED(matrix_ccs) **p_mtx, const jmtxf_matrix_ccs *in,
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
    if (in->base.type != JMTX_TYPE_CCS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }

    return JMTX_NAME_TYPED(matrix_ccs_from_float(p_mtx, in, allocator_callbacks);
}

#ifndef _MSC_BUILD
/***********************************************************************************************************************
 *                                                                                                                     *
 *                                          FLOAT <-> COMPLEX FLOAT                                                    *
 *                                                                                                                     *
 **********************************************************************************************************************/

/**
 * Creates a new CCS matrix with single precision from the real part of a complex CCS matrix with single precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxf_matrix_ccs_from_cfloat_real(jmtxf_matrix_ccs **p_mtx, const jmtxc_matrix_ccs *in,
                                              const jmtx_allocator_callbacks *allocator_callbacks)
{
    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }

    JMTX_INDEX_T *offsets = NULL;
    JMTX_INDEX_T *indices = NULL;

    jmtxf_matrix_ccs *mtx = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*mtx));
    if (!mtx)
    {
        return JMTX_RESULT_BAD_ALLOC;
    }

    float *values = allocator_callbacks->alloc(allocator_callbacks->state, (in->n_entries) * sizeof(*values));
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
        values[i] = crealf(in->values[i]);
    }

    mtx->base.cols = in->base.cols;
    mtx->base.type = JMTX_TYPE_CCS;
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
 * Creates a new CCS matrix with single precision from the imaginary part of a complex CCS matrix with single precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxf_matrix_ccs_from_cfloat_imag(jmtxf_matrix_ccs **p_mtx, const jmtxc_matrix_ccs *in,
                                              const jmtx_allocator_callbacks *allocator_callbacks)
{
    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }

    JMTX_INDEX_T *offsets = NULL;
    JMTX_INDEX_T *indices = NULL;

    jmtxf_matrix_ccs *mtx = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*mtx));
    if (!mtx)
    {
        return JMTX_RESULT_BAD_ALLOC;
    }

    float *values = allocator_callbacks->alloc(allocator_callbacks->state, (in->n_entries) * sizeof(*values));
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
        values[i] = cimagf(in->values[i]);
    }

    mtx->base.cols = in->base.cols;
    mtx->base.type = JMTX_TYPE_CCS;
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
 * Creates a new CCS matrix with single precision from the real part of a complex CCS matrix with single precision.
 * Requires no memory allocation by reusing the memory of the initial matrix. Can not fail if the input matrix is valid.
 * @param in matrix which to convert (will be invalid if function succeeds)
 * @return converted matrix
 */
jmtxf_matrix_ccs *jmtxf_matrix_ccs_from_cfloat_real_inplace(jmtxc_matrix_ccs *in)
{
    float *const values = (float *)in->values;

    for (JMTX_FAST_INT_T i = 0; i < in->n_entries; ++i)
    {
        values[i] = crealf(in->values[i]);
    }

    in->capacity = (in->capacity * sizeof(*in->values)) / sizeof(*values);
    in->base.type = JMTX_TYPE_CCS;

    return (jmtxf_matrix_ccs *)in;
}

/**
 * Creates a new CCS matrix with single precision from the imaginary part of a complex CCS matrix with single precision.
 * Requires no memory allocation by reusing the memory of the initial matrix. Can not fail if the input matrix is valid.
 * @param in matrix which to convert (will be invalid if function succeeds)
 * @return converted matrix
 */
jmtxf_matrix_ccs *jmtxf_matrix_ccs_from_cfloat_imag_inplace(jmtxc_matrix_ccs *in)
{
    float *const values = (float *)in->values;

    for (JMTX_FAST_INT_T i = 0; i < in->n_entries; ++i)
    {
        values[i] = cimagf(in->values[i]);
    }

    in->capacity = (in->capacity * sizeof(*in->values)) / sizeof(*values);
    in->base.type = JMTX_TYPE_CCS;

    return (jmtxf_matrix_ccs *)in;
}

/**
 * Creates a new CCS matrix with single precision from the real part of a complex CCS matrix with single precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxs_matrix_ccs_from_cfloat_real(jmtxf_matrix_ccs **p_mtx, const jmtxc_matrix_ccs *in,
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
    if (in->base.type != JMTXC_TYPE_CCS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }

    return jmtxf_matrix_ccs_from_cfloat_real(p_mtx, in, allocator_callbacks);
}

/**
 * Creates a new CCS matrix with single precision from the imaginary part of a complex CCS matrix with single precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxs_matrix_ccs_from_cfloat_imag(jmtxf_matrix_ccs **p_mtx, const jmtxc_matrix_ccs *in,
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
    if (in->base.type != JMTXC_TYPE_CCS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }

    return jmtxf_matrix_ccs_from_cfloat_imag(p_mtx, in, allocator_callbacks);
}

/**
 * Creates a new complex CCS matrix with single precision from a CCS matrix with single precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in_real matrix to use as the real component of the new matrix (if NULL real part is zero)
 * @param in_imag matrix to use as the real component of the new matrix (if NULL imaginary part is zero)
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxc_matrix_ccs_from_float(jmtxc_matrix_ccs **p_mtx, const jmtxf_matrix_ccs *in_real,
                                        const jmtxf_matrix_ccs *in_imag,
                                        const jmtx_allocator_callbacks *allocator_callbacks)
{
    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }

    JMTX_INDEX_T *offsets = NULL;
    JMTX_INDEX_T *indices = NULL;
    JMTX_INDEX_T r, c;

    assert(in_real || in_imag);

    if (in_real)
    {
        r = in_real->base.rows;
        c = in_real->base.cols;
    }
    else // if (in_imag)
    {
        r = in_imag->base.rows;
        c = in_imag->base.cols;
    }
    jmtxc_matrix_ccs *mtx = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*mtx));
    if (!mtx)
    {
        return JMTX_RESULT_BAD_ALLOC;
    }

    if (in_real && in_imag)
    {
        const JMTX_INDEX_T needed_capacity = in_real->n_entries + in_imag->n_entries;

        _Complex float *values =
            allocator_callbacks->alloc(allocator_callbacks->state, (needed_capacity) * sizeof(*values));
        if (!values)
        {
            allocator_callbacks->free(allocator_callbacks->state, mtx);
            return JMTX_RESULT_BAD_ALLOC;
        }

        indices = allocator_callbacks->alloc(allocator_callbacks->state, (needed_capacity) * sizeof(*indices));
        if (!indices)
        {
            allocator_callbacks->free(allocator_callbacks->state, indices);
            allocator_callbacks->free(allocator_callbacks->state, mtx);
            return JMTX_RESULT_BAD_ALLOC;
        }

        offsets = allocator_callbacks->alloc(allocator_callbacks->state, (c) * sizeof(*offsets));
        if (!offsets)
        {
            allocator_callbacks->free(allocator_callbacks->state, offsets);
            allocator_callbacks->free(allocator_callbacks->state, indices);
            allocator_callbacks->free(allocator_callbacks->state, mtx);
            return JMTX_RESULT_BAD_ALLOC;
        }
        JMTX_INDEX_T i, j;
        for (i = 0, j = 0; i < c; ++i)
        {
            float *v_r, *v_i;
            JMTX_INDEX_T *i_r, *i_i;
            const JMTX_INDEX_T c_r = jmtxf_matrix_ccs_get_col(in_real, i, &i_r, &v_r);
            const JMTX_INDEX_T c_i = jmtxf_matrix_ccs_get_col(in_imag, i, &i_i, &v_i);
            JMTX_INDEX_T k_r, k_i;
            for (k_r = 0, k_i = 0; k_r < c_r && k_i < c_i; ++j)
            {
                if (i_r[k_r] < i_i[k_i])
                {
                    values[j] = v_r[k_r];
                    indices[j] = i_r[k_r];
                    k_r += 1;
                }
                else if (i_r[k_r] > i_i[k_i])
                {
                    values[j] = v_i[k_i] * _Complex_I;
                    indices[j] = i_i[k_r];
                    k_i += 1;
                }
                else // if (i_r[k_r] == i_i[k_i])
                {
                    values[j] = v_i[k_i] * _Complex_I + v_r[k_r];
                    indices[j] = i_i[k_r];
                    k_r += 1;
                    k_i += 1;
                }
            }

            if (k_i < c_i)
            {
                assert(k_r == c_r);
                while (k_i < c_i)
                {
                    values[j] = v_i[k_i] * _Complex_I;
                    indices[j] = i_i[k_i];
                    j += 1;
                    k_i += 1;
                }
            }
            else
            {
                assert(k_i == c_i);
                while (k_r < c_r)
                {
                    values[j] = v_r[k_r];
                    indices[j] = i_r[k_r];
                    j += 1;
                    k_r += 1;
                }
            }

            offsets[i] = j;
        }

        mtx->indices = indices;
        mtx->values = values;
        mtx->capacity = needed_capacity;
        mtx->n_entries = j;
        mtx->end_of_column_offsets = offsets;
    }
    else if (in_imag)
    {
        _Complex float *values =
            allocator_callbacks->alloc(allocator_callbacks->state, (in_imag->n_entries) * sizeof(*values));
        if (!values)
        {
            allocator_callbacks->free(allocator_callbacks->state, mtx);
            return JMTX_RESULT_BAD_ALLOC;
        }

        indices = allocator_callbacks->alloc(allocator_callbacks->state, (in_imag->n_entries) * sizeof(*indices));
        if (!indices)
        {
            allocator_callbacks->free(allocator_callbacks->state, indices);
            allocator_callbacks->free(allocator_callbacks->state, mtx);
            return JMTX_RESULT_BAD_ALLOC;
        }

        offsets = allocator_callbacks->alloc(allocator_callbacks->state, (in_imag->base.cols) * sizeof(*offsets));
        if (!offsets)
        {
            allocator_callbacks->free(allocator_callbacks->state, offsets);
            allocator_callbacks->free(allocator_callbacks->state, indices);
            allocator_callbacks->free(allocator_callbacks->state, mtx);
            return JMTX_RESULT_BAD_ALLOC;
        }
        memcpy(offsets, in_imag->end_of_column_offsets, (in_imag->base.cols) * sizeof(*offsets));
        memcpy(indices, in_imag->indices, (in_imag->n_entries) * sizeof(*indices));

        for (JMTX_FAST_INT_T i = 0; i < in_imag->n_entries; ++i)
        {
            values[i] = _Complex_I * in_imag->values[i];
        }

        mtx->indices = indices;
        mtx->values = values;
        mtx->capacity = in_imag->n_entries;
        mtx->n_entries = in_imag->n_entries;
        mtx->end_of_column_offsets = offsets;
    }
    else // if (in_real)
    {
        _Complex float *values =
            allocator_callbacks->alloc(allocator_callbacks->state, (in_real->n_entries) * sizeof(*values));
        if (!values)
        {
            allocator_callbacks->free(allocator_callbacks->state, mtx);
            return JMTX_RESULT_BAD_ALLOC;
        }

        indices = allocator_callbacks->alloc(allocator_callbacks->state, (in_real->n_entries) * sizeof(*indices));
        if (!indices)
        {
            allocator_callbacks->free(allocator_callbacks->state, indices);
            allocator_callbacks->free(allocator_callbacks->state, mtx);
            return JMTX_RESULT_BAD_ALLOC;
        }

        offsets = allocator_callbacks->alloc(allocator_callbacks->state, (in_real->base.cols) * sizeof(*offsets));
        if (!offsets)
        {
            allocator_callbacks->free(allocator_callbacks->state, offsets);
            allocator_callbacks->free(allocator_callbacks->state, indices);
            allocator_callbacks->free(allocator_callbacks->state, mtx);
            return JMTX_RESULT_BAD_ALLOC;
        }
        memcpy(offsets, in_real->end_of_column_offsets, (in_real->base.cols) * sizeof(*offsets));
        memcpy(indices, in_real->indices, (in_real->n_entries) * sizeof(*indices));

        for (JMTX_FAST_INT_T i = 0; i < in_real->n_entries; ++i)
        {
            values[i] = _Complex_I * in_real->values[i];
        }

        mtx->indices = indices;
        mtx->values = values;
        mtx->capacity = in_real->n_entries;
        mtx->n_entries = in_real->n_entries;
        mtx->end_of_column_offsets = offsets;
    }

    mtx->base.cols = c;
    mtx->base.type = JMTXC_TYPE_CCS;
    mtx->base.rows = r;
    mtx->base.allocator_callbacks = *allocator_callbacks;
    *p_mtx = mtx;

    return JMTX_RESULT_SUCCESS;
}

/**
 * Creates a new complex CCS matrix with single precision from a CCS matrix with single precision as its real part.
 * Only one memory reallocation.
 * may be needed. Can not fail if the input matrix is valid.
 * @param in matrix which to convert
 * @return converted matrix, or NULL in case of allocation failure
 */
jmtxc_matrix_ccs *jmtxc_matrix_ccs_from_float_real_inplace(jmtxf_matrix_ccs *in)
{
    _Complex float *vals = NULL;
    JMTX_INDEX_T capacity_for_complex = in->capacity * sizeof(*in->values) / sizeof(*vals);
    if (capacity_for_complex < in->n_entries)
    {
        capacity_for_complex = in->n_entries;
        vals = in->base.allocator_callbacks.realloc(in->base.allocator_callbacks.state, in->values,
                                                    sizeof(*vals) * capacity_for_complex);
        if (!vals)
        {
            return NULL;
        }
        in->values = (float *)vals;
    }

    vals = (_Complex float *)in->values;
    in->capacity = capacity_for_complex;

    for (JMTX_INDEX_T i = 0; i < in->n_entries; ++i)
    {
        vals[in->n_entries - 1 - i] = (_Complex float)in->values[in->n_entries - 1 - i];
    }
    in->base.type = JMTXC_TYPE_CCS;

    return (jmtxc_matrix_ccs *)in;
}

/**
 * Creates a new CCS matrix with single precision from a CCS matrix with single precision as its imaginary part.
 * Only one memory reallocation.
 * may be needed. Can not fail if the input matrix is valid.
 * @param in matrix which to convert
 * @return converted matrix, or NULL in case of allocation failure
 */
jmtxc_matrix_ccs *jmtxc_matrix_ccs_from_float_imag_inplace(jmtxf_matrix_ccs *in)
{
    _Complex float *vals = NULL;
    JMTX_INDEX_T capacity_for_complex = in->capacity * sizeof(*in->values) / sizeof(*vals);
    if (capacity_for_complex < in->n_entries)
    {
        capacity_for_complex = in->n_entries;
        vals = in->base.allocator_callbacks.realloc(in->base.allocator_callbacks.state, in->values,
                                                    sizeof(*vals) * capacity_for_complex);
        if (!vals)
        {
            return NULL;
        }
        in->values = (float *)vals;
    }

    vals = (_Complex float *)in->values;
    in->capacity = capacity_for_complex;

    for (JMTX_INDEX_T i = 0; i < in->n_entries; ++i)
    {
        vals[in->n_entries - 1 - i] = (_Complex float)in->values[in->n_entries - 1 - i] * _Complex_I;
    }
    in->base.type = JMTXC_TYPE_CCS;

    return (jmtxc_matrix_ccs *)in;
}

/**
 * Creates a new complex CCS matrix with single precision from a CCS matrix with single precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in_real matrix to use as the real component of the new matrix (if NULL real part is zero)
 * @param in_imag matrix to use as the real component of the new matrix (if NULL imaginary part is zero)
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxcs_matrix_ccs_from_float(jmtxc_matrix_ccs **p_mtx, const jmtxf_matrix_ccs *in_real,
                                         const jmtxf_matrix_ccs *in_imag,
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
    if ((in_real && in_real->base.type != JMTX_TYPE_CCS) || (in_imag && in_imag->base.type != JMTX_TYPE_CCS))
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if ((in_real && in_imag) && (in_real->base.rows != in_imag->base.rows || in_real->base.cols != in_imag->base.cols))
    {
        return JMTX_RESULT_BAD_MATRIX;
    }

    return jmtxc_matrix_ccs_from_float(p_mtx, in_real, in_imag, allocator_callbacks);
}

/***********************************************************************************************************************
 *                                                                                                                     *
 *                                          DOUBLE <-> COMPLEX DOUBLE                                                  *
 *                                                                                                                     *
 **********************************************************************************************************************/

/**
 * Creates a new CCS matrix with JMTX_SCALAR_T precision from the real part of a complex CCS matrix with JMTX_SCALAR_T precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result JMTX_NAME_TYPED(matrix_ccs_from_cdouble_real(JMTX_NAME_TYPED(matrix_ccs) **p_mtx, const jmtxz_matrix_ccs *in,
                                               const jmtx_allocator_callbacks *allocator_callbacks)
{
    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }

    JMTX_INDEX_T *offsets = NULL;
    JMTX_INDEX_T *indices = NULL;

    JMTX_NAME_TYPED(matrix_ccs) *mtx = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*mtx));
    if (!mtx)
    {
        return JMTX_RESULT_BAD_ALLOC;
    }

    JMTX_SCALAR_T *values = allocator_callbacks->alloc(allocator_callbacks->state, (in->n_entries) * sizeof(*values));
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
        values[i] = creal(in->values[i]);
    }

    mtx->base.cols = in->base.cols;
    mtx->base.type = JMTXD_TYPE_CCS;
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
 * Creates a new CCS matrix with JMTX_SCALAR_T precision from the imaginary part of a complex CCS matrix with JMTX_SCALAR_T precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result JMTX_NAME_TYPED(matrix_ccs_from_cdouble_imag(JMTX_NAME_TYPED(matrix_ccs) **p_mtx, const jmtxz_matrix_ccs *in,
                                               const jmtx_allocator_callbacks *allocator_callbacks)
{
    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }

    JMTX_INDEX_T *offsets = NULL;
    JMTX_INDEX_T *indices = NULL;

    JMTX_NAME_TYPED(matrix_ccs) *mtx = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*mtx));
    if (!mtx)
    {
        return JMTX_RESULT_BAD_ALLOC;
    }

    JMTX_SCALAR_T *values = allocator_callbacks->alloc(allocator_callbacks->state, (in->n_entries) * sizeof(*values));
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
        values[i] = cimag(in->values[i]);
    }

    mtx->base.cols = in->base.cols;
    mtx->base.type = JMTXD_TYPE_CCS;
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
 * Creates a new CCS matrix with JMTX_SCALAR_T precision from the real part of a complex CCS matrix with JMTX_SCALAR_T precision.
 * Requires no memory allocation by reusing the memory of the initial matrix. Can not fail if the input matrix is valid.
 * @param in matrix which to convert (will be invalid if function succeeds)
 * @return converted matrix
 */
JMTX_NAME_TYPED(matrix_ccs) *JMTX_NAME_TYPED(matrix_ccs_from_cdouble_real_inplace(jmtxz_matrix_ccs *in)
{
    JMTX_SCALAR_T *const values = (JMTX_SCALAR_T *)in->values;

    for (JMTX_FAST_INT_T i = 0; i < in->n_entries; ++i)
    {
        values[i] = creal(in->values[i]);
    }

    in->capacity = (in->capacity * sizeof(*in->values)) / sizeof(*values);
    in->base.type = JMTXD_TYPE_CCS;

    return (JMTX_NAME_TYPED(matrix_ccs) *)in;
}

/**
 * Creates a new CCS matrix with JMTX_SCALAR_T precision from the imaginary part of a complex CCS matrix with JMTX_SCALAR_T precision.
 * Requires no memory allocation by reusing the memory of the initial matrix. Can not fail if the input matrix is valid.
 * @param in matrix which to convert (will be invalid if function succeeds)
 * @return converted matrix
 */
JMTX_NAME_TYPED(matrix_ccs) *JMTX_NAME_TYPED(matrix_ccs_from_cdouble_imag_inplace(jmtxz_matrix_ccs *in)
{
    JMTX_SCALAR_T *const values = (JMTX_SCALAR_T *)in->values;

    for (JMTX_FAST_INT_T i = 0; i < in->n_entries; ++i)
    {
        values[i] = cimag(in->values[i]);
    }

    in->capacity = (in->capacity * sizeof(*in->values)) / sizeof(*values);
    in->base.type = JMTX_TYPE_CCS;

    return (JMTX_NAME_TYPED(matrix_ccs) *)in;
}

/**
 * Creates a new CCS matrix with JMTX_SCALAR_T precision from the real part of a complex CCS matrix with JMTX_SCALAR_T precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result .,.,.,.,..matrix_ccs_from_cdouble_real(JMTX_NAME_TYPED(matrix_ccs) **p_mtx, const jmtxz_matrix_ccs *in,
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
    if (in->base.type != JMTXZ_TYPE_CCS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }

    return JMTX_NAME_TYPED(matrix_ccs_from_cdouble_real(p_mtx, in, allocator_callbacks);
}

/**
 * Creates a new CCS matrix with JMTX_SCALAR_T precision from the imaginary part of a complex CCS matrix with JMTX_SCALAR_T precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result .,.,.,.,..matrix_ccs_from_cdouble_imag(JMTX_NAME_TYPED(matrix_ccs) **p_mtx, const jmtxz_matrix_ccs *in,
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
    if (in->base.type != JMTXZ_TYPE_CCS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }

    return JMTX_NAME_TYPED(matrix_ccs_from_cdouble_imag(p_mtx, in, allocator_callbacks);
}

/**
 * Creates a new complex CCS matrix with JMTX_SCALAR_T precision from a CCS matrix with JMTX_SCALAR_T precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in_real matrix to use as the real component of the new matrix (if NULL real part is zero)
 * @param in_imag matrix to use as the real component of the new matrix (if NULL imaginary part is zero)
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxz_matrix_ccs_from_double(jmtxz_matrix_ccs **p_mtx, const JMTX_NAME_TYPED(matrix_ccs) *in_real,
                                         const JMTX_NAME_TYPED(matrix_ccs) *in_imag,
                                         const jmtx_allocator_callbacks *allocator_callbacks)
{
    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }

    JMTX_INDEX_T *offsets = NULL;
    JMTX_INDEX_T *indices = NULL;
    JMTX_INDEX_T r, c;

    assert(in_real || in_imag);

    if (in_real)
    {
        r = in_real->base.rows;
        c = in_real->base.cols;
    }
    else // if (in_imag)
    {
        r = in_imag->base.rows;
        c = in_imag->base.cols;
    }
    jmtxz_matrix_ccs *mtx = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*mtx));
    if (!mtx)
    {
        return JMTX_RESULT_BAD_ALLOC;
    }

    if (in_real && in_imag)
    {
        const JMTX_INDEX_T needed_capacity = in_real->n_entries + in_imag->n_entries;

        _Complex JMTX_SCALAR_T *values =
            allocator_callbacks->alloc(allocator_callbacks->state, (needed_capacity) * sizeof(*values));
        if (!values)
        {
            allocator_callbacks->free(allocator_callbacks->state, mtx);
            return JMTX_RESULT_BAD_ALLOC;
        }

        indices = allocator_callbacks->alloc(allocator_callbacks->state, (needed_capacity) * sizeof(*indices));
        if (!indices)
        {
            allocator_callbacks->free(allocator_callbacks->state, indices);
            allocator_callbacks->free(allocator_callbacks->state, mtx);
            return JMTX_RESULT_BAD_ALLOC;
        }

        offsets = allocator_callbacks->alloc(allocator_callbacks->state, (c) * sizeof(*offsets));
        if (!offsets)
        {
            allocator_callbacks->free(allocator_callbacks->state, offsets);
            allocator_callbacks->free(allocator_callbacks->state, indices);
            allocator_callbacks->free(allocator_callbacks->state, mtx);
            return JMTX_RESULT_BAD_ALLOC;
        }
        JMTX_INDEX_T i, j;
        for (i = 0, j = 0; i < c; ++i)
        {
            JMTX_SCALAR_T *v_r, *v_i;
            JMTX_INDEX_T *i_r, *i_i;
            const JMTX_INDEX_T c_r = JMTX_NAME_TYPED(matrix_ccs_get_col(in_real, i, &i_r, &v_r);
            const JMTX_INDEX_T c_i = JMTX_NAME_TYPED(matrix_ccs_get_col(in_imag, i, &i_i, &v_i);
            JMTX_INDEX_T k_r, k_i;
            for (k_r = 0, k_i = 0; k_r < c_r && k_i < c_i; ++j)
            {
                if (i_r[k_r] < i_i[k_i])
                {
                    values[j] = v_r[k_r];
                    indices[j] = i_r[k_r];
                    k_r += 1;
                }
                else if (i_r[k_r] > i_i[k_i])
                {
                    values[j] = v_i[k_i] * _Complex_I;
                    indices[j] = i_i[k_r];
                    k_i += 1;
                }
                else // if (i_r[k_r] == i_i[k_i])
                {
                    values[j] = v_i[k_i] * _Complex_I + v_r[k_r];
                    indices[j] = i_i[k_r];
                    k_r += 1;
                    k_i += 1;
                }
            }

            if (k_i < c_i)
            {
                assert(k_r == c_r);
                while (k_i < c_i)
                {
                    values[j] = v_i[k_i] * _Complex_I;
                    indices[j] = i_i[k_i];
                    j += 1;
                    k_i += 1;
                }
            }
            else
            {
                assert(k_i == c_i);
                while (k_r < c_r)
                {
                    values[j] = v_r[k_r];
                    indices[j] = i_r[k_r];
                    j += 1;
                    k_r += 1;
                }
            }

            offsets[i] = j;
        }

        mtx->indices = indices;
        mtx->values = values;
        mtx->capacity = needed_capacity;
        mtx->n_entries = j;
        mtx->end_of_column_offsets = offsets;
    }
    else if (in_imag)
    {
        _Complex JMTX_SCALAR_T *values =
            allocator_callbacks->alloc(allocator_callbacks->state, (in_imag->n_entries) * sizeof(*values));
        if (!values)
        {
            allocator_callbacks->free(allocator_callbacks->state, mtx);
            return JMTX_RESULT_BAD_ALLOC;
        }

        indices = allocator_callbacks->alloc(allocator_callbacks->state, (in_imag->n_entries) * sizeof(*indices));
        if (!indices)
        {
            allocator_callbacks->free(allocator_callbacks->state, indices);
            allocator_callbacks->free(allocator_callbacks->state, mtx);
            return JMTX_RESULT_BAD_ALLOC;
        }

        offsets = allocator_callbacks->alloc(allocator_callbacks->state, (in_imag->base.cols) * sizeof(*offsets));
        if (!offsets)
        {
            allocator_callbacks->free(allocator_callbacks->state, offsets);
            allocator_callbacks->free(allocator_callbacks->state, indices);
            allocator_callbacks->free(allocator_callbacks->state, mtx);
            return JMTX_RESULT_BAD_ALLOC;
        }
        memcpy(offsets, in_imag->end_of_column_offsets, (in_imag->base.cols) * sizeof(*offsets));
        memcpy(indices, in_imag->indices, (in_imag->n_entries) * sizeof(*indices));

        for (JMTX_FAST_INT_T i = 0; i < in_imag->n_entries; ++i)
        {
            values[i] = _Complex_I * in_imag->values[i];
        }

        mtx->indices = indices;
        mtx->values = values;
        mtx->capacity = in_imag->n_entries;
        mtx->n_entries = in_imag->n_entries;
        mtx->end_of_column_offsets = offsets;
    }
    else // if (in_real)
    {
        _Complex JMTX_SCALAR_T *values =
            allocator_callbacks->alloc(allocator_callbacks->state, (in_real->n_entries) * sizeof(*values));
        if (!values)
        {
            allocator_callbacks->free(allocator_callbacks->state, mtx);
            return JMTX_RESULT_BAD_ALLOC;
        }

        indices = allocator_callbacks->alloc(allocator_callbacks->state, (in_real->n_entries) * sizeof(*indices));
        if (!indices)
        {
            allocator_callbacks->free(allocator_callbacks->state, indices);
            allocator_callbacks->free(allocator_callbacks->state, mtx);
            return JMTX_RESULT_BAD_ALLOC;
        }

        offsets = allocator_callbacks->alloc(allocator_callbacks->state, (in_real->base.cols) * sizeof(*offsets));
        if (!offsets)
        {
            allocator_callbacks->free(allocator_callbacks->state, offsets);
            allocator_callbacks->free(allocator_callbacks->state, indices);
            allocator_callbacks->free(allocator_callbacks->state, mtx);
            return JMTX_RESULT_BAD_ALLOC;
        }
        memcpy(offsets, in_real->end_of_column_offsets, (in_real->base.cols) * sizeof(*offsets));
        memcpy(indices, in_real->indices, (in_real->n_entries) * sizeof(*indices));

        for (JMTX_FAST_INT_T i = 0; i < in_real->n_entries; ++i)
        {
            values[i] = _Complex_I * in_real->values[i];
        }

        mtx->indices = indices;
        mtx->values = values;
        mtx->capacity = in_real->n_entries;
        mtx->n_entries = in_real->n_entries;
        mtx->end_of_column_offsets = offsets;
    }

    mtx->base.cols = c;
    mtx->base.type = JMTXZ_TYPE_CCS;
    mtx->base.rows = r;
    mtx->base.allocator_callbacks = *allocator_callbacks;
    *p_mtx = mtx;

    return JMTX_RESULT_SUCCESS;
}

/**
 * Creates a new complex CCS matrix with JMTX_SCALAR_T precision from a CCS matrix with JMTX_SCALAR_T precision as its real part.
 * Only one memory reallocation.
 * may be needed. Can not fail if the input matrix is valid.
 * @param in matrix which to convert
 * @return converted matrix, or NULL in case of allocation failure
 */
jmtxz_matrix_ccs *jmtxz_matrix_ccs_from_double_real_inplace(JMTX_NAME_TYPED(matrix_ccs) *in)
{
    _Complex JMTX_SCALAR_T *vals = NULL;
    JMTX_INDEX_T capacity_for_complex = in->capacity * sizeof(*in->values) / sizeof(*vals);
    if (capacity_for_complex < in->n_entries)
    {
        capacity_for_complex = in->n_entries;
        vals = in->base.allocator_callbacks.realloc(in->base.allocator_callbacks.state, in->values,
                                                    sizeof(*vals) * capacity_for_complex);
        if (!vals)
        {
            return NULL;
        }
        in->values = (JMTX_SCALAR_T *)vals;
    }

    vals = (_Complex JMTX_SCALAR_T *)in->values;
    in->capacity = capacity_for_complex;

    for (JMTX_INDEX_T i = 0; i < in->n_entries; ++i)
    {
        vals[in->n_entries - 1 - i] = (_Complex double)in->values[in->n_entries - 1 - i];
    }
    in->base.type = JMTXZ_TYPE_CCS;

    return (jmtxz_matrix_ccs *)in;
}

/**
 * Creates a new CCS matrix with JMTX_SCALAR_T precision from a CCS matrix with JMTX_SCALAR_T precision as its imaginary part.
 * Only one memory reallocation.
 * may be needed. Can not fail if the input matrix is valid.
 * @param in matrix which to convert
 * @return converted matrix, or NULL in case of allocation failure
 */
jmtxz_matrix_ccs *jmtxz_matrix_ccs_from_double_imag_inplace(JMTX_NAME_TYPED(matrix_ccs) *in)
{
    _Complex JMTX_SCALAR_T *vals = NULL;
    JMTX_INDEX_T capacity_for_complex = in->capacity * sizeof(*in->values) / sizeof(*vals);
    if (capacity_for_complex < in->n_entries)
    {
        capacity_for_complex = in->n_entries;
        vals = in->base.allocator_callbacks.realloc(in->base.allocator_callbacks.state, in->values,
                                                    sizeof(*vals) * capacity_for_complex);
        if (!vals)
        {
            return NULL;
        }
        in->values = (JMTX_SCALAR_T *)vals;
    }

    vals = (_Complex JMTX_SCALAR_T *)in->values;
    in->capacity = capacity_for_complex;

    for (JMTX_INDEX_T i = 0; i < in->n_entries; ++i)
    {
        vals[in->n_entries - 1 - i] = (_Complex double)in->values[in->n_entries - 1 - i] * _Complex_I;
    }
    in->base.type = JMTXZ_TYPE_CCS;

    return (jmtxz_matrix_ccs *)in;
}

/**
 * Creates a new complex CCS matrix with JMTX_SCALAR_T precision from a CCS matrix with JMTX_SCALAR_T precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in_real matrix to use as the real component of the new matrix (if NULL real part is zero)
 * @param in_imag matrix to use as the real component of the new matrix (if NULL imaginary part is zero)
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxzs_matrix_ccs_from_double(jmtxz_matrix_ccs **p_mtx, const JMTX_NAME_TYPED(matrix_ccs) *in_real,
                                          const JMTX_NAME_TYPED(matrix_ccs) *in_imag,
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
    if ((in_real && in_real->base.type != JMTXD_TYPE_CCS) || (in_imag && in_imag->base.type != JMTXD_TYPE_CCS))
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if ((in_real && in_imag) && (in_real->base.rows != in_imag->base.rows || in_real->base.cols != in_imag->base.cols))
    {
        return JMTX_RESULT_BAD_MATRIX;
    }

    return jmtxz_matrix_ccs_from_double(p_mtx, in_real, in_imag, allocator_callbacks);
}

/***********************************************************************************************************************
 *                                                                                                                     *
 *                                 COMPLEX FLOAT <->  COMPLEX DOUBLE                                                   *
 *                                                                                                                     *
 **********************************************************************************************************************/
/**
 * Creates a new complex CCS matrix with single precision from a complex CCS matrix with JMTX_SCALAR_T precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxc_matrix_ccs_from_cdouble(jmtxc_matrix_ccs **p_mtx, const jmtxz_matrix_ccs *in,
                                          const jmtx_allocator_callbacks *allocator_callbacks)
{
    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }

    JMTX_INDEX_T *offsets = NULL;
    JMTX_INDEX_T *indices = NULL;

    jmtxc_matrix_ccs *mtx = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*mtx));
    if (!mtx)
    {
        return JMTX_RESULT_BAD_ALLOC;
    }

    _Complex float *values = allocator_callbacks->alloc(allocator_callbacks->state, (in->n_entries) * sizeof(*values));
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
        values[i] = (_Complex float)in->values[i];
    }

    mtx->base.cols = in->base.cols;
    mtx->base.type = JMTXC_TYPE_CCS;
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
 * Creates a new complex CCS matrix with single precision from a complex CCS matrix with JMTX_SCALAR_T precision. Requires no
 * memory allocation by reusing the memory of the initial matrix. Can not fail if the input matrix is valid.
 * @param in matrix which to convert (will be invalid if function succeeds)
 * @return converted matrix
 */
jmtxc_matrix_ccs *jmtxc_matrix_ccs_from_cdouble_inplace(jmtxz_matrix_ccs *in)
{
    _Complex float *const values = (_Complex float *)in->values;

    for (JMTX_FAST_INT_T i = 0; i < in->n_entries; ++i)
    {
        values[i] = (_Complex float)in->values[i];
    }

    in->capacity = (in->capacity * sizeof(*in->values)) / sizeof(*values);
    in->base.type = JMTXC_TYPE_CCS;

    return (jmtxc_matrix_ccs *)in;
}

/**
 * Creates a new complex CCS matrix with single precision from a complex CCS matrix with JMTX_SCALAR_T precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxcs_matrix_ccs_from_cdouble(jmtxc_matrix_ccs **p_mtx, const jmtxz_matrix_ccs *in,
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
    if (in->base.type != JMTXZ_TYPE_CCS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }

    return jmtxc_matrix_ccs_from_cdouble(p_mtx, in, allocator_callbacks);
}

/**
 * Creates a new complex CCS matrix with JMTX_SCALAR_T precision from a complex CCS matrix with single precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxz_matrix_ccs_from_cfloat(jmtxz_matrix_ccs **p_mtx, const jmtxc_matrix_ccs *in,
                                         const jmtx_allocator_callbacks *allocator_callbacks)
{
    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }

    JMTX_INDEX_T *offsets = NULL;
    JMTX_INDEX_T *indices = NULL;

    jmtxz_matrix_ccs *mtx = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*mtx));
    if (!mtx)
    {
        return JMTX_RESULT_BAD_ALLOC;
    }

    _Complex JMTX_SCALAR_T *values =
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
        values[i] = (_Complex double)in->values[i];
    }

    mtx->base.cols = in->base.cols;
    mtx->base.type = JMTXZ_TYPE_CCS;
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
 * Creates a new complex CCS matrix with JMTX_SCALAR_T precision from a complex CCS matrix with single precision. Only one
 * memory reallocation may be needed. Can not fail if the input matrix is valid.
 * @param in matrix which to convert
 * @return converted matrix, or NULL in case of allocation failure
 */
jmtxz_matrix_ccs *jmtxz_matrix_ccs_from_cfloat_inplace(jmtxc_matrix_ccs *in)
{
    _Complex JMTX_SCALAR_T *vals = NULL;
    JMTX_INDEX_T capacity_for_doubles = in->capacity * sizeof(*in->values) / sizeof(*vals);
    if (capacity_for_doubles < in->n_entries)
    {
        capacity_for_doubles = in->n_entries;
        vals = in->base.allocator_callbacks.realloc(in->base.allocator_callbacks.state, in->values,
                                                    sizeof(*vals) * capacity_for_doubles);
        if (!vals)
        {
            return NULL;
        }
        in->values = (_Complex float *)vals;
    }

    vals = (_Complex JMTX_SCALAR_T *)in->values;
    in->capacity = capacity_for_doubles;

    for (JMTX_INDEX_T i = 0; i < in->n_entries; ++i)
    {
        vals[in->n_entries - 1 - i] = (_Complex double)in->values[in->n_entries - 1 - i];
    }
    in->base.type = JMTXD_TYPE_CCS;

    return (jmtxz_matrix_ccs *)in;
}

/**
 * Creates a new complex CCS matrix with JMTX_SCALAR_T precision from a complex CCS matrix with single precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxzs_matrix_ccs_from_cfloat(jmtxz_matrix_ccs **p_mtx, const jmtxc_matrix_ccs *in,
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
    if (in->base.type != JMTXC_TYPE_CCS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }

    return jmtxz_matrix_ccs_from_cfloat(p_mtx, in, allocator_callbacks);
}

#endif //!_MSC_BUILD
