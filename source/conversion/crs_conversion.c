//
// Created by jan on 1.12.2023.
//

#include <complex.h>
#include <assert.h>

#include "../matrix_base_internal.h"
#include "../float/matrices/sparse_row_compressed_internal.h"
#include "../double/matrices/sparse_row_compressed_internal.h"
#include "../cfloat/matrices/sparse_row_compressed_internal.h"
#include "../cdouble/matrices/sparse_row_compressed_internal.h"
#include "../../include/jmtx/conversion/crs_conversion.h"

/***********************************************************************************************************************
 *                                                                                                                     *
 *                                          FLOAT <-> DOUBLE                                                           *
 *                                                                                                                     *
 **********************************************************************************************************************/

/**
 * Creates a new CRS matrix with single precision from a CRS matrix with double precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtx_matrix_crs_from_double(jmtx_matrix_crs** p_mtx, const jmtxd_matrix_crs* in,
                                        const jmtx_allocator_callbacks* allocator_callbacks)
{
    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }

    uint32_t* offsets = NULL;
    uint32_t* indices = NULL;

    jmtx_matrix_crs* mtx = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*mtx));
    if (!mtx)
    {
        return JMTX_RESULT_BAD_ALLOC;
    }

    float* values = allocator_callbacks->alloc(allocator_callbacks->state, (in->n_entries) * sizeof(*values));
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

    offsets = allocator_callbacks->alloc(allocator_callbacks->state, (in->base.rows) * sizeof(*offsets));
    if (!offsets)
    {
        allocator_callbacks->free(allocator_callbacks->state, offsets);
        allocator_callbacks->free(allocator_callbacks->state, indices);
        allocator_callbacks->free(allocator_callbacks->state, mtx);
        return JMTX_RESULT_BAD_ALLOC;
    }
    memcpy(offsets, in->end_of_row_offsets, (in->base.rows) * sizeof(*offsets));
    memcpy(indices, in->indices, (in->n_entries) * sizeof(*indices));

    for (uint_fast32_t i = 0; i < in->n_entries; ++i)
    {
        values[i] = (float)in->values[i];
    }

    mtx->base.cols = in->base.cols;
    mtx->base.type = JMTX_TYPE_CRS;
    mtx->base.rows = in->base.rows;
    mtx->base.allocator_callbacks = *allocator_callbacks;
    mtx->indices = indices;
    mtx->values = values;
    mtx->capacity = in->n_entries;
    mtx->n_entries = in->n_entries;
    mtx->end_of_row_offsets = offsets;
    *p_mtx = mtx;


    return JMTX_RESULT_SUCCESS;
}

/**
 * Creates a new CRS matrix with single precision from a CRS matrix with double precision. Requires minimum memory
 * allocation by reusing the memory of the initial matrix.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert (will be invalid if function succeeds)
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_matrix_crs* jmtx_matrix_crs_from_double_inplace(jmtxd_matrix_crs* in)
{
    float* const values = (float*)in->values;

    for (uint_fast32_t i = 0; i < in->n_entries; ++i)
    {
        values[i] = (float)in->values[i];
    }

    in->capacity = (in->capacity * sizeof(*in->values)) / sizeof(*values);
    in->base.type = JMTX_TYPE_CRS;

    return (jmtx_matrix_crs*)in;
}

/**
 * Creates a new CRS matrix with single precision from a CRS matrix with double precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxs_matrix_crs_from_double(jmtx_matrix_crs** p_mtx, const jmtxd_matrix_crs* in,
                                         const jmtx_allocator_callbacks* allocator_callbacks)
{
    if (!p_mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!in)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (allocator_callbacks && (!allocator_callbacks->free || !allocator_callbacks->alloc || !allocator_callbacks->realloc))
    {
        return JMTX_RESULT_BAD_PARAM;
    }
    if (in->base.type != JMTXD_TYPE_CRS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }

    return jmtx_matrix_crs_from_double(p_mtx, in, allocator_callbacks);
}

/**
 * Creates a new CRS matrix with double precision from a CRS matrix with single precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxd_matrix_crs_from_float(jmtxd_matrix_crs** p_mtx, const jmtx_matrix_crs* in,
                                        const jmtx_allocator_callbacks* allocator_callbacks)
{
    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }

    uint32_t* offsets = NULL;
    uint32_t* indices = NULL;

    jmtxd_matrix_crs* mtx = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*mtx));
    if (!mtx)
    {
        return JMTX_RESULT_BAD_ALLOC;
    }

    double* values = allocator_callbacks->alloc(allocator_callbacks->state, (in->n_entries) * sizeof(*values));
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

    offsets = allocator_callbacks->alloc(allocator_callbacks->state, (in->base.rows) * sizeof(*offsets));
    if (!offsets)
    {
        allocator_callbacks->free(allocator_callbacks->state, offsets);
        allocator_callbacks->free(allocator_callbacks->state, indices);
        allocator_callbacks->free(allocator_callbacks->state, mtx);
        return JMTX_RESULT_BAD_ALLOC;
    }
    memcpy(offsets, in->end_of_row_offsets, (in->base.rows) * sizeof(*offsets));
    memcpy(indices, in->indices, (in->n_entries) * sizeof(*indices));

    for (uint_fast32_t i = 0; i < in->n_entries; ++i)
    {
        values[i] = (double)in->values[i];
    }

    mtx->base.cols = in->base.cols;
    mtx->base.type = JMTXD_TYPE_CRS;
    mtx->base.rows = in->base.rows;
    mtx->base.allocator_callbacks = *allocator_callbacks;
    mtx->indices = indices;
    mtx->values = values;
    mtx->capacity = in->n_entries;
    mtx->n_entries = in->n_entries;
    mtx->end_of_row_offsets = offsets;
    *p_mtx = mtx;

    return JMTX_RESULT_SUCCESS;
}

/**
 * Creates a new CRS matrix with double precision from a CRS matrix with single precision. Only one memory reallocation
 * may be needed. Can not fail if the input matrix is valid.
 * @param in matrix which to convert
 * @return converted matrix, or NULL in case of allocation failure
 */
jmtxd_matrix_crs* jmtxd_matrix_crs_from_float_inplace(jmtx_matrix_crs* in)
{
    double* vals = NULL;
    uint32_t capacity_for_doubles = in->capacity * sizeof(*in->values) / sizeof(*vals);
    if (capacity_for_doubles < in->n_entries)
    {
        capacity_for_doubles = in->n_entries;
        vals = in->base.allocator_callbacks.realloc(in->base.allocator_callbacks.state, in->values, sizeof(*vals) * capacity_for_doubles);
        if (!vals)
        {
            return NULL;
        }
        in->values = (float*)vals;
    }

    vals = (double*)in->values;
    in->capacity = capacity_for_doubles;

    for (uint32_t i = 0; i < in->n_entries; ++i)
    {
        vals[in->n_entries - 1 - i] = (double)in->values[in->n_entries - 1 - i];
    }
    in->base.type = JMTXD_TYPE_CRS;

    return (jmtxd_matrix_crs*)in;
}

/**
 * Creates a new CRS matrix with double precision from a CRS matrix with single precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxds_matrix_crs_from_float(jmtxd_matrix_crs** p_mtx, const jmtx_matrix_crs* in,
                                         const jmtx_allocator_callbacks* allocator_callbacks)
{
    if (!p_mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!in)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (allocator_callbacks && (!allocator_callbacks->free || !allocator_callbacks->alloc || !allocator_callbacks->realloc))
    {
        return JMTX_RESULT_BAD_PARAM;
    }
    if (in->base.type != JMTX_TYPE_CRS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }

    return jmtxd_matrix_crs_from_float(p_mtx, in, allocator_callbacks);
}


/***********************************************************************************************************************
 *                                                                                                                     *
 *                                          FLOAT <-> COMPLEX FLOAT                                                    *
 *                                                                                                                     *
 **********************************************************************************************************************/

/**
 * Creates a new CRS matrix with single precision from the real part of a complex CRS matrix with single precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtx_matrix_crs_from_cfloat_real(jmtx_matrix_crs** p_mtx, const jmtxc_matrix_crs* in,
                                             const jmtx_allocator_callbacks* allocator_callbacks)
{
    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }

    uint32_t* offsets = NULL;
    uint32_t* indices = NULL;

    jmtx_matrix_crs* mtx = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*mtx));
    if (!mtx)
    {
        return JMTX_RESULT_BAD_ALLOC;
    }

    float* values = allocator_callbacks->alloc(allocator_callbacks->state, (in->n_entries) * sizeof(*values));
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

    offsets = allocator_callbacks->alloc(allocator_callbacks->state, (in->base.rows) * sizeof(*offsets));
    if (!offsets)
    {
        allocator_callbacks->free(allocator_callbacks->state, offsets);
        allocator_callbacks->free(allocator_callbacks->state, indices);
        allocator_callbacks->free(allocator_callbacks->state, mtx);
        return JMTX_RESULT_BAD_ALLOC;
    }
    memcpy(offsets, in->end_of_row_offsets, (in->base.rows) * sizeof(*offsets));
    memcpy(indices, in->indices, (in->n_entries) * sizeof(*indices));

    for (uint_fast32_t i = 0; i < in->n_entries; ++i)
    {
        values[i] = crealf(in->values[i]);
    }

    mtx->base.cols = in->base.cols;
    mtx->base.type = JMTX_TYPE_CRS;
    mtx->base.rows = in->base.rows;
    mtx->base.allocator_callbacks = *allocator_callbacks;
    mtx->indices = indices;
    mtx->values = values;
    mtx->capacity = in->n_entries;
    mtx->n_entries = in->n_entries;
    mtx->end_of_row_offsets = offsets;
    *p_mtx = mtx;


    return JMTX_RESULT_SUCCESS;
}
/**
 * Creates a new CRS matrix with single precision from the imaginary part of a complex CRS matrix with single precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtx_matrix_crs_from_cfloat_imag(jmtx_matrix_crs** p_mtx, const jmtxc_matrix_crs* in,
                                             const jmtx_allocator_callbacks* allocator_callbacks)
{
    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }

    uint32_t* offsets = NULL;
    uint32_t* indices = NULL;

    jmtx_matrix_crs* mtx = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*mtx));
    if (!mtx)
    {
        return JMTX_RESULT_BAD_ALLOC;
    }

    float* values = allocator_callbacks->alloc(allocator_callbacks->state, (in->n_entries) * sizeof(*values));
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

    offsets = allocator_callbacks->alloc(allocator_callbacks->state, (in->base.rows) * sizeof(*offsets));
    if (!offsets)
    {
        allocator_callbacks->free(allocator_callbacks->state, offsets);
        allocator_callbacks->free(allocator_callbacks->state, indices);
        allocator_callbacks->free(allocator_callbacks->state, mtx);
        return JMTX_RESULT_BAD_ALLOC;
    }
    memcpy(offsets, in->end_of_row_offsets, (in->base.rows) * sizeof(*offsets));
    memcpy(indices, in->indices, (in->n_entries) * sizeof(*indices));

    for (uint_fast32_t i = 0; i < in->n_entries; ++i)
    {
        values[i] = cimagf(in->values[i]);
    }

    mtx->base.cols = in->base.cols;
    mtx->base.type = JMTX_TYPE_CRS;
    mtx->base.rows = in->base.rows;
    mtx->base.allocator_callbacks = *allocator_callbacks;
    mtx->indices = indices;
    mtx->values = values;
    mtx->capacity = in->n_entries;
    mtx->n_entries = in->n_entries;
    mtx->end_of_row_offsets = offsets;
    *p_mtx = mtx;


    return JMTX_RESULT_SUCCESS;
}

/**
 * Creates a new CRS matrix with single precision from the real part of a complex CRS matrix with single precision.
 * Requires no memory allocation by reusing the memory of the initial matrix. Can not fail if the input matrix is valid.
 * @param in matrix which to convert (will be invalid if function succeeds)
 * @return converted matrix
 */
jmtx_matrix_crs* jmtx_matrix_crs_from_cfloat_real_inplace(jmtxc_matrix_crs* in)
{
    float* const values = (float*)in->values;

    for (uint_fast32_t i = 0; i < in->n_entries; ++i)
    {
        values[i] = crealf(in->values[i]);
    }

    in->capacity = (in->capacity * sizeof(*in->values)) / sizeof(*values);
    in->base.type = JMTX_TYPE_CRS;

    return (jmtx_matrix_crs*)in;
}

/**
 * Creates a new CRS matrix with single precision from the imaginary part of a complex CRS matrix with single precision.
 * Requires no memory allocation by reusing the memory of the initial matrix. Can not fail if the input matrix is valid.
 * @param in matrix which to convert (will be invalid if function succeeds)
 * @return converted matrix
 */
jmtx_matrix_crs* jmtx_matrix_crs_from_cfloat_imag_inplace(jmtxc_matrix_crs* in)
{
    float* const values = (float*)in->values;

    for (uint_fast32_t i = 0; i < in->n_entries; ++i)
    {
        values[i] = cimagf(in->values[i]);
    }

    in->capacity = (in->capacity * sizeof(*in->values)) / sizeof(*values);
    in->base.type = JMTX_TYPE_CRS;

    return (jmtx_matrix_crs*)in;
}

/**
 * Creates a new CRS matrix with single precision from the real part of a complex CRS matrix with single precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxs_matrix_crs_from_cfloat_real(jmtx_matrix_crs** p_mtx, const jmtxc_matrix_crs* in,
                                              const jmtx_allocator_callbacks* allocator_callbacks)
{
    if (!p_mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!in)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (allocator_callbacks && (!allocator_callbacks->free || !allocator_callbacks->alloc || !allocator_callbacks->realloc))
    {
        return JMTX_RESULT_BAD_PARAM;
    }
    if (in->base.type != JMTXC_TYPE_CRS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }

    return jmtx_matrix_crs_from_cfloat_real(p_mtx, in, allocator_callbacks);
}

/**
 * Creates a new CRS matrix with single precision from the imaginary part of a complex CRS matrix with single precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxs_matrix_crs_from_cfloat_imag(jmtx_matrix_crs** p_mtx, const jmtxc_matrix_crs* in,
                                              const jmtx_allocator_callbacks* allocator_callbacks)
{
    if (!p_mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!in)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (allocator_callbacks && (!allocator_callbacks->free || !allocator_callbacks->alloc || !allocator_callbacks->realloc))
    {
        return JMTX_RESULT_BAD_PARAM;
    }
    if (in->base.type != JMTXC_TYPE_CRS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }

    return jmtx_matrix_crs_from_cfloat_imag(p_mtx, in, allocator_callbacks);
}

/**
 * Creates a new complex CRS matrix with single precision from a CRS matrix with single precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in_real matrix to use as the real component of the new matrix (if NULL real part is zero)
 * @param in_imag matrix to use as the real component of the new matrix (if NULL imaginary part is zero)
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxc_matrix_crs_from_float(jmtxc_matrix_crs** p_mtx, const jmtx_matrix_crs* in_real,
                                        const jmtx_matrix_crs* in_imag,
                                        const jmtx_allocator_callbacks* allocator_callbacks)
{
    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }

    uint32_t* offsets = NULL;
    uint32_t* indices = NULL;
    uint32_t r, c;

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
    jmtxc_matrix_crs* mtx = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*mtx));
    if (!mtx)
    {
        return JMTX_RESULT_BAD_ALLOC;
    }

    if (in_real && in_imag)
    {
        const uint32_t needed_capacity = in_real->n_entries + in_imag->n_entries;


        _Complex float* values = allocator_callbacks->alloc(allocator_callbacks->state, (needed_capacity) * sizeof(*values));
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

        offsets = allocator_callbacks->alloc(allocator_callbacks->state, (r) * sizeof(*offsets));
        if (!offsets)
        {
            allocator_callbacks->free(allocator_callbacks->state, offsets);
            allocator_callbacks->free(allocator_callbacks->state, indices);
            allocator_callbacks->free(allocator_callbacks->state, mtx);
            return JMTX_RESULT_BAD_ALLOC;
        }
        uint32_t i, j;
        for (i = 0, j = 0; i < r; ++i)
        {
            float* v_r, *v_i;
            uint32_t* i_r, *i_i;
            const uint32_t c_r = jmtx_matrix_crs_get_row(in_real, i, &i_r, &v_r);
            const uint32_t c_i = jmtx_matrix_crs_get_row(in_imag, i, &i_i, &v_i);
            uint32_t k_r, k_i;
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
                else //if (i_r[k_r] == i_i[k_i])
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
        mtx->end_of_row_offsets = offsets;

    }
    else if (in_imag)
    {
        _Complex float* values = allocator_callbacks->alloc(allocator_callbacks->state, (in_imag->n_entries) * sizeof(*values));
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

        offsets = allocator_callbacks->alloc(allocator_callbacks->state, (in_imag->base.rows) * sizeof(*offsets));
        if (!offsets)
        {
            allocator_callbacks->free(allocator_callbacks->state, offsets);
            allocator_callbacks->free(allocator_callbacks->state, indices);
            allocator_callbacks->free(allocator_callbacks->state, mtx);
            return JMTX_RESULT_BAD_ALLOC;
        }
        memcpy(offsets, in_imag->end_of_row_offsets, (in_imag->base.rows) * sizeof(*offsets));
        memcpy(indices, in_imag->indices, (in_imag->n_entries) * sizeof(*indices));

        for (uint_fast32_t i = 0; i < in_imag->n_entries; ++i)
        {
            values[i] = _Complex_I * in_imag->values[i];
        }

        mtx->indices = indices;
        mtx->values = values;
        mtx->capacity = in_imag->n_entries;
        mtx->n_entries = in_imag->n_entries;
        mtx->end_of_row_offsets = offsets;
    }
    else // if (in_real)
    {
        _Complex float* values = allocator_callbacks->alloc(allocator_callbacks->state, (in_real->n_entries) * sizeof(*values));
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

        offsets = allocator_callbacks->alloc(allocator_callbacks->state, (in_real->base.rows) * sizeof(*offsets));
        if (!offsets)
        {
            allocator_callbacks->free(allocator_callbacks->state, offsets);
            allocator_callbacks->free(allocator_callbacks->state, indices);
            allocator_callbacks->free(allocator_callbacks->state, mtx);
            return JMTX_RESULT_BAD_ALLOC;
        }
        memcpy(offsets, in_real->end_of_row_offsets, (in_real->base.rows) * sizeof(*offsets));
        memcpy(indices, in_real->indices, (in_real->n_entries) * sizeof(*indices));

        for (uint_fast32_t i = 0; i < in_real->n_entries; ++i)
        {
            values[i] = _Complex_I * in_real->values[i];
        }

        mtx->indices = indices;
        mtx->values = values;
        mtx->capacity = in_real->n_entries;
        mtx->n_entries = in_real->n_entries;
        mtx->end_of_row_offsets = offsets;

    }
    
    mtx->base.cols = c;
    mtx->base.type = JMTXC_TYPE_CRS;
    mtx->base.rows = r;
    mtx->base.allocator_callbacks = *allocator_callbacks;
    *p_mtx = mtx;

    return JMTX_RESULT_SUCCESS;
}

/**
 * Creates a new complex CRS matrix with single precision from a CRS matrix with single precision as its real part.
 * Only one memory reallocation.
 * may be needed. Can not fail if the input matrix is valid.
 * @param in matrix which to convert
 * @return converted matrix, or NULL in case of allocation failure
 */
jmtxc_matrix_crs* jmtxc_matrix_crs_from_float_real_inplace(jmtx_matrix_crs* in)
{
    _Complex float* vals = NULL;
    uint32_t capacity_for_complex = in->capacity * sizeof(*in->values) / sizeof(*vals);
    if (capacity_for_complex < in->n_entries)
    {
        capacity_for_complex = in->n_entries;
        vals = in->base.allocator_callbacks.realloc(in->base.allocator_callbacks.state, in->values, sizeof(*vals) * capacity_for_complex);
        if (!vals)
        {
            return NULL;
        }
        in->values = (float*)vals;
    }

    vals = (_Complex float*)in->values;
    in->capacity = capacity_for_complex;

    for (uint32_t i = 0; i < in->n_entries; ++i)
    {
        vals[in->n_entries - 1 - i] = (_Complex float)in->values[in->n_entries - 1 - i];
    }
    in->base.type = JMTXC_TYPE_CRS;

    return (jmtxc_matrix_crs*)in;
}

/**
 * Creates a new CRS matrix with single precision from a CRS matrix with single precision as its imaginary part.
 * Only one memory reallocation.
 * may be needed. Can not fail if the input matrix is valid.
 * @param in matrix which to convert
 * @return converted matrix, or NULL in case of allocation failure
 */
jmtxc_matrix_crs* jmtxc_matrix_crs_from_float_imag_inplace(jmtx_matrix_crs* in)
{
    _Complex float* vals = NULL;
    uint32_t capacity_for_complex = in->capacity * sizeof(*in->values) / sizeof(*vals);
    if (capacity_for_complex < in->n_entries)
    {
        capacity_for_complex = in->n_entries;
        vals = in->base.allocator_callbacks.realloc(in->base.allocator_callbacks.state, in->values, sizeof(*vals) * capacity_for_complex);
        if (!vals)
        {
            return NULL;
        }
        in->values = (float*)vals;
    }

    vals = (_Complex float*)in->values;
    in->capacity = capacity_for_complex;

    for (uint32_t i = 0; i < in->n_entries; ++i)
    {
        vals[in->n_entries - 1 - i] = (_Complex float)in->values[in->n_entries - 1 - i] * _Complex_I;
    }
    in->base.type = JMTXC_TYPE_CRS;

    return (jmtxc_matrix_crs*)in;
}

/**
 * Creates a new complex CRS matrix with single precision from a CRS matrix with single precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in_real matrix to use as the real component of the new matrix (if NULL real part is zero)
 * @param in_imag matrix to use as the real component of the new matrix (if NULL imaginary part is zero)
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxcs_matrix_crs_from_float(jmtxc_matrix_crs** p_mtx, const jmtx_matrix_crs* in_real,
                                         const jmtx_matrix_crs* in_imag,
                                         const jmtx_allocator_callbacks* allocator_callbacks)
{
    if (!p_mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!in_real && !in_imag)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (allocator_callbacks && (!allocator_callbacks->free || !allocator_callbacks->alloc || !allocator_callbacks->realloc))
    {
        return JMTX_RESULT_BAD_PARAM;
    }
    if ((in_real && in_real->base.type != JMTX_TYPE_CRS) ||
        (in_imag && in_imag->base.type != JMTX_TYPE_CRS))
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if ((in_real && in_imag) && (in_real->base.rows != in_imag->base.rows || in_real->base.cols != in_imag->base.cols))
    {
        return JMTX_RESULT_BAD_MATRIX;
    }

    return jmtxc_matrix_crs_from_float(p_mtx, in_real, in_imag, allocator_callbacks);
}

/***********************************************************************************************************************
 *                                                                                                                     *
 *                                          DOUBLE <-> COMPLEX DOUBLE                                                  *
 *                                                                                                                     *
 **********************************************************************************************************************/

/**
 * Creates a new CRS matrix with double precision from the real part of a complex CRS matrix with double precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxd_matrix_crs_from_cdouble_real(jmtxd_matrix_crs** p_mtx, const jmtxz_matrix_crs* in,
                                               const jmtx_allocator_callbacks* allocator_callbacks)
{
    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }

    uint32_t* offsets = NULL;
    uint32_t* indices = NULL;

    jmtxd_matrix_crs* mtx = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*mtx));
    if (!mtx)
    {
        return JMTX_RESULT_BAD_ALLOC;
    }

    double* values = allocator_callbacks->alloc(allocator_callbacks->state, (in->n_entries) * sizeof(*values));
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

    offsets = allocator_callbacks->alloc(allocator_callbacks->state, (in->base.rows) * sizeof(*offsets));
    if (!offsets)
    {
        allocator_callbacks->free(allocator_callbacks->state, offsets);
        allocator_callbacks->free(allocator_callbacks->state, indices);
        allocator_callbacks->free(allocator_callbacks->state, mtx);
        return JMTX_RESULT_BAD_ALLOC;
    }
    memcpy(offsets, in->end_of_row_offsets, (in->base.rows) * sizeof(*offsets));
    memcpy(indices, in->indices, (in->n_entries) * sizeof(*indices));

    for (uint_fast32_t i = 0; i < in->n_entries; ++i)
    {
        values[i] = creal(in->values[i]);
    }

    mtx->base.cols = in->base.cols;
    mtx->base.type = JMTXD_TYPE_CRS;
    mtx->base.rows = in->base.rows;
    mtx->base.allocator_callbacks = *allocator_callbacks;
    mtx->indices = indices;
    mtx->values = values;
    mtx->capacity = in->n_entries;
    mtx->n_entries = in->n_entries;
    mtx->end_of_row_offsets = offsets;
    *p_mtx = mtx;


    return JMTX_RESULT_SUCCESS;
}
/**
 * Creates a new CRS matrix with double precision from the imaginary part of a complex CRS matrix with double precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxd_matrix_crs_from_cdouble_imag(jmtxd_matrix_crs** p_mtx, const jmtxz_matrix_crs* in,
                                               const jmtx_allocator_callbacks* allocator_callbacks)
{
    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }

    uint32_t* offsets = NULL;
    uint32_t* indices = NULL;

    jmtxd_matrix_crs* mtx = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*mtx));
    if (!mtx)
    {
        return JMTX_RESULT_BAD_ALLOC;
    }

    double* values = allocator_callbacks->alloc(allocator_callbacks->state, (in->n_entries) * sizeof(*values));
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

    offsets = allocator_callbacks->alloc(allocator_callbacks->state, (in->base.rows) * sizeof(*offsets));
    if (!offsets)
    {
        allocator_callbacks->free(allocator_callbacks->state, offsets);
        allocator_callbacks->free(allocator_callbacks->state, indices);
        allocator_callbacks->free(allocator_callbacks->state, mtx);
        return JMTX_RESULT_BAD_ALLOC;
    }
    memcpy(offsets, in->end_of_row_offsets, (in->base.rows) * sizeof(*offsets));
    memcpy(indices, in->indices, (in->n_entries) * sizeof(*indices));

    for (uint_fast32_t i = 0; i < in->n_entries; ++i)
    {
        values[i] = cimag(in->values[i]);
    }

    mtx->base.cols = in->base.cols;
    mtx->base.type = JMTXD_TYPE_CRS;
    mtx->base.rows = in->base.rows;
    mtx->base.allocator_callbacks = *allocator_callbacks;
    mtx->indices = indices;
    mtx->values = values;
    mtx->capacity = in->n_entries;
    mtx->n_entries = in->n_entries;
    mtx->end_of_row_offsets = offsets;
    *p_mtx = mtx;


    return JMTX_RESULT_SUCCESS;
}

/**
 * Creates a new CRS matrix with double precision from the real part of a complex CRS matrix with double precision.
 * Requires no memory allocation by reusing the memory of the initial matrix. Can not fail if the input matrix is valid.
 * @param in matrix which to convert (will be invalid if function succeeds)
 * @return converted matrix
 */
jmtxd_matrix_crs* jmtxd_matrix_crs_from_cdouble_real_inplace(jmtxz_matrix_crs* in)
{
    double* const values = (double*)in->values;

    for (uint_fast32_t i = 0; i < in->n_entries; ++i)
    {
        values[i] = creal(in->values[i]);
    }

    in->capacity = (in->capacity * sizeof(*in->values)) / sizeof(*values);
    in->base.type = JMTXD_TYPE_CRS;

    return (jmtxd_matrix_crs*)in;
}

/**
 * Creates a new CRS matrix with double precision from the imaginary part of a complex CRS matrix with double precision.
 * Requires no memory allocation by reusing the memory of the initial matrix. Can not fail if the input matrix is valid.
 * @param in matrix which to convert (will be invalid if function succeeds)
 * @return converted matrix
 */
jmtxd_matrix_crs* jmtxd_matrix_crs_from_cdouble_imag_inplace(jmtxz_matrix_crs* in)
{
    double* const values = (double*)in->values;

    for (uint_fast32_t i = 0; i < in->n_entries; ++i)
    {
        values[i] = cimag(in->values[i]);
    }

    in->capacity = (in->capacity * sizeof(*in->values)) / sizeof(*values);
    in->base.type = JMTX_TYPE_CRS;

    return (jmtxd_matrix_crs*)in;
}

/**
 * Creates a new CRS matrix with double precision from the real part of a complex CRS matrix with double precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxds_matrix_crs_from_cdouble_real(jmtxd_matrix_crs** p_mtx, const jmtxz_matrix_crs* in,
                                               const jmtx_allocator_callbacks* allocator_callbacks)
{
    if (!p_mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!in)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (allocator_callbacks && (!allocator_callbacks->free || !allocator_callbacks->alloc || !allocator_callbacks->realloc))
    {
        return JMTX_RESULT_BAD_PARAM;
    }
    if (in->base.type != JMTXZ_TYPE_CRS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }

    return jmtxd_matrix_crs_from_cdouble_real(p_mtx, in, allocator_callbacks);
}

/**
 * Creates a new CRS matrix with double precision from the imaginary part of a complex CRS matrix with double precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxds_matrix_crs_from_cdouble_imag(jmtxd_matrix_crs** p_mtx, const jmtxz_matrix_crs* in,
                                               const jmtx_allocator_callbacks* allocator_callbacks)
{
    if (!p_mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!in)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (allocator_callbacks && (!allocator_callbacks->free || !allocator_callbacks->alloc || !allocator_callbacks->realloc))
    {
        return JMTX_RESULT_BAD_PARAM;
    }
    if (in->base.type != JMTXZ_TYPE_CRS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }

    return jmtxd_matrix_crs_from_cdouble_imag(p_mtx, in, allocator_callbacks);
}

/**
 * Creates a new complex CRS matrix with double precision from a CRS matrix with double precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in_real matrix to use as the real component of the new matrix (if NULL real part is zero)
 * @param in_imag matrix to use as the real component of the new matrix (if NULL imaginary part is zero)
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxz_matrix_crs_from_double(jmtxz_matrix_crs** p_mtx, const jmtxd_matrix_crs* in_real,
                                         const jmtxd_matrix_crs* in_imag,
                                         const jmtx_allocator_callbacks* allocator_callbacks)
{
    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }

    uint32_t* offsets = NULL;
    uint32_t* indices = NULL;
    uint32_t r, c;

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
    jmtxz_matrix_crs* mtx = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*mtx));
    if (!mtx)
    {
        return JMTX_RESULT_BAD_ALLOC;
    }

    if (in_real && in_imag)
    {
        const uint32_t needed_capacity = in_real->n_entries + in_imag->n_entries;


        _Complex double* values = allocator_callbacks->alloc(allocator_callbacks->state, (needed_capacity) * sizeof(*values));
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

        offsets = allocator_callbacks->alloc(allocator_callbacks->state, (r) * sizeof(*offsets));
        if (!offsets)
        {
            allocator_callbacks->free(allocator_callbacks->state, offsets);
            allocator_callbacks->free(allocator_callbacks->state, indices);
            allocator_callbacks->free(allocator_callbacks->state, mtx);
            return JMTX_RESULT_BAD_ALLOC;
        }
        uint32_t i, j;
        for (i = 0, j = 0; i < r; ++i)
        {
            double* v_r, *v_i;
            uint32_t* i_r, *i_i;
            const uint32_t c_r = jmtxd_matrix_crs_get_row(in_real, i, &i_r, &v_r);
            const uint32_t c_i = jmtxd_matrix_crs_get_row(in_imag, i, &i_i, &v_i);
            uint32_t k_r, k_i;
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
                else //if (i_r[k_r] == i_i[k_i])
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
        mtx->end_of_row_offsets = offsets;

    }
    else if (in_imag)
    {
        _Complex double* values = allocator_callbacks->alloc(allocator_callbacks->state, (in_imag->n_entries) * sizeof(*values));
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

        offsets = allocator_callbacks->alloc(allocator_callbacks->state, (in_imag->base.rows) * sizeof(*offsets));
        if (!offsets)
        {
            allocator_callbacks->free(allocator_callbacks->state, offsets);
            allocator_callbacks->free(allocator_callbacks->state, indices);
            allocator_callbacks->free(allocator_callbacks->state, mtx);
            return JMTX_RESULT_BAD_ALLOC;
        }
        memcpy(offsets, in_imag->end_of_row_offsets, (in_imag->base.rows) * sizeof(*offsets));
        memcpy(indices, in_imag->indices, (in_imag->n_entries) * sizeof(*indices));

        for (uint_fast32_t i = 0; i < in_imag->n_entries; ++i)
        {
            values[i] = _Complex_I * in_imag->values[i];
        }

        mtx->indices = indices;
        mtx->values = values;
        mtx->capacity = in_imag->n_entries;
        mtx->n_entries = in_imag->n_entries;
        mtx->end_of_row_offsets = offsets;
    }
    else // if (in_real)
    {
        _Complex double* values = allocator_callbacks->alloc(allocator_callbacks->state, (in_real->n_entries) * sizeof(*values));
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

        offsets = allocator_callbacks->alloc(allocator_callbacks->state, (in_real->base.rows) * sizeof(*offsets));
        if (!offsets)
        {
            allocator_callbacks->free(allocator_callbacks->state, offsets);
            allocator_callbacks->free(allocator_callbacks->state, indices);
            allocator_callbacks->free(allocator_callbacks->state, mtx);
            return JMTX_RESULT_BAD_ALLOC;
        }
        memcpy(offsets, in_real->end_of_row_offsets, (in_real->base.rows) * sizeof(*offsets));
        memcpy(indices, in_real->indices, (in_real->n_entries) * sizeof(*indices));

        for (uint_fast32_t i = 0; i < in_real->n_entries; ++i)
        {
            values[i] = _Complex_I * in_real->values[i];
        }

        mtx->indices = indices;
        mtx->values = values;
        mtx->capacity = in_real->n_entries;
        mtx->n_entries = in_real->n_entries;
        mtx->end_of_row_offsets = offsets;

    }

    mtx->base.cols = c;
    mtx->base.type = JMTXZ_TYPE_CRS;
    mtx->base.rows = r;
    mtx->base.allocator_callbacks = *allocator_callbacks;
    *p_mtx = mtx;

    return JMTX_RESULT_SUCCESS;
}

/**
 * Creates a new complex CRS matrix with double precision from a CRS matrix with double precision as its real part.
 * Only one memory reallocation.
 * may be needed. Can not fail if the input matrix is valid.
 * @param in matrix which to convert
 * @return converted matrix, or NULL in case of allocation failure
 */
jmtxz_matrix_crs* jmtxz_matrix_crs_from_double_real_inplace(jmtxd_matrix_crs* in)
{
    _Complex double* vals = NULL;
    uint32_t capacity_for_complex = in->capacity * sizeof(*in->values) / sizeof(*vals);
    if (capacity_for_complex < in->n_entries)
    {
        capacity_for_complex = in->n_entries;
        vals = in->base.allocator_callbacks.realloc(in->base.allocator_callbacks.state, in->values, sizeof(*vals) * capacity_for_complex);
        if (!vals)
        {
            return NULL;
        }
        in->values = (double*)vals;
    }

    vals = (_Complex double*)in->values;
    in->capacity = capacity_for_complex;

    for (uint32_t i = 0; i < in->n_entries; ++i)
    {
        vals[in->n_entries - 1 - i] = (_Complex double)in->values[in->n_entries - 1 - i];
    }
    in->base.type = JMTXZ_TYPE_CRS;

    return (jmtxz_matrix_crs*)in;
}

/**
 * Creates a new CRS matrix with double precision from a CRS matrix with double precision as its imaginary part.
 * Only one memory reallocation.
 * may be needed. Can not fail if the input matrix is valid.
 * @param in matrix which to convert
 * @return converted matrix, or NULL in case of allocation failure
 */
jmtxz_matrix_crs* jmtxz_matrix_crs_from_double_imag_inplace(jmtxd_matrix_crs* in)
{
    _Complex double* vals = NULL;
    uint32_t capacity_for_complex = in->capacity * sizeof(*in->values) / sizeof(*vals);
    if (capacity_for_complex < in->n_entries)
    {
        capacity_for_complex = in->n_entries;
        vals = in->base.allocator_callbacks.realloc(in->base.allocator_callbacks.state, in->values, sizeof(*vals) * capacity_for_complex);
        if (!vals)
        {
            return NULL;
        }
        in->values = (double*)vals;
    }

    vals = (_Complex double*)in->values;
    in->capacity = capacity_for_complex;

    for (uint32_t i = 0; i < in->n_entries; ++i)
    {
        vals[in->n_entries - 1 - i] = (_Complex double)in->values[in->n_entries - 1 - i] * _Complex_I;
    }
    in->base.type = JMTXZ_TYPE_CRS;

    return (jmtxz_matrix_crs*)in;
}

/**
 * Creates a new complex CRS matrix with double precision from a CRS matrix with double precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in_real matrix to use as the real component of the new matrix (if NULL real part is zero)
 * @param in_imag matrix to use as the real component of the new matrix (if NULL imaginary part is zero)
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxzs_matrix_crs_from_double(jmtxz_matrix_crs** p_mtx, const jmtxd_matrix_crs* in_real,
                                          const jmtxd_matrix_crs* in_imag,
                                          const jmtx_allocator_callbacks* allocator_callbacks)
{
    if (!p_mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!in_real && !in_imag)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (allocator_callbacks && (!allocator_callbacks->free || !allocator_callbacks->alloc || !allocator_callbacks->realloc))
    {
        return JMTX_RESULT_BAD_PARAM;
    }
    if ((in_real && in_real->base.type != JMTXD_TYPE_CRS) ||
        (in_imag && in_imag->base.type != JMTXD_TYPE_CRS))
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if ((in_real && in_imag) && (in_real->base.rows != in_imag->base.rows || in_real->base.cols != in_imag->base.cols))
    {
        return JMTX_RESULT_BAD_MATRIX;
    }

    return jmtxz_matrix_crs_from_double(p_mtx, in_real, in_imag, allocator_callbacks);
}

/***********************************************************************************************************************
 *                                                                                                                     *
 *                                 COMPLEX FLOAT <->  COMPLEX DOUBLE                                                   *
 *                                                                                                                     *
 **********************************************************************************************************************/
/**
 * Creates a new complex CRS matrix with single precision from a complex CRS matrix with double precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxc_matrix_crs_from_cdouble(jmtxc_matrix_crs** p_mtx, const jmtxz_matrix_crs* in,
                                          const jmtx_allocator_callbacks* allocator_callbacks)
{
    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }

    uint32_t* offsets = NULL;
    uint32_t* indices = NULL;

    jmtxc_matrix_crs* mtx = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*mtx));
    if (!mtx)
    {
        return JMTX_RESULT_BAD_ALLOC;
    }

    _Complex float* values = allocator_callbacks->alloc(allocator_callbacks->state, (in->n_entries) * sizeof(*values));
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

    offsets = allocator_callbacks->alloc(allocator_callbacks->state, (in->base.rows) * sizeof(*offsets));
    if (!offsets)
    {
        allocator_callbacks->free(allocator_callbacks->state, offsets);
        allocator_callbacks->free(allocator_callbacks->state, indices);
        allocator_callbacks->free(allocator_callbacks->state, mtx);
        return JMTX_RESULT_BAD_ALLOC;
    }
    memcpy(offsets, in->end_of_row_offsets, (in->base.rows) * sizeof(*offsets));
    memcpy(indices, in->indices, (in->n_entries) * sizeof(*indices));

    for (uint_fast32_t i = 0; i < in->n_entries; ++i)
    {
        values[i] = (_Complex float)in->values[i];
    }

    mtx->base.cols = in->base.cols;
    mtx->base.type = JMTXC_TYPE_CRS;
    mtx->base.rows = in->base.rows;
    mtx->base.allocator_callbacks = *allocator_callbacks;
    mtx->indices = indices;
    mtx->values = values;
    mtx->capacity = in->n_entries;
    mtx->n_entries = in->n_entries;
    mtx->end_of_row_offsets = offsets;
    *p_mtx = mtx;


    return JMTX_RESULT_SUCCESS;
}
/**
 * Creates a new complex CRS matrix with single precision from a complex CRS matrix with double precision. Requires no memory
 * allocation by reusing the memory of the initial matrix. Can not fail if the input matrix is valid.
 * @param in matrix which to convert (will be invalid if function succeeds)
 * @return converted matrix
 */
jmtxc_matrix_crs* jmtxc_matrix_crs_from_cdouble_inplace(jmtxz_matrix_crs* in)
{
    _Complex float* const values = (_Complex float*)in->values;

    for (uint_fast32_t i = 0; i < in->n_entries; ++i)
    {
        values[i] = (_Complex float)in->values[i];
    }

    in->capacity = (in->capacity * sizeof(*in->values)) / sizeof(*values);
    in->base.type = JMTXC_TYPE_CRS;

    return (jmtxc_matrix_crs*)in;
}

/**
 * Creates a new complex CRS matrix with single precision from a complex CRS matrix with double precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxcs_matrix_crs_from_cdouble(jmtxc_matrix_crs** p_mtx, const jmtxz_matrix_crs* in,
                                           const jmtx_allocator_callbacks* allocator_callbacks)
{
    if (!p_mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!in)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (allocator_callbacks && (!allocator_callbacks->free || !allocator_callbacks->alloc || !allocator_callbacks->realloc))
    {
        return JMTX_RESULT_BAD_PARAM;
    }
    if (in->base.type != JMTXZ_TYPE_CRS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }

    return jmtxc_matrix_crs_from_cdouble(p_mtx, in, allocator_callbacks);
}

/**
 * Creates a new complex CRS matrix with double precision from a complex CRS matrix with single precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxz_matrix_crs_from_cfloat(jmtxz_matrix_crs** p_mtx, const jmtxc_matrix_crs* in,
                                         const jmtx_allocator_callbacks* allocator_callbacks)
{
    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }

    uint32_t* offsets = NULL;
    uint32_t* indices = NULL;

    jmtxz_matrix_crs* mtx = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*mtx));
    if (!mtx)
    {
        return JMTX_RESULT_BAD_ALLOC;
    }

    _Complex double* values = allocator_callbacks->alloc(allocator_callbacks->state, (in->n_entries) * sizeof(*values));
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

    offsets = allocator_callbacks->alloc(allocator_callbacks->state, (in->base.rows) * sizeof(*offsets));
    if (!offsets)
    {
        allocator_callbacks->free(allocator_callbacks->state, offsets);
        allocator_callbacks->free(allocator_callbacks->state, indices);
        allocator_callbacks->free(allocator_callbacks->state, mtx);
        return JMTX_RESULT_BAD_ALLOC;
    }
    memcpy(offsets, in->end_of_row_offsets, (in->base.rows) * sizeof(*offsets));
    memcpy(indices, in->indices, (in->n_entries) * sizeof(*indices));

    for (uint_fast32_t i = 0; i < in->n_entries; ++i)
    {
        values[i] = (_Complex double)in->values[i];
    }

    mtx->base.cols = in->base.cols;
    mtx->base.type = JMTXZ_TYPE_CRS;
    mtx->base.rows = in->base.rows;
    mtx->base.allocator_callbacks = *allocator_callbacks;
    mtx->indices = indices;
    mtx->values = values;
    mtx->capacity = in->n_entries;
    mtx->n_entries = in->n_entries;
    mtx->end_of_row_offsets = offsets;
    *p_mtx = mtx;

    return JMTX_RESULT_SUCCESS;
}

/**
 * Creates a new complex CRS matrix with double precision from a complex CRS matrix with single precision. Only one memory reallocation
 * may be needed. Can not fail if the input matrix is valid.
 * @param in matrix which to convert
 * @return converted matrix, or NULL in case of allocation failure
 */
jmtxz_matrix_crs* jmtxz_matrix_crs_from_cfloat_inplace(jmtxc_matrix_crs* in)
{
    _Complex double* vals = NULL;
    uint32_t capacity_for_doubles = in->capacity * sizeof(*in->values) / sizeof(*vals);
    if (capacity_for_doubles < in->n_entries)
    {
        capacity_for_doubles = in->n_entries;
        vals = in->base.allocator_callbacks.realloc(in->base.allocator_callbacks.state, in->values, sizeof(*vals) * capacity_for_doubles);
        if (!vals)
        {
            return NULL;
        }
        in->values = (_Complex float*)vals;
    }

    vals = (_Complex double*)in->values;
    in->capacity = capacity_for_doubles;

    for (uint32_t i = 0; i < in->n_entries; ++i)
    {
        vals[in->n_entries - 1 - i] = (_Complex double)in->values[in->n_entries - 1 - i];
    }
    in->base.type = JMTXD_TYPE_CRS;

    return (jmtxz_matrix_crs*)in;
}

/**
 * Creates a new complex CRS matrix with double precision from a complex CRS matrix with single precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxzs_matrix_crs_from_cfloat(jmtxz_matrix_crs** p_mtx, const jmtxc_matrix_crs* in,
                                          const jmtx_allocator_callbacks* allocator_callbacks)
{
    if (!p_mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!in)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (allocator_callbacks && (!allocator_callbacks->free || !allocator_callbacks->alloc || !allocator_callbacks->realloc))
    {
        return JMTX_RESULT_BAD_PARAM;
    }
    if (in->base.type != JMTXC_TYPE_CRS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }

    return jmtxz_matrix_crs_from_cfloat(p_mtx, in, allocator_callbacks);
}

