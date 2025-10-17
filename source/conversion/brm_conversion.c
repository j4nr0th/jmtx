//
// Created by jan on 1.12.2023.
//

#include <assert.h>
#include <complex.h>

#include "../double/matrices/band_row_major.h"
#include "../float/matrices/band_row_major.h"
#include "../matrix_base.h"
#ifndef _MSC_BUILD
#    include "../cdouble/matrices/band_row_major.h"
#    include "../cfloat/matrices/band_row_major.h"
#endif
#include "brm_conversion.h"

static inline JMTX_FAST_INT_T brm_row_offset(const jmtxf_matrix_brm mtx[JMTX_ARRAY_ATTRIB(const static 1)],
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

static inline JMTX_FAST_INT_T brm_row_offsetd(const JMTX_NAME_TYPED(matrix_brm) mtx[JMTX_ARRAY_ATTRIB(const static 1)],
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

#ifndef _MSC_BUILD

static inline JMTX_FAST_INT_T brm_row_offsetc(const jmtxc_matrix_brm mtx[JMTX_ARRAY_ATTRIB(const static 1)],
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

static inline JMTX_FAST_INT_T brm_row_offsetz(const jmtxz_matrix_brm mtx[JMTX_ARRAY_ATTRIB(const static 1)],
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

#endif

/***********************************************************************************************************************
 *                                                                                                                     *
 *                                          FLOAT <-> DOUBLE                                                           *
 *                                                                                                                     *
 **********************************************************************************************************************/

/**
 * Creates a new BRM matrix with single precision from a BRM matrix with JMTX_SCALAR_T precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxf_matrix_brm_from_double(jmtxf_matrix_brm **p_mtx, const JMTX_NAME_TYPED(matrix_brm) * in,
                                         const jmtx_allocator_callbacks *allocator_callbacks)
{
    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }

    jmtxf_matrix_brm *mtx = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*mtx));
    if (!mtx)
    {
        return JMTX_RESULT_BAD_ALLOC;
    }

    const JMTX_FAST_INT_T count = brm_row_offsetd(in, in->base.rows);
    float *values = allocator_callbacks->alloc(allocator_callbacks->state, count * sizeof(*values));
    if (!values)
    {
        allocator_callbacks->free(allocator_callbacks->state, mtx);
        return JMTX_RESULT_BAD_ALLOC;
    }

    for (JMTX_FAST_INT_T i = 0; i < count; ++i)
    {
        values[i] = (float)in->values[i];
    }

    mtx->base.cols = in->base.cols;
    mtx->base.type = JMTX_TYPE_BRM;
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
jmtxf_matrix_brm *jmtxf_matrix_brm_from_double_inplace(JMTX_NAME_TYPED(matrix_brm) * in)
{
    const JMTX_FAST_INT_T count = brm_row_offsetd(in, in->base.rows);
    float *const values = (float *)in->values;

    for (JMTX_FAST_INT_T i = 0; i < count; ++i)
    {
        values[i] = (float)in->values[i];
    }

    float *const new_ptr =
        in->base.allocator_callbacks.realloc(in->base.allocator_callbacks.state, values, sizeof(*new_ptr) * count);
    if (new_ptr)
    {
        in->values = (JMTX_SCALAR_T *)new_ptr;
    }

    in->base.type = JMTX_TYPE_BRM;
    return (jmtxf_matrix_brm *)in;
}

/**
 * Creates a new BRM matrix with single precision from a BRM matrix with JMTX_SCALAR_T precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxs_matrix_brm_from_double(jmtxf_matrix_brm **p_mtx, const JMTX_NAME_TYPED(matrix_brm) * in,
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
    if (in->base.type != JMTXD_TYPE_BRM)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }

    return jmtxf_matrix_brm_from_double(p_mtx, in, allocator_callbacks);
}

/**
 * Creates a new BRM matrix with JMTX_SCALAR_T precision from a BRM matrix with single precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result JMTX_NAME_TYPED(matrix_brm_from_float(JMTX_NAME_TYPED(matrix_brm) **p_mtx, const jmtxf_matrix_brm *in,
                                        const jmtx_allocator_callbacks *allocator_callbacks)
{
    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }

    JMTX_NAME_TYPED(matrix_brm) *mtx = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*mtx));
    if (!mtx)
    {
        return JMTX_RESULT_BAD_ALLOC;
    }

    const JMTX_FAST_INT_T count = brm_row_offset(in, in->base.rows);
    JMTX_SCALAR_T *values = allocator_callbacks->alloc(allocator_callbacks->state, count * sizeof(*values));
    if (!values)
    {
        allocator_callbacks->free(allocator_callbacks->state, mtx);
        return JMTX_RESULT_BAD_ALLOC;
    }

    for (JMTX_FAST_INT_T i = 0; i < count; ++i)
    {
        values[i] = (float)in->values[i];
    }

    mtx->base.cols = in->base.cols;
    mtx->base.type = JMTXD_TYPE_BRM;
    mtx->base.rows = in->base.rows;
    mtx->base.allocator_callbacks = *allocator_callbacks;
    mtx->values = values;
    mtx->lower_bandwidth = in->lower_bandwidth;
    mtx->upper_bandwidth = in->upper_bandwidth;

    *p_mtx = mtx;

    return JMTX_RESULT_SUCCESS;
}

/**
 * Creates a new BRM matrix with JMTX_SCALAR_T precision from a BRM matrix with single precision. Only one memory reallocation
 * may be needed. Can not fail if the input matrix is valid.
 * @param in matrix which to convert
 * @return converted matrix, or NULL in case of allocation failure
 */
JMTX_NAME_TYPED(matrix_brm) *JMTX_NAME_TYPED(matrix_brm_from_float_inplace(jmtxf_matrix_brm *in)
{
    const JMTX_FAST_INT_T count = brm_row_offset(in, in->base.rows);
    JMTX_SCALAR_T *const values =
        in->base.allocator_callbacks.realloc(in->base.allocator_callbacks.state, in->values, sizeof(*values) * count);
    if (!values)
    {
        return NULL;
    }
    in->values = (float *)values;

    for (JMTX_FAST_INT_T i = 0; i < count; ++i)
    {
        values[count - 1 - i] = (JMTX_SCALAR_T)in->values[count - 1 - i];
    }

    in->base.type = JMTXD_TYPE_BRM;
    return (JMTX_NAME_TYPED(matrix_brm) *)in;
}

/**
 * Creates a new BRM matrix with JMTX_SCALAR_T precision from a BRM matrix with single precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result .,.,.,.,..matrix_brm_from_float(JMTX_NAME_TYPED(matrix_brm) **p_mtx, const jmtxf_matrix_brm *in,
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
    if (in->base.type != JMTX_TYPE_BRM)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }

    return JMTX_NAME_TYPED(matrix_brm_from_float(p_mtx, in, allocator_callbacks);
}

#ifndef _MSC_BUILD
/***********************************************************************************************************************
 *                                                                                                                     *
 *                                          FLOAT <-> COMPLEX FLOAT                                                    *
 *                                                                                                                     *
 **********************************************************************************************************************/

/**
 * Creates a new BRM matrix with single precision from the real part of a complex BRM matrix with single precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxf_matrix_brm_from_cfloat_real(jmtxf_matrix_brm **p_mtx, const jmtxc_matrix_brm *in,
                                              const jmtx_allocator_callbacks *allocator_callbacks)
{
    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }

    jmtxf_matrix_brm *mtx = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*mtx));
    if (!mtx)
    {
        return JMTX_RESULT_BAD_ALLOC;
    }

    const JMTX_FAST_INT_T count = brm_row_offsetc(in, in->base.rows);
    float *values = allocator_callbacks->alloc(allocator_callbacks->state, count * sizeof(*values));
    if (!values)
    {
        allocator_callbacks->free(allocator_callbacks->state, mtx);
        return JMTX_RESULT_BAD_ALLOC;
    }

    for (JMTX_FAST_INT_T i = 0; i < count; ++i)
    {
        values[i] = crealf(in->values[i]);
    }

    mtx->base.cols = in->base.cols;
    mtx->base.type = JMTX_TYPE_BRM;
    mtx->base.rows = in->base.rows;
    mtx->base.allocator_callbacks = *allocator_callbacks;
    mtx->values = values;
    mtx->lower_bandwidth = in->lower_bandwidth;
    mtx->upper_bandwidth = in->upper_bandwidth;

    *p_mtx = mtx;

    return JMTX_RESULT_SUCCESS;
}
/**
 * Creates a new BRM matrix with single precision from the imaginary part of a complex BRM matrix with single precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxf_matrix_brm_from_cfloat_imag(jmtxf_matrix_brm **p_mtx, const jmtxc_matrix_brm *in,
                                              const jmtx_allocator_callbacks *allocator_callbacks)
{
    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }

    jmtxf_matrix_brm *mtx = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*mtx));
    if (!mtx)
    {
        return JMTX_RESULT_BAD_ALLOC;
    }

    const JMTX_FAST_INT_T count = brm_row_offsetc(in, in->base.rows);
    float *values = allocator_callbacks->alloc(allocator_callbacks->state, count * sizeof(*values));
    if (!values)
    {
        allocator_callbacks->free(allocator_callbacks->state, mtx);
        return JMTX_RESULT_BAD_ALLOC;
    }

    for (JMTX_FAST_INT_T i = 0; i < count; ++i)
    {
        values[i] = cimagf(in->values[i]);
    }

    mtx->base.cols = in->base.cols;
    mtx->base.type = JMTX_TYPE_BRM;
    mtx->base.rows = in->base.rows;
    mtx->base.allocator_callbacks = *allocator_callbacks;
    mtx->values = values;
    mtx->lower_bandwidth = in->lower_bandwidth;
    mtx->upper_bandwidth = in->upper_bandwidth;

    *p_mtx = mtx;

    return JMTX_RESULT_SUCCESS;
}

/**
 * Creates a new BRM matrix with single precision from the real part of a complex BRM matrix with single precision.
 * Requires no memory allocation by reusing the memory of the initial matrix. Can not fail if the input matrix is valid.
 * @param in matrix which to convert (will be invalid if function succeeds)
 * @return converted matrix
 */
jmtxf_matrix_brm *jmtxf_matrix_brm_from_cfloat_real_inplace(jmtxc_matrix_brm *in)
{
    const JMTX_FAST_INT_T count = brm_row_offsetc(in, in->base.rows);
    float *const values = (float *)in->values;

    for (JMTX_FAST_INT_T i = 0; i < count; ++i)
    {
        values[i] = crealf(in->values[i]);
    }

    float *const new_ptr =
        in->base.allocator_callbacks.realloc(in->base.allocator_callbacks.state, values, sizeof(*new_ptr) * count);
    if (new_ptr)
    {
        in->values = (_Complex float *)new_ptr;
    }

    in->base.type = JMTX_TYPE_BRM;
    return (jmtxf_matrix_brm *)in;
}

/**
 * Creates a new BRM matrix with single precision from the imaginary part of a complex BRM matrix with single precision.
 * Requires no memory allocation by reusing the memory of the initial matrix. Can not fail if the input matrix is valid.
 * @param in matrix which to convert (will be invalid if function succeeds)
 * @return converted matrix
 */
jmtxf_matrix_brm *jmtxf_matrix_brm_from_cfloat_imag_inplace(jmtxc_matrix_brm *in)
{
    const JMTX_FAST_INT_T count = brm_row_offsetc(in, in->base.rows);
    float *const values = (float *)in->values;

    for (JMTX_FAST_INT_T i = 0; i < count; ++i)
    {
        values[i] = cimagf(in->values[i]);
    }

    float *const new_ptr =
        in->base.allocator_callbacks.realloc(in->base.allocator_callbacks.state, values, sizeof(*new_ptr) * count);
    if (new_ptr)
    {
        in->values = (_Complex float *)new_ptr;
    }

    in->base.type = JMTX_TYPE_BRM;
    return (jmtxf_matrix_brm *)in;
}

/**
 * Creates a new BRM matrix with single precision from the real part of a complex BRM matrix with single precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxs_matrix_brm_from_cfloat_real(jmtxf_matrix_brm **p_mtx, const jmtxc_matrix_brm *in,
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
    if (in->base.type != JMTXC_TYPE_BRM)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }

    return jmtxf_matrix_brm_from_cfloat_real(p_mtx, in, allocator_callbacks);
}

/**
 * Creates a new BRM matrix with single precision from the imaginary part of a complex BRM matrix with single precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxs_matrix_brm_from_cfloat_imag(jmtxf_matrix_brm **p_mtx, const jmtxc_matrix_brm *in,
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
    if (in->base.type != JMTXC_TYPE_BRM)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }

    return jmtxf_matrix_brm_from_cfloat_imag(p_mtx, in, allocator_callbacks);
}

/**
 * Creates a new complex BRM matrix with single precision from a BRM matrix with single precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in_real matrix to use as the real component of the new matrix (if NULL real part is zero)
 * @param in_imag matrix to use as the real component of the new matrix (if NULL imaginary part is zero)
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxc_matrix_brm_from_float(jmtxc_matrix_brm **p_mtx, const jmtxf_matrix_brm *in_real,
                                        const jmtxf_matrix_brm *in_imag,
                                        const jmtx_allocator_callbacks *allocator_callbacks)
{
    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }
    assert(in_imag || in_real);
    assert((!in_imag || !in_real) ||
           ((in_imag->lower_bandwidth == in_real->lower_bandwidth) &&
            (in_imag->upper_bandwidth == in_real->upper_bandwidth) && (in_imag->base.rows == in_real->base.rows) &&
            (in_imag->base.cols == in_real->base.cols)));

    jmtxc_matrix_brm *mtx = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*mtx));
    if (!mtx)
    {
        return JMTX_RESULT_BAD_ALLOC;
    }

    JMTX_INDEX_T r, c, lbw, ubw;
    if (in_real)
    {
        r = in_real->base.rows;
        c = in_real->base.cols;
        lbw = in_real->lower_bandwidth;
        ubw = in_real->upper_bandwidth;
    }
    else
    {
        r = in_imag->base.rows;
        c = in_imag->base.cols;
        lbw = in_imag->lower_bandwidth;
        ubw = in_imag->upper_bandwidth;
    }

    const JMTX_FAST_INT_T count = in_real ? brm_row_offset(in_real, r) : brm_row_offset(in_imag, c);
    _Complex float *values = allocator_callbacks->alloc(allocator_callbacks->state, count * sizeof(*values));
    if (!values)
    {
        allocator_callbacks->free(allocator_callbacks->state, mtx);
        return JMTX_RESULT_BAD_ALLOC;
    }

    if (in_real && in_imag)
    {
        for (JMTX_FAST_INT_T i = 0; i < count; ++i)
        {
            values[i] = in_real->values[i] + in_imag->values[i] * _Complex_I;
        }
    }
    else if (in_real)
    {
        for (JMTX_FAST_INT_T i = 0; i < count; ++i)
        {
            values[i] = (_Complex float)in_real->values[i];
        }
    }
    else // if (in_imag)
    {
        for (JMTX_FAST_INT_T i = 0; i < count; ++i)
        {
            values[i] = _Complex_I * in_imag->values[i];
        }
    }

    mtx->base.cols = c;
    mtx->base.type = JMTXC_TYPE_BRM;
    mtx->base.rows = r;
    mtx->base.allocator_callbacks = *allocator_callbacks;
    mtx->values = values;
    mtx->lower_bandwidth = lbw;
    mtx->upper_bandwidth = ubw;

    *p_mtx = mtx;

    return JMTX_RESULT_SUCCESS;
}

/**
 * Creates a new complex BRM matrix with single precision from a BRM matrix with single precision as its real part.
 * Only one memory reallocation.
 * may be needed. Can not fail if the input matrix is valid.
 * @param in matrix which to convert
 * @return converted matrix, or NULL in case of allocation failure
 */
jmtxc_matrix_brm *jmtxc_matrix_brm_from_float_real_inplace(jmtxf_matrix_brm *in)
{
    const JMTX_FAST_INT_T count = brm_row_offset(in, in->base.rows);
    _Complex float *const values =
        in->base.allocator_callbacks.realloc(in->base.allocator_callbacks.state, in->values, sizeof(*values) * count);
    if (!values)
    {
        return NULL;
    }
    in->values = (float *)values;

    for (JMTX_FAST_INT_T i = 0; i < count; ++i)
    {
        values[count - 1 - i] = (_Complex float)in->values[count - 1 - i];
    }

    in->base.type = JMTXC_TYPE_BRM;
    return (jmtxc_matrix_brm *)in;
}

/**
 * Creates a new BRM matrix with single precision from a BRM matrix with single precision as its imaginary part.
 * Only one memory reallocation.
 * may be needed. Can not fail if the input matrix is valid.
 * @param in matrix which to convert
 * @return converted matrix, or NULL in case of allocation failure
 */
jmtxc_matrix_brm *jmtxc_matrix_brm_from_float_imag_inplace(jmtxf_matrix_brm *in)
{
    const JMTX_FAST_INT_T count = brm_row_offset(in, in->base.rows);
    _Complex float *const values =
        in->base.allocator_callbacks.realloc(in->base.allocator_callbacks.state, in->values, sizeof(*values) * count);
    if (!values)
    {
        return NULL;
    }
    in->values = (float *)values;

    for (JMTX_FAST_INT_T i = 0; i < count; ++i)
    {
        values[count - 1 - i] = _Complex_I * in->values[count - 1 - i];
    }

    in->base.type = JMTXC_TYPE_BRM;
    return (jmtxc_matrix_brm *)in;
}

/**
 * Creates a new complex BRM matrix with single precision from a BRM matrix with single precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in_real matrix to use as the real component of the new matrix (if NULL real part is zero)
 * @param in_imag matrix to use as the real component of the new matrix (if NULL imaginary part is zero)
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxcs_matrix_brm_from_float(jmtxc_matrix_brm **p_mtx, const jmtxf_matrix_brm *in_real,
                                         const jmtxf_matrix_brm *in_imag,
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
    if ((in_real && in_real->base.type != JMTX_TYPE_BRM) || (in_imag && in_imag->base.type != JMTX_TYPE_BRM))
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if ((in_real && in_imag) &&
        (in_real->base.rows != in_imag->base.rows || in_real->base.cols != in_imag->base.cols ||
         in_real->lower_bandwidth != in_imag->lower_bandwidth || in_real->upper_bandwidth != in_imag->upper_bandwidth))
    {
        return JMTX_RESULT_BAD_MATRIX;
    }

    return jmtxc_matrix_brm_from_float(p_mtx, in_real, in_imag, allocator_callbacks);
}

/***********************************************************************************************************************
 *                                                                                                                     *
 *                                          DOUBLE <-> COMPLEX DOUBLE                                                  *
 *                                                                                                                     *
 **********************************************************************************************************************/
/**
 * Creates a new BRM matrix with JMTX_SCALAR_T precision from the real part of a complex BRM matrix with JMTX_SCALAR_T precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result JMTX_NAME_TYPED(matrix_brm_from_cdouble_real(JMTX_NAME_TYPED(matrix_brm) **p_mtx, const jmtxz_matrix_brm *in,
                                               const jmtx_allocator_callbacks *allocator_callbacks)
{
    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }

    JMTX_NAME_TYPED(matrix_brm) *mtx = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*mtx));
    if (!mtx)
    {
        return JMTX_RESULT_BAD_ALLOC;
    }

    const JMTX_FAST_INT_T count = brm_row_offsetz(in, in->base.rows);
    JMTX_SCALAR_T *values = allocator_callbacks->alloc(allocator_callbacks->state, count * sizeof(*values));
    if (!values)
    {
        allocator_callbacks->free(allocator_callbacks->state, mtx);
        return JMTX_RESULT_BAD_ALLOC;
    }

    for (JMTX_FAST_INT_T i = 0; i < count; ++i)
    {
        values[i] = crealf(in->values[i]);
    }

    mtx->base.cols = in->base.cols;
    mtx->base.type = JMTXD_TYPE_BRM;
    mtx->base.rows = in->base.rows;
    mtx->base.allocator_callbacks = *allocator_callbacks;
    mtx->values = values;
    mtx->lower_bandwidth = in->lower_bandwidth;
    mtx->upper_bandwidth = in->upper_bandwidth;

    *p_mtx = mtx;

    return JMTX_RESULT_SUCCESS;
}
/**
 * Creates a new BRM matrix with JMTX_SCALAR_T precision from the imaginary part of a complex BRM matrix with JMTX_SCALAR_T precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result JMTX_NAME_TYPED(matrix_brm_from_cdouble_imag(JMTX_NAME_TYPED(matrix_brm) **p_mtx, const jmtxz_matrix_brm *in,
                                               const jmtx_allocator_callbacks *allocator_callbacks)
{
    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }

    JMTX_NAME_TYPED(matrix_brm) *mtx = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*mtx));
    if (!mtx)
    {
        return JMTX_RESULT_BAD_ALLOC;
    }

    const JMTX_FAST_INT_T count = brm_row_offsetz(in, in->base.rows);
    JMTX_SCALAR_T *values = allocator_callbacks->alloc(allocator_callbacks->state, count * sizeof(*values));
    if (!values)
    {
        allocator_callbacks->free(allocator_callbacks->state, mtx);
        return JMTX_RESULT_BAD_ALLOC;
    }

    for (JMTX_FAST_INT_T i = 0; i < count; ++i)
    {
        values[i] = cimagf(in->values[i]);
    }

    mtx->base.cols = in->base.cols;
    mtx->base.type = JMTXD_TYPE_BRM;
    mtx->base.rows = in->base.rows;
    mtx->base.allocator_callbacks = *allocator_callbacks;
    mtx->values = values;
    mtx->lower_bandwidth = in->lower_bandwidth;
    mtx->upper_bandwidth = in->upper_bandwidth;

    *p_mtx = mtx;

    return JMTX_RESULT_SUCCESS;
}

/**
 * Creates a new BRM matrix with JMTX_SCALAR_T precision from the real part of a complex BRM matrix with JMTX_SCALAR_T precision.
 * Requires no memory allocation by reusing the memory of the initial matrix. Can not fail if the input matrix is valid.
 * @param in matrix which to convert (will be invalid if function succeeds)
 * @return converted matrix
 */
JMTX_NAME_TYPED(matrix_brm) *JMTX_NAME_TYPED(matrix_brm_from_cdouble_real_inplace(jmtxz_matrix_brm *in)
{
    const JMTX_FAST_INT_T count = brm_row_offsetz(in, in->base.rows);
    JMTX_SCALAR_T *const values = (JMTX_SCALAR_T *)in->values;

    for (JMTX_FAST_INT_T i = 0; i < count; ++i)
    {
        values[i] = creal(in->values[i]);
    }

    JMTX_SCALAR_T *const new_ptr =
        in->base.allocator_callbacks.realloc(in->base.allocator_callbacks.state, values, sizeof(*new_ptr) * count);
    if (new_ptr)
    {
        in->values = (_Complex JMTX_SCALAR_T *)new_ptr;
    }

    in->base.type = JMTXD_TYPE_BRM;
    return (JMTX_NAME_TYPED(matrix_brm) *)in;
}

/**
 * Creates a new BRM matrix with JMTX_SCALAR_T precision from the imaginary part of a complex BRM matrix with JMTX_SCALAR_T precision.
 * Requires no memory allocation by reusing the memory of the initial matrix. Can not fail if the input matrix is valid.
 * @param in matrix which to convert (will be invalid if function succeeds)
 * @return converted matrix
 */
JMTX_NAME_TYPED(matrix_brm) *JMTX_NAME_TYPED(matrix_brm_from_cdouble_imag_inplace(jmtxz_matrix_brm *in)
{
    const JMTX_FAST_INT_T count = brm_row_offsetz(in, in->base.rows);
    JMTX_SCALAR_T *const values = (JMTX_SCALAR_T *)in->values;

    for (JMTX_FAST_INT_T i = 0; i < count; ++i)
    {
        values[i] = cimag(in->values[i]);
    }

    JMTX_SCALAR_T *const new_ptr =
        in->base.allocator_callbacks.realloc(in->base.allocator_callbacks.state, values, sizeof(*new_ptr) * count);
    if (new_ptr)
    {
        in->values = (_Complex JMTX_SCALAR_T *)new_ptr;
    }

    in->base.type = JMTXD_TYPE_BRM;
    return (JMTX_NAME_TYPED(matrix_brm) *)in;
}

/**
 * Creates a new BRM matrix with JMTX_SCALAR_T precision from the real part of a complex BRM matrix with JMTX_SCALAR_T precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result .,.,.,.,..matrix_brm_from_cdouble_real(JMTX_NAME_TYPED(matrix_brm) **p_mtx, const jmtxz_matrix_brm *in,
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
    if (in->base.type != JMTXZ_TYPE_BRM)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }

    return JMTX_NAME_TYPED(matrix_brm_from_cdouble_real(p_mtx, in, allocator_callbacks);
}

/**
 * Creates a new BRM matrix with JMTX_SCALAR_T precision from the imaginary part of a complex BRM matrix with JMTX_SCALAR_T precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result .,.,.,.,..matrix_brm_from_cdouble_imag(JMTX_NAME_TYPED(matrix_brm) **p_mtx, const jmtxz_matrix_brm *in,
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
    if (in->base.type != JMTXZ_TYPE_BRM)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }

    return JMTX_NAME_TYPED(matrix_brm_from_cdouble_imag(p_mtx, in, allocator_callbacks);
}

/**
 * Creates a new complex BRM matrix with JMTX_SCALAR_T precision from a BRM matrix with JMTX_SCALAR_T precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in_real matrix to use as the real component of the new matrix (if NULL real part is zero)
 * @param in_imag matrix to use as the real component of the new matrix (if NULL imaginary part is zero)
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxz_matrix_brm_from_double(jmtxz_matrix_brm **p_mtx, const JMTX_NAME_TYPED(matrix_brm) *in_real,
                                         const JMTX_NAME_TYPED(matrix_brm) *in_imag,
                                         const jmtx_allocator_callbacks *allocator_callbacks)
{
    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }
    assert(in_imag || in_real);
    assert((!in_imag || !in_real) ||
           ((in_imag->lower_bandwidth == in_real->lower_bandwidth) &&
            (in_imag->upper_bandwidth == in_real->upper_bandwidth) && (in_imag->base.rows == in_real->base.rows) &&
            (in_imag->base.cols == in_real->base.cols)));

    jmtxz_matrix_brm *mtx = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*mtx));
    if (!mtx)
    {
        return JMTX_RESULT_BAD_ALLOC;
    }

    JMTX_INDEX_T r, c, lbw, ubw;
    if (in_real)
    {
        r = in_real->base.rows;
        c = in_real->base.cols;
        lbw = in_real->lower_bandwidth;
        ubw = in_real->upper_bandwidth;
    }
    else
    {
        r = in_imag->base.rows;
        c = in_imag->base.cols;
        lbw = in_imag->lower_bandwidth;
        ubw = in_imag->upper_bandwidth;
    }

    const JMTX_FAST_INT_T count = in_real ? brm_row_offsetd(in_real, r) : brm_row_offsetd(in_imag, c);
    _Complex JMTX_SCALAR_T *values = allocator_callbacks->alloc(allocator_callbacks->state, count * sizeof(*values));
    if (!values)
    {
        allocator_callbacks->free(allocator_callbacks->state, mtx);
        return JMTX_RESULT_BAD_ALLOC;
    }

    if (in_real && in_imag)
    {
        for (JMTX_FAST_INT_T i = 0; i < count; ++i)
        {
            values[i] = in_real->values[i] + in_imag->values[i] * _Complex_I;
        }
    }
    else if (in_real)
    {
        for (JMTX_FAST_INT_T i = 0; i < count; ++i)
        {
            values[i] = (_Complex double)in_real->values[i];
        }
    }
    else // if (in_imag)
    {
        for (JMTX_FAST_INT_T i = 0; i < count; ++i)
        {
            values[i] = _Complex_I * in_imag->values[i];
        }
    }

    mtx->base.cols = c;
    mtx->base.type = JMTXZ_TYPE_BRM;
    mtx->base.rows = r;
    mtx->base.allocator_callbacks = *allocator_callbacks;
    mtx->values = values;
    mtx->lower_bandwidth = lbw;
    mtx->upper_bandwidth = ubw;

    *p_mtx = mtx;

    return JMTX_RESULT_SUCCESS;
}

/**
 * Creates a new complex BRM matrix with JMTX_SCALAR_T precision from a BRM matrix with JMTX_SCALAR_T precision as its real part.
 * Only one memory reallocation.
 * may be needed. Can not fail if the input matrix is valid.
 * @param in matrix which to convert
 * @return converted matrix, or NULL in case of allocation failure
 */
jmtxz_matrix_brm *jmtxz_matrix_brm_from_double_real_inplace(JMTX_NAME_TYPED(matrix_brm) *in)
{
    const JMTX_FAST_INT_T count = brm_row_offsetd(in, in->base.rows);
    _Complex JMTX_SCALAR_T *const values =
        in->base.allocator_callbacks.realloc(in->base.allocator_callbacks.state, in->values, sizeof(*values) * count);
    if (!values)
    {
        return NULL;
    }
    in->values = (JMTX_SCALAR_T *)values;

    for (JMTX_FAST_INT_T i = 0; i < count; ++i)
    {
        values[count - 1 - i] = (_Complex double)in->values[count - 1 - i];
    }

    in->base.type = JMTXZ_TYPE_BRM;
    return (jmtxz_matrix_brm *)in;
}

/**
 * Creates a new BRM matrix with JMTX_SCALAR_T precision from a BRM matrix with JMTX_SCALAR_T precision as its imaginary part.
 * Only one memory reallocation.
 * may be needed. Can not fail if the input matrix is valid.
 * @param in matrix which to convert
 * @return converted matrix, or NULL in case of allocation failure
 */
jmtxz_matrix_brm *jmtxz_matrix_brm_from_double_imag_inplace(JMTX_NAME_TYPED(matrix_brm) *in)
{
    const JMTX_FAST_INT_T count = brm_row_offsetd(in, in->base.rows);
    _Complex JMTX_SCALAR_T *const values =
        in->base.allocator_callbacks.realloc(in->base.allocator_callbacks.state, in->values, sizeof(*values) * count);
    if (!values)
    {
        return NULL;
    }
    in->values = (JMTX_SCALAR_T *)values;

    for (JMTX_FAST_INT_T i = 0; i < count; ++i)
    {
        values[count - 1 - i] = _Complex_I * in->values[count - 1 - i];
    }

    in->base.type = JMTXZ_TYPE_BRM;
    return (jmtxz_matrix_brm *)in;
}

/**
 * Creates a new complex BRM matrix with JMTX_SCALAR_T precision from a BRM matrix with JMTX_SCALAR_T precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in_real matrix to use as the real component of the new matrix (if NULL real part is zero)
 * @param in_imag matrix to use as the real component of the new matrix (if NULL imaginary part is zero)
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxzs_matrix_brm_from_double(jmtxz_matrix_brm **p_mtx, const JMTX_NAME_TYPED(matrix_brm) *in_real,
                                          const JMTX_NAME_TYPED(matrix_brm) *in_imag,
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
    if ((in_real && in_real->base.type != JMTXD_TYPE_BRM) || (in_imag && in_imag->base.type != JMTXD_TYPE_BRM))
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if ((in_real && in_imag) &&
        (in_real->base.rows != in_imag->base.rows || in_real->base.cols != in_imag->base.cols ||
         in_real->lower_bandwidth != in_imag->lower_bandwidth || in_real->upper_bandwidth != in_imag->upper_bandwidth))
    {
        return JMTX_RESULT_BAD_MATRIX;
    }

    return jmtxz_matrix_brm_from_double(p_mtx, in_real, in_imag, allocator_callbacks);
}

/***********************************************************************************************************************
 *                                                                                                                     *
 *                                 COMPLEX FLOAT <->  COMPLEX DOUBLE                                                   *
 *                                                                                                                     *
 **********************************************************************************************************************/
/**
 * Creates a new complex BRM matrix with single precision from a complex BRM matrix with JMTX_SCALAR_T precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxc_matrix_brm_from_cdouble(jmtxc_matrix_brm **p_mtx, const jmtxz_matrix_brm *in,
                                          const jmtx_allocator_callbacks *allocator_callbacks)
{
    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }

    jmtxc_matrix_brm *mtx = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*mtx));
    if (!mtx)
    {
        return JMTX_RESULT_BAD_ALLOC;
    }

    const JMTX_FAST_INT_T count = brm_row_offsetz(in, in->base.rows);
    _Complex float *values = allocator_callbacks->alloc(allocator_callbacks->state, count * sizeof(*values));
    if (!values)
    {
        allocator_callbacks->free(allocator_callbacks->state, mtx);
        return JMTX_RESULT_BAD_ALLOC;
    }

    for (JMTX_FAST_INT_T i = 0; i < count; ++i)
    {
        values[i] = (_Complex float)in->values[i];
    }

    mtx->base.cols = in->base.cols;
    mtx->base.type = JMTXC_TYPE_BRM;
    mtx->base.rows = in->base.rows;
    mtx->base.allocator_callbacks = *allocator_callbacks;
    mtx->values = values;
    mtx->lower_bandwidth = in->lower_bandwidth;
    mtx->upper_bandwidth = in->upper_bandwidth;

    *p_mtx = mtx;

    return JMTX_RESULT_SUCCESS;
}
/**
 * Creates a new complex BRM matrix with single precision from a complex BRM matrix with JMTX_SCALAR_T precision. Requires no
 * memory allocation by reusing the memory of the initial matrix. Can not fail if the input matrix is valid.
 * @param in matrix which to convert (will be invalid if function succeeds)
 * @return converted matrix
 */
jmtxc_matrix_brm *jmtxc_matrix_brm_from_cdouble_inplace(jmtxz_matrix_brm *in)
{
    const JMTX_FAST_INT_T count = brm_row_offsetz(in, in->base.rows);
    _Complex float *const values = (_Complex float *)in->values;

    for (JMTX_FAST_INT_T i = 0; i < count; ++i)
    {
        values[i] = (float)in->values[i];
    }

    _Complex float *const new_ptr =
        in->base.allocator_callbacks.realloc(in->base.allocator_callbacks.state, values, sizeof(*new_ptr) * count);
    if (new_ptr)
    {
        in->values = (_Complex JMTX_SCALAR_T *)new_ptr;
    }

    in->base.type = JMTXC_TYPE_BRM;
    return (jmtxc_matrix_brm *)in;
}

/**
 * Creates a new complex BRM matrix with single precision from a complex BRM matrix with JMTX_SCALAR_T precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxcs_matrix_brm_from_cdouble(jmtxc_matrix_brm **p_mtx, const jmtxz_matrix_brm *in,
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
    if (in->base.type != JMTXZ_TYPE_BRM)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }

    return jmtxc_matrix_brm_from_cdouble(p_mtx, in, allocator_callbacks);
}

/**
 * Creates a new complex BRM matrix with JMTX_SCALAR_T precision from a complex BRM matrix with single precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxz_matrix_brm_from_cfloat(jmtxz_matrix_brm **p_mtx, const jmtxc_matrix_brm *in,
                                         const jmtx_allocator_callbacks *allocator_callbacks)
{
    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }

    jmtxz_matrix_brm *mtx = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*mtx));
    if (!mtx)
    {
        return JMTX_RESULT_BAD_ALLOC;
    }

    const JMTX_FAST_INT_T count = brm_row_offsetc(in, in->base.rows);
    _Complex JMTX_SCALAR_T *values = allocator_callbacks->alloc(allocator_callbacks->state, count * sizeof(*values));
    if (!values)
    {
        allocator_callbacks->free(allocator_callbacks->state, mtx);
        return JMTX_RESULT_BAD_ALLOC;
    }

    for (JMTX_FAST_INT_T i = 0; i < count; ++i)
    {
        values[i] = (_Complex float)in->values[i];
    }

    mtx->base.cols = in->base.cols;
    mtx->base.type = JMTXZ_TYPE_BRM;
    mtx->base.rows = in->base.rows;
    mtx->base.allocator_callbacks = *allocator_callbacks;
    mtx->values = values;
    mtx->lower_bandwidth = in->lower_bandwidth;
    mtx->upper_bandwidth = in->upper_bandwidth;

    *p_mtx = mtx;

    return JMTX_RESULT_SUCCESS;
}

/**
 * Creates a new complex BRM matrix with JMTX_SCALAR_T precision from a complex BRM matrix with single precision. Only one
 * memory reallocation may be needed. Can not fail if the input matrix is valid.
 * @param in matrix which to convert
 * @return converted matrix, or NULL in case of allocation failure
 */
jmtxz_matrix_brm *jmtxz_matrix_brm_from_cfloat_inplace(jmtxc_matrix_brm *in)
{
    const JMTX_FAST_INT_T count = brm_row_offsetc(in, in->base.rows);
    _Complex JMTX_SCALAR_T *const values =
        in->base.allocator_callbacks.realloc(in->base.allocator_callbacks.state, in->values, sizeof(*values) * count);
    if (!values)
    {
        return NULL;
    }
    in->values = (_Complex float *)values;

    for (JMTX_FAST_INT_T i = 0; i < count; ++i)
    {
        values[count - 1 - i] = (_Complex double)in->values[count - 1 - i];
    }

    in->base.type = JMTXZ_TYPE_BRM;
    return (jmtxz_matrix_brm *)in;
}

/**
 * Creates a new complex BRM matrix with JMTX_SCALAR_T precision from a complex BRM matrix with single precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxzs_matrix_brm_from_cfloat(jmtxz_matrix_brm **p_mtx, const jmtxc_matrix_brm *in,
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
    if (in->base.type != JMTXC_TYPE_BRM)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }

    return jmtxz_matrix_brm_from_cfloat(p_mtx, in, allocator_callbacks);
}
#endif //!_MSC_BUILD
