//
// Created by jan on 1.12.2023.
//

#include "../../include/jmtx/conversion/brm_conversion.h"
#include "../matrix_base_internal.h"
#include "../float/matrices/band_row_major_internal.h"
#include "../double/matrices/band_row_major_internal.h"

static inline uint_fast32_t brm_row_offset(const jmtx_matrix_brm mtx[const static 1], const uint_fast32_t row)
{
    uint_fast32_t offset = (mtx->lower_bandwidth + 1 + mtx->upper_bandwidth) * row;
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
        const uint_fast32_t offset_row = row - (mtx->base.rows - 1 - mtx->upper_bandwidth);   //  rows offset from
        offset -= offset_row * (offset_row - 1) / 2;
    }
    return offset;
}

static inline uint_fast32_t brm_row_offsetd(const jmtxd_matrix_brm mtx[const static 1], const uint_fast32_t row)
{
    uint_fast32_t offset = (mtx->lower_bandwidth + 1 + mtx->upper_bandwidth) * row;
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
        const uint_fast32_t offset_row = row - (mtx->base.rows - 1 - mtx->upper_bandwidth);   //  rows offset from
        offset -= offset_row * (offset_row - 1) / 2;
    }
    return offset;
}
/**
 * Creates a new BRM matrix with single precision from a BRM matrix with double precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtx_matrix_brm_from_double(jmtx_matrix_brm** p_mtx, const jmtxd_matrix_brm* in,
                                        const jmtx_allocator_callbacks* allocator_callbacks)
{
    if (!allocator_callbacks)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }


    jmtx_matrix_brm* mtx = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*mtx));
    if (!mtx)
    {
        return JMTX_RESULT_BAD_ALLOC;
    }

    const uint_fast32_t count = brm_row_offsetd(in, in->base.rows);
    float* values = allocator_callbacks->alloc(allocator_callbacks->state, count * sizeof(*values));
    if (!values)
    {
        allocator_callbacks->free(allocator_callbacks->state, mtx);
        return JMTX_RESULT_BAD_ALLOC;
    }
    
    for (uint_fast32_t i = 0; i < count; ++i)
    {
        values[i] = (float)in->values[i];
    }

    mtx->base.cols = in->base.cols;
    mtx->base.type = JMTX_TYPE_BRM;
    mtx->base.rows = in->base.cols;
    mtx->base.allocator_callbacks = *allocator_callbacks;
    mtx->values = values;
    mtx->lower_bandwidth = in->lower_bandwidth;
    mtx->upper_bandwidth = in->upper_bandwidth;

    *p_mtx = mtx;


    return JMTX_RESULT_SUCCESS;
}

/**
 * Creates a new BRM matrix with single precision from a BRM matrix with double precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxs_matrix_brm_from_double(jmtx_matrix_brm** p_mtx, const jmtxd_matrix_brm* in,
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
    if (in->base.type != JMTXD_TYPE_BRM)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }

    return jmtx_matrix_brm_from_double(p_mtx, in, allocator_callbacks);
}

/**
 * Creates a new BRM matrix with double precision from a BRM matrix with single precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxd_matrix_brm_from_float(jmtxd_matrix_brm** p_mtx, const jmtx_matrix_brm* in,
                                        const jmtx_allocator_callbacks* allocator_callbacks)
{
    if (!allocator_callbacks)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }


    jmtxd_matrix_brm* mtx = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*mtx));
    if (!mtx)
    {
        return JMTX_RESULT_BAD_ALLOC;
    }

    const uint_fast32_t count = brm_row_offset(in, in->base.rows);
    double* values = allocator_callbacks->alloc(allocator_callbacks->state, count * sizeof(*values));
    if (!values)
    {
        allocator_callbacks->free(allocator_callbacks->state, mtx);
        return JMTX_RESULT_BAD_ALLOC;
    }

    for (uint_fast32_t i = 0; i < count; ++i)
    {
        values[i] = (float)in->values[i];
    }

    mtx->base.cols = in->base.cols;
    mtx->base.type = JMTXD_TYPE_BRM;
    mtx->base.rows = in->base.cols;
    mtx->base.allocator_callbacks = *allocator_callbacks;
    mtx->values = values;
    mtx->lower_bandwidth = in->lower_bandwidth;
    mtx->upper_bandwidth = in->upper_bandwidth;

    *p_mtx = mtx;


    return JMTX_RESULT_SUCCESS;
}

/**
 * Creates a new BRM matrix with double precision from a BRM matrix with single precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxds_matrix_brm_from_float(jmtxd_matrix_brm** p_mtx, const jmtx_matrix_brm* in,
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
    if (in->base.type != JMTX_TYPE_BRM)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }

    return jmtxd_matrix_brm_from_float(p_mtx, in, allocator_callbacks);
}
