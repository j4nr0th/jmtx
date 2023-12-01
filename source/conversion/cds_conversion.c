//
// Created by jan on 1.12.2023.
//

#include "../../include/jmtx/conversion/cds_conversion.h"
#include "../matrix_base_internal.h"
#include "../float/matrices/sparse_diagonal_compressed_internal.h"
#include "../double/matrices/sparse_diagonal_compressed_internal.h"

/**
 * Creates a new CDS matrix with single precision from a CDS matrix with double precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtx_matrix_cds_from_double(jmtx_matrix_cds** p_mtx, const jmtxd_matrix_cds* in,
                                        const jmtx_allocator_callbacks* allocator_callbacks)
{

    jmtx_matrix_cds* mtx;
    jmtx_result res = jmtx_matrix_cds_new(&mtx, in->base.cols, in->base.rows, 0, (int32_t[]) {0}, allocator_callbacks);
    if (res != JMTX_RESULT_SUCCESS)
    {
        return res;
    }
    
    for (uint_fast32_t i = 0; i < in->sub_diagonals.count; ++i)
    {
        const int32_t dia_idx = -(int32_t)in->sub_diagonals.indices[i];
        const double* const dia_ptr = in->sub_diagonals.diagonals[i];
        uint32_t len;
        float* const ptr = jmtx_matrix_cds_allocate_diagonal(mtx, dia_idx, &len);
        for (uint_fast32_t j = 0; j < len; ++j)
        {
            ptr[j] = (float)dia_ptr[j];
        }
    }
    if (in->main_diagonal)
    {
        uint32_t len;
        float* const ptr = jmtx_matrix_cds_allocate_diagonal(mtx, 0, &len);
        for (uint_fast32_t j = 0; j < len; ++j)
        {
            ptr[j] = (float)in->main_diagonal[j];
        }
    }
    for (uint_fast32_t i = 0; i < in->super_diagonals.count; ++i)
    {
        const int32_t dia_idx = (int32_t)in->super_diagonals.indices[i];
        const double* const dia_ptr = in->super_diagonals.diagonals[i];
        uint32_t len;
        float* const ptr = jmtx_matrix_cds_allocate_diagonal(mtx, dia_idx, &len);
        for (uint_fast32_t j = 0; j < len; ++j)
        {
            ptr[j] = (float)dia_ptr[j];
        }
    }

 

    mtx->base.cols = in->base.cols;
    mtx->base.type = JMTX_TYPE_CDS;
    mtx->base.rows = in->base.cols;
    mtx->base.allocator_callbacks = *allocator_callbacks;

    *p_mtx = mtx;

    return JMTX_RESULT_SUCCESS;
}

/**
 * Creates a new CDS matrix with single precision from a CDS matrix with double precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxs_matrix_cds_from_double(jmtx_matrix_cds** p_mtx, const jmtxd_matrix_cds* in,
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
    if (in->base.type != JMTXD_TYPE_CDS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }

    return jmtx_matrix_cds_from_double(p_mtx, in, allocator_callbacks);
}

/**
 * Creates a new CDS matrix with double precision from a CDS matrix with single precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxd_matrix_cds_from_float(jmtxd_matrix_cds** p_mtx, const jmtx_matrix_cds* in,
                                        const jmtx_allocator_callbacks* allocator_callbacks)
{

    jmtxd_matrix_cds* mtx;
    jmtx_result res = jmtxd_matrix_cds_new(&mtx, in->base.cols, in->base.rows, 0, (int32_t[]) {0}, allocator_callbacks);
    if (res != JMTX_RESULT_SUCCESS)
    {
        return res;
    }

    for (uint_fast32_t i = 0; i < in->sub_diagonals.count; ++i)
    {
        const int32_t dia_idx = -(int32_t)in->sub_diagonals.indices[i];
        const float* const dia_ptr = in->sub_diagonals.diagonals[i];
        uint32_t len;
        double * const ptr = jmtxd_matrix_cds_allocate_diagonal(mtx, dia_idx, &len);
        for (uint_fast32_t j = 0; j < len; ++j)
        {
            ptr[j] = (double)dia_ptr[j];
        }
    }
    if (in->main_diagonal)
    {
        uint32_t len;
        double* const ptr = jmtxd_matrix_cds_allocate_diagonal(mtx, 0, &len);
        for (uint_fast32_t j = 0; j < len; ++j)
        {
            ptr[j] = (double)in->main_diagonal[j];
        }
    }
    for (uint_fast32_t i = 0; i < in->super_diagonals.count; ++i)
    {
        const int32_t dia_idx = (int32_t)in->super_diagonals.indices[i];
        const float * const dia_ptr = in->super_diagonals.diagonals[i];
        uint32_t len;
        double* const ptr = jmtxd_matrix_cds_allocate_diagonal(mtx, dia_idx, &len);
        for (uint_fast32_t j = 0; j < len; ++j)
        {
            ptr[j] = (double)dia_ptr[j];
        }
    }



    mtx->base.cols = in->base.cols;
    mtx->base.type = JMTXD_TYPE_CDS;
    mtx->base.rows = in->base.cols;
    mtx->base.allocator_callbacks = *allocator_callbacks;

    *p_mtx = mtx;

    return JMTX_RESULT_SUCCESS;
}

/**
 * Creates a new CDS matrix with double precision from a CDS matrix with single precision.
 * @param p_mtx Pointer which receives the pointer to the new matrix
 * @param in matrix which to convert
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxds_matrix_cds_from_float(jmtxd_matrix_cds** p_mtx, const jmtx_matrix_cds* in,
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
    if (in->base.type != JMTX_TYPE_CDS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }

    return jmtxd_matrix_cds_from_float(p_mtx, in, allocator_callbacks);
}
