// Automatically generated from source/float/matrices/sparse_conversion.c on Thu Dec 14 17:58:20 2023
//
// Created by jan on 2.11.2023.
//

#include <assert.h>
#include "sparse_row_compressed_internal.h"
#include "sparse_column_compressed_internal.h"
#include "sparse_diagonal_compressed_internal.h"
#include "band_row_major_internal.h"
#include "../../../include/jmtx/double/matrices/sparse_row_compressed_safe.h"
#include "../../../include/jmtx/double/matrices/sparse_column_compressed_safe.h"
#include "../../../include/jmtx/double/matrices/sparse_conversion.h"
#include "../../../include/jmtx/double/matrices/sparse_conversion_safe.h"
#include "../../../include/jmtx/double/matrices/sparse_diagonal_compressed_safe.h"
#include "../../../include/jmtx/double/matrices/band_row_major_safe.h"

jmtxd_matrix_ccs* jmtxd_convert_crs_to_ccs_inplace_transpose(jmtxd_matrix_crs* in)
{
    //  Lol, lmao even
    in->base.type = JMTXD_TYPE_CCS;

    return (jmtxd_matrix_ccs*)in;
}

jmtxd_matrix_crs* jmtxd_convert_ccs_to_crs_inplace_transpose(jmtxd_matrix_ccs* in)
{
    //  Lol, lmao even
    in->base.type = JMTXD_TYPE_CRS;
    return (jmtxd_matrix_crs*)in;
}

/**
 * Converts a CRS matrix into the CCS format. Input matrix remains untouched.
 * @param in CRS matrix to convert
 * @param p_out pointer which receives the converted matrix
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on allocation failure
 */
jmtx_result jmtxd_convert_crs_to_ccs(
        const jmtxd_matrix_crs* in, jmtxd_matrix_ccs** p_out, const jmtx_allocator_callbacks* allocator_callbacks)
{
    jmtxd_matrix_crs* cpy;
    jmtx_result res = jmtxd_matrix_crs_transpose(in, &cpy, allocator_callbacks);
    if (res != JMTX_RESULT_SUCCESS)
    {
        return res;
    }
    *p_out = jmtxd_convert_crs_to_ccs_inplace_transpose(cpy);
    return JMTX_RESULT_SUCCESS;
}

/**
 * Converts a CRS matrix into the CCS format. Input matrix remains untouched.
 * @param in CRS matrix to convert
 * @param p_out pointer which receives the converted matrix
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on allocation failure
 */
jmtx_result jmtxds_convert_crs_to_ccs(const jmtxd_matrix_crs* in, jmtxd_matrix_ccs** p_out,
                                     const jmtx_allocator_callbacks* allocator_callbacks)
{
    if (!in)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!p_out)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (in->base.type != JMTXD_TYPE_CRS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }
    else if (allocator_callbacks->alloc == NULL || allocator_callbacks->free == NULL)
    {
        return JMTX_RESULT_BAD_PARAM;
    }
    return jmtxd_convert_crs_to_ccs(in, p_out, allocator_callbacks);
}
/**
 * Converts a CCS matrix into the CRS format. Input matrix remains untouched.
 * @param in CCS matrix to convert
 * @param p_out pointer which receives the converted matrix
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on allocation failure
 */
jmtx_result jmtxd_convert_ccs_to_crs(
        const jmtxd_matrix_ccs* in, jmtxd_matrix_crs** p_out, const jmtx_allocator_callbacks* allocator_callbacks)
{
    jmtxd_matrix_ccs* cpy;
    jmtx_result res = jmtxds_matrix_ccs_transpose(in, &cpy, allocator_callbacks);
    if (res != JMTX_RESULT_SUCCESS)
    {
        return res;
    }
    *p_out = jmtxd_convert_ccs_to_crs_inplace_transpose(cpy);
    return JMTX_RESULT_SUCCESS;
}

/**
 * Converts a CCS matrix into the CRS format. Input matrix remains untouched.
 * @param in CCS matrix to convert
 * @param p_out pointer which receives the converted matrix
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on allocation failure
 */
jmtx_result jmtxds_convert_ccs_to_crs(const jmtxd_matrix_ccs* in, jmtxd_matrix_crs** p_out,
                                    const jmtx_allocator_callbacks* allocator_callbacks)
{
    if (!in)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!p_out)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (in->base.type != JMTXD_TYPE_CCS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }
    else if (allocator_callbacks->alloc == NULL || allocator_callbacks->free == NULL)
    {
        return JMTX_RESULT_BAD_PARAM;
    }
    return jmtxd_convert_ccs_to_crs(in, p_out, allocator_callbacks);
}
/**
 * Converts a CDS matrix into the CRS format. Input matrix remains untouched.
 * @param in CDS matrix to convert
 * @param p_out pointer which receives the converted matrix
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on allocation failure
 */
jmtx_result jmtxd_convert_cds_to_crs(
        const jmtxd_matrix_cds* in, jmtxd_matrix_crs** p_out, const jmtx_allocator_callbacks* allocator_callbacks)
{
    const uint_fast32_t max_dia_len = in->base.rows > in->base.cols ? in->base.cols : in->base.rows;
    jmtxd_matrix_crs* mtx;
    const uint_fast32_t max_entries = max_dia_len * jmtxd_matrix_cds_diagonal_count(in);
    jmtx_result res = jmtxd_matrix_crs_new(
        &mtx, in->base.rows, in->base.cols, max_entries, allocator_callbacks);
    if (res != JMTX_RESULT_SUCCESS)
    {
        return res;
    }
    for (uint_fast32_t i = 0, p = 0; i < in->base.rows; ++i)
    {
        assert(p <= max_entries);
        p += jmtxd_matrix_cds_get_row(in, i, max_entries - p, mtx->values + p, mtx->indices + p);
        assert(p <= max_entries);
        mtx->end_of_row_offsets[i] = p;
    }
    mtx->n_entries = mtx->end_of_row_offsets[mtx->base.rows - 1];
    (void)jmtxd_matrix_crs_shrink(mtx);
    *p_out = mtx;
    return JMTX_RESULT_SUCCESS;
}

/**
 * Converts a CDS matrix into the CRS format. Input matrix remains untouched.
 * @param in CDS matrix to convert
 * @param p_out pointer which receives the converted matrix
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result jmtxds_convert_cds_to_crs(const jmtxd_matrix_cds* in, jmtxd_matrix_crs** p_out,
                                     const jmtx_allocator_callbacks* allocator_callbacks)
{
    if (!in)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!p_out)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (in->base.type != JMTXD_TYPE_CDS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }
    else if (allocator_callbacks->alloc == NULL || allocator_callbacks->free == NULL)
    {
        return JMTX_RESULT_BAD_PARAM;
    }
    return jmtxd_convert_cds_to_crs(in, p_out, allocator_callbacks);
}

/**
 * Converts a CDS matrix into the CCS format. Input matrix remains untouched.
 * @param in CDS matrix to convert
 * @param p_out pointer which receives the converted matrix
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on allocation failure
 */
jmtx_result jmtxd_convert_cds_to_ccs(const jmtxd_matrix_cds* in, jmtxd_matrix_ccs** p_out,
                                    const jmtx_allocator_callbacks* allocator_callbacks)
{
    const uint_fast32_t max_dia_len = in->base.rows > in->base.cols ? in->base.cols : in->base.rows;
    jmtxd_matrix_ccs* mtx;
    const uint_fast32_t max_entries = max_dia_len * jmtxd_matrix_cds_diagonal_count(in);
    jmtx_result res = jmtxd_matrix_ccs_new(
        &mtx, in->base.rows, in->base.cols, max_entries, allocator_callbacks);
    if (res != JMTX_RESULT_SUCCESS)
    {
        return res;
    }
    for (uint_fast32_t i = 0, p = 0; i < in->base.cols; ++i)
    {
        assert(p <= max_entries);
        p += jmtxd_matrix_cds_get_col(in, i, max_entries - p, mtx->values + p, mtx->indices + p);
        assert(p <= max_entries);
        mtx->end_of_column_offsets[i] = p;
    }
    mtx->n_entries = mtx->end_of_column_offsets[mtx->base.cols - 1];
    (void)jmtxd_matrix_ccs_shrink(mtx);
    *p_out = mtx;
    return JMTX_RESULT_SUCCESS;
}
/**
 * Converts a CDS matrix into the CCS format. Input matrix remains untouched.
 * @param in CDS matrix to convert
 * @param p_out pointer which receives the converted matrix
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result jmtxds_convert_cds_to_ccs(const jmtxd_matrix_cds* in, jmtxd_matrix_ccs** p_out,
                                     const jmtx_allocator_callbacks* allocator_callbacks)
{
    if (!in)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!p_out)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (in->base.type != JMTXD_TYPE_CDS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }
    else if (allocator_callbacks->alloc == NULL || allocator_callbacks->free == NULL)
    {
        return JMTX_RESULT_BAD_PARAM;
    }
    return jmtxd_convert_cds_to_ccs(in, p_out, allocator_callbacks);
}

/**
 * Converts a CRS matrix into the CRS format. Input matrix remains untouched.
 * @param in CRS matrix to convert
 * @param p_out pointer which receives the converted matrix
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on allocation failure
 */
jmtx_result jmtxd_convert_crs_to_cds(const jmtxd_matrix_crs* in, jmtxd_matrix_cds** p_out,
                                    const jmtx_allocator_callbacks* allocator_callbacks)
{
    jmtxd_matrix_cds* mtx;
    jmtx_result res = jmtxd_matrix_cds_new(
        &mtx, in->base.rows, in->base.cols, 0, (const int32_t[]) {0}, allocator_callbacks);
    if (res != JMTX_RESULT_SUCCESS)
    {
        return res;
    }
    for (uint_fast32_t i = 0; i < in->base.rows; ++i)
    {
        uint32_t* indices;
        double* values;
        const uint32_t count = jmtxd_matrix_crs_get_row(in, i, &indices, &values);
        res = jmtxd_matrix_cds_set_row(mtx, i, count, values, indices);
        if (res != JMTX_RESULT_SUCCESS)
        {
            jmtxd_matrix_cds_destroy(mtx);
            return res;
        }
    }
    *p_out = mtx;
    
    return JMTX_RESULT_SUCCESS;
}
/**
 * Converts a CRS matrix into the CDS format. Input matrix remains untouched.
 * @param in CRS matrix to convert
 * @param p_out pointer which receives the converted matrix
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result jmtxds_convert_crs_to_cds(const jmtxd_matrix_crs* in, jmtxd_matrix_cds** p_out,
                                     const jmtx_allocator_callbacks* allocator_callbacks)
{
    if (!in)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!p_out)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (in->base.type != JMTXD_TYPE_CRS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }
    else if (allocator_callbacks->alloc == NULL || allocator_callbacks->free == NULL)
    {
        return JMTX_RESULT_BAD_PARAM;
    }
    return jmtxd_convert_crs_to_cds(in, p_out, allocator_callbacks);
}

/**
 * Converts a CCS matrix into the CDS format. Input matrix remains untouched.
 * @param in CCS matrix to convert
 * @param p_out pointer which receives the converted matrix
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on allocation failure
 */
jmtx_result jmtxd_convert_ccs_to_cds(const jmtxd_matrix_ccs* in, jmtxd_matrix_cds** p_out,
                                    const jmtx_allocator_callbacks* allocator_callbacks)
{
    jmtxd_matrix_cds* mtx;
    jmtx_result res = jmtxd_matrix_cds_new(
        &mtx, in->base.rows, in->base.cols, 0, (const int32_t[]) {0}, allocator_callbacks);
    if (res != JMTX_RESULT_SUCCESS)
    {
        return res;
    }
    for (uint_fast32_t i = 0; i < in->base.cols; ++i)
    {
        uint32_t* indices;
        double* values;
        const uint32_t count = jmtxd_matrix_ccs_get_col(in, i, &indices, &values);
        res = jmtxd_matrix_cds_set_col(mtx, i, count, values, indices);
        if (res != JMTX_RESULT_SUCCESS)
        {
            jmtxd_matrix_cds_destroy(mtx);
            return res;
        }
    }
    *p_out = mtx;

    return JMTX_RESULT_SUCCESS;
}

/**
 * Converts a CCS matrix into the CDS format. Input matrix remains untouched.
 * @param in CDS matrix to convert
 * @param p_out pointer which receives the converted matrix
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result jmtxds_convert_ccs_to_cds(const jmtxd_matrix_ccs* in, jmtxd_matrix_cds** p_out,
                                     const jmtx_allocator_callbacks* allocator_callbacks)
{
    if (!in)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!p_out)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (in->base.type != JMTXD_TYPE_CCS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }
    else if (allocator_callbacks->alloc == NULL || allocator_callbacks->free == NULL)
    {
        return JMTX_RESULT_BAD_PARAM;
    }
    return jmtxd_convert_ccs_to_cds(in, p_out, allocator_callbacks);
}

/**
 * Converts a BRM matrix into the CRS format. Input matrix remains untouched.
 * @param in BRM matrix to convert
 * @param p_out pointer which receives the converted matrix
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on allocation failure
 */
jmtx_result jmtxd_convert_brm_to_crs(const jmtxd_matrix_brm* in, jmtxd_matrix_crs** p_out,
                                    const jmtx_allocator_callbacks* allocator_callbacks)
{
    const uint_fast32_t max_entries = (in->lower_bandwidth + in->upper_bandwidth + 1) * in->base.rows;
    jmtxd_matrix_crs* mtx;
    jmtx_result res = jmtxd_matrix_crs_new(&mtx, in->base.rows, in->base.cols, max_entries, allocator_callbacks);
    if (res != JMTX_RESULT_SUCCESS)
    {
        return res;
    }
    for (uint_fast32_t i = 0, p = 0; i < in->base.rows; ++i)
    {
        double* elements;
        const uint_fast32_t count = jmtxd_matrix_brm_get_row(in, i, &elements);
        memcpy(mtx->values + p, elements, sizeof(*elements) * count);
        const uint_fast32_t first = jmtxd_matrix_brm_first_pos_in_row(in, i);
        for (uint_fast32_t j = 0; j < count; ++j)
        {
            mtx->indices[i + j] = first + j;
        }
        assert(p <= max_entries);
        p += count;
        assert(p <= max_entries);
        mtx->end_of_row_offsets[i] = p;
    }
    mtx->n_entries = mtx->end_of_row_offsets[mtx->base.rows - 1];
    jmtxd_matrix_crs_shrink(mtx);
    *p_out = mtx;
    return JMTX_RESULT_SUCCESS;
}
/**
 * Converts a BRM matrix into the CRS format. Input matrix remains untouched.
 * @param in BRM matrix to convert
 * @param p_out pointer which receives the converted matrix
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result jmtxds_convert_brm_to_crs(const jmtxd_matrix_brm* in, jmtxd_matrix_crs** p_out,
                                     const jmtx_allocator_callbacks* allocator_callbacks)
{
    if (!in)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!p_out)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (in->base.type != JMTXD_TYPE_BRM)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }
    else if (allocator_callbacks->alloc == NULL || allocator_callbacks->free == NULL)
    {
        return JMTX_RESULT_BAD_PARAM;
    }
    return jmtxd_convert_brm_to_crs(in, p_out, allocator_callbacks);
}

/**
 * Converts a BRM matrix into the CCS format. Input matrix remains untouched.
 * @param in BRM matrix to convert
 * @param p_out pointer which receives the converted matrix
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on allocation failure
 */
jmtx_result jmtxd_convert_brm_to_ccs(const jmtxd_matrix_brm* in, jmtxd_matrix_ccs** p_out,
                                    const jmtx_allocator_callbacks* allocator_callbacks)
{
    const uint_fast32_t max_entries = (in->lower_bandwidth + in->upper_bandwidth + 1) * in->base.rows;
    jmtxd_matrix_ccs* mtx;
    jmtx_result res = jmtxd_matrix_ccs_new(&mtx, in->base.rows, in->base.cols, max_entries, allocator_callbacks);
    if (res != JMTX_RESULT_SUCCESS)
    {
        return res;
    }
    for (uint_fast32_t i = 0, p = 0; i < in->base.rows; ++i)
    {
        const uint_fast32_t count = jmtxd_matrix_brm_get_col(in, i, mtx->values + p);
        const uint_fast32_t first = jmtxd_matrix_brm_first_pos_in_col(in, i);
        for (uint_fast32_t j = 0; j < count; ++j)
        {
            mtx->indices[i + j] = first + j;
        }
        assert(p <= max_entries);
        p += count;
        assert(p <= max_entries);
        mtx->end_of_column_offsets[i] = p;
    }
    mtx->n_entries = mtx->end_of_column_offsets[mtx->base.cols - 1];
    jmtxd_matrix_ccs_shrink(mtx);
    *p_out = mtx;
    return JMTX_RESULT_SUCCESS;
}
/**
 * Converts a BRM matrix into the CCS format. Input matrix remains untouched.
 * @param in BRM matrix to convert
 * @param p_out pointer which receives the converted matrix
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result jmtxds_convert_brm_to_ccs(const jmtxd_matrix_brm* in, jmtxd_matrix_ccs** p_out,
                                     const jmtx_allocator_callbacks* allocator_callbacks)
{
    if (!in)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!p_out)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (in->base.type != JMTXD_TYPE_BRM)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }
    else if (allocator_callbacks->alloc == NULL || allocator_callbacks->free == NULL)
    {
        return JMTX_RESULT_BAD_PARAM;
    }
    return jmtxd_convert_brm_to_ccs(in, p_out, allocator_callbacks);
}
/**
 * Converts a BRM matrix into the CDS format. Input matrix remains untouched.
 * @param in BRM matrix to convert
 * @param p_out pointer which receives the converted matrix
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on allocation failure
 */
jmtx_result jmtxd_convert_brm_to_cds(const jmtxd_matrix_brm* in, jmtxd_matrix_cds** p_out,
                                    const jmtx_allocator_callbacks* allocator_callbacks)
{
    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }
    const uint_fast32_t max_entries = (in->lower_bandwidth + in->upper_bandwidth + 1) * in->base.rows;
    jmtxd_matrix_cds* mtx;
    uint32_t* const index_array = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*index_array) * max_entries);
    if (!index_array)
    {
        return JMTX_RESULT_BAD_ALLOC;
    }
    double* const value_array = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*value_array) * max_entries);
    if (!value_array)
    {
        allocator_callbacks->free(allocator_callbacks->state, index_array);
        return JMTX_RESULT_BAD_ALLOC;
    }
    jmtx_result res = jmtxd_matrix_cds_new(&mtx, in->base.rows, in->base.cols, 0, (const int32_t[]){0}, allocator_callbacks);
    if (res != JMTX_RESULT_SUCCESS)
    {
        allocator_callbacks->free(allocator_callbacks->state, value_array);
        allocator_callbacks->free(allocator_callbacks->state, index_array);
        return res;
    }
    for (uint_fast32_t i = 0, p = 0; i < in->base.rows; ++i)
    {
        double* elements;
        const uint_fast32_t count = jmtxd_matrix_brm_get_row(in, i, &elements);
        const uint_fast32_t first = jmtxd_matrix_brm_first_pos_in_row(in, i);
        uint_fast32_t k = 0;
        for (uint_fast32_t j = 0; j < count; ++j)
        {
            if (elements[j] != 0)
            {
                index_array[k] = j + first;
                value_array[k] = elements[j];
                k += 1;
            }
        }
        jmtxd_matrix_cds_set_row(mtx, i, count, value_array, index_array);
        assert(p <= max_entries);
        p += count;
        assert(p <= max_entries);

    }
    allocator_callbacks->free(allocator_callbacks->state, value_array);
    allocator_callbacks->free(allocator_callbacks->state, index_array);

    *p_out = mtx;
    return JMTX_RESULT_SUCCESS;
}
/**
 * Converts a BRM matrix into the CDS format. Input matrix remains untouched.
 * @param in BRM matrix to convert
 * @param p_out pointer which receives the converted matrix
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result jmtxds_convert_brm_to_cds(const jmtxd_matrix_brm* in, jmtxd_matrix_cds** p_out,
                                     const jmtx_allocator_callbacks* allocator_callbacks)
{
    if (!in)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!p_out)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (in->base.type != JMTXD_TYPE_BRM)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }
    else if (allocator_callbacks->alloc == NULL || allocator_callbacks->free == NULL)
    {
        return JMTX_RESULT_BAD_PARAM;
    }
    return jmtxd_convert_brm_to_cds(in, p_out, allocator_callbacks);
}
/**
 * Converts a CRS matrix into the BRM format. Input matrix remains untouched.
 * @param in CRS matrix to convert
 * @param p_out pointer which receives the converted matrix
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on allocation failure
 */
jmtx_result jmtxd_convert_crs_to_brm(const jmtxd_matrix_crs* in, jmtxd_matrix_brm** p_out,
                                    const jmtx_allocator_callbacks* allocator_callbacks)
{
    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }
    jmtxd_matrix_brm* mtx;
    const uint32_t lbw = jmtxd_matrix_crs_find_lower_bandwidth(in);
    const uint32_t ubw = jmtxd_matrix_crs_find_upper_bandwidth(in);
    jmtx_result res = jmtxd_matrix_brm_new(&mtx, in->base.rows, in->base.cols, ubw, lbw, NULL, allocator_callbacks);
    if (res != JMTX_RESULT_SUCCESS)
    {
        return res;
    }
    for (uint_fast32_t i = 0, p = 0; i < in->base.rows; ++i)
    {
        uint32_t* indices;
        double* elements;
        const uint_fast32_t count = jmtxd_matrix_crs_get_row(in, i, &indices, &elements);
        const uint_fast32_t row_len = jmtxd_matrix_brm_length_of_row(mtx, i);
        const uint_fast32_t first = jmtxd_matrix_brm_first_pos_in_row(mtx, i);
        for (uint32_t k = 0, l = 0; k < row_len; ++k)
        {
            if (l < count && indices[l] == k + first)
            {
                mtx->values[p + k] = elements[l];
                l += 1;
            }
            else
            {
                mtx->values[p + k] = 0;
            }
        }
        p += row_len;
    }

    *p_out = mtx;
    return JMTX_RESULT_SUCCESS;
}

/**
 * Converts a CRS matrix into the BRM format. Input matrix remains untouched.
 * @param in CRS matrix to convert
 * @param p_out pointer which receives the converted matrix
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result jmtxds_convert_crs_to_brm(const jmtxd_matrix_crs* in, jmtxd_matrix_brm** p_out,
                                     const jmtx_allocator_callbacks* allocator_callbacks)
{
    if (!in)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!p_out)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (in->base.type != JMTXD_TYPE_CRS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }
    else if (allocator_callbacks->alloc == NULL || allocator_callbacks->free == NULL)
    {
        return JMTX_RESULT_BAD_PARAM;
    }

    return jmtxd_convert_crs_to_brm(in, p_out, allocator_callbacks);
}

/**
 * Converts a CCS matrix into the BRM format. Input matrix remains untouched.
 * @param in CCS matrix to convert
 * @param p_out pointer which receives the converted matrix
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on allocation failure
 */
jmtx_result jmtxd_convert_ccs_to_brm(const jmtxd_matrix_ccs* in, jmtxd_matrix_brm** p_out,
                                    const jmtx_allocator_callbacks* allocator_callbacks)
{
    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }
    jmtxd_matrix_brm* mtx;
    const uint32_t lbw = jmtxd_matrix_ccs_find_lower_bandwidth(in);
    const uint32_t ubw = jmtxd_matrix_ccs_find_upper_bandwidth(in);
    double* const values = allocator_callbacks->alloc(allocator_callbacks->state, (lbw + ubw + 1) * sizeof(*values));
    if (!values)
    {
        return JMTX_RESULT_BAD_ALLOC;
    }
    jmtx_result res = jmtxd_matrix_brm_new(&mtx, in->base.rows, in->base.cols, ubw, lbw, NULL, allocator_callbacks);
    if (res != JMTX_RESULT_SUCCESS)
    {
        allocator_callbacks->free(allocator_callbacks->state, values);
        return res;
    }
    for (uint_fast32_t i = 0; i < in->base.cols; ++i)
    {
        uint32_t* indices;
        double* elements;
        const uint_fast32_t count = jmtxd_matrix_ccs_get_col(in, i, &indices, &elements);
        const uint_fast32_t first = jmtxd_matrix_brm_first_pos_in_col(mtx, i);
        const uint_fast32_t last = jmtxd_matrix_brm_last_pos_in_col(mtx, i);
        for (uint32_t k = first, l = 0; k <= last; ++k)
        {
            if (l < count && indices[l] == k)
            {
                values[k - first] = elements[l];
                l += 1;
            }
            else
            {
                values[k - first] = 0;
            }
        }
        jmtxd_matrix_brm_set_col(mtx, i, values);
    }

    allocator_callbacks->free(allocator_callbacks->state, values);
    *p_out = mtx;
    return JMTX_RESULT_SUCCESS;
}

/**
 * Converts a CCS matrix into the BRM format. Input matrix remains untouched.
 * @param in CCS matrix to convert
 * @param p_out pointer which receives the converted matrix
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result jmtxds_convert_ccs_to_brm(const jmtxd_matrix_ccs* in, jmtxd_matrix_brm** p_out,
                                     const jmtx_allocator_callbacks* allocator_callbacks)
{
    if (!in)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!p_out)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (in->base.type != JMTXD_TYPE_CCS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }
    else if (allocator_callbacks->alloc == NULL || allocator_callbacks->free == NULL)
    {
        return JMTX_RESULT_BAD_PARAM;
    }
    return jmtxd_convert_ccs_to_brm(in, p_out, allocator_callbacks);
}

/**
 * Converts a CDS matrix into the BRM format. Input matrix remains untouched.
 * @param in CDS matrix to convert
 * @param p_out pointer which receives the converted matrix
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on allocation failure
 */
jmtx_result jmtxd_convert_cds_to_brm(const jmtxd_matrix_cds* in, jmtxd_matrix_brm** p_out,
                                    const jmtx_allocator_callbacks* allocator_callbacks)
{
    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }
    jmtxd_matrix_brm* mtx;
    const uint32_t lbw = in->sub_diagonals.count == 0 ? 0 : in->sub_diagonals.indices[in->sub_diagonals.count - 1];
    const uint32_t ubw = in->super_diagonals.count == 0 ? 0 : in->super_diagonals.indices[in->super_diagonals.count - 1];
    double* const values = allocator_callbacks->alloc(allocator_callbacks->state, (lbw + ubw + 1) * sizeof(*values));
    if (!values)
    {
        return JMTX_RESULT_BAD_ALLOC;
    }
    double* const values_in = allocator_callbacks->alloc(allocator_callbacks->state, (lbw + ubw + 1) * sizeof(*values_in));
    if (!values_in)
    {
        allocator_callbacks->free(allocator_callbacks->state, values);
        return JMTX_RESULT_BAD_ALLOC;
    }
    uint32_t* const indices = allocator_callbacks->alloc(allocator_callbacks->state, (lbw + ubw + 1) * sizeof(*indices));
    if (!indices)
    {
        allocator_callbacks->free(allocator_callbacks->state, values_in);
        allocator_callbacks->free(allocator_callbacks->state, values);
        return JMTX_RESULT_BAD_ALLOC;
    }
    jmtx_result res = jmtxd_matrix_brm_new(&mtx, in->base.rows, in->base.cols, ubw, lbw, NULL, allocator_callbacks);
    if (res != JMTX_RESULT_SUCCESS)
    {
        allocator_callbacks->free(allocator_callbacks->state, indices);
        allocator_callbacks->free(allocator_callbacks->state, values_in);
        allocator_callbacks->free(allocator_callbacks->state, values);
        return res;
    }
    for (uint_fast32_t i = 0; i < in->base.cols; ++i)
    {
        const uint_fast32_t count = jmtxd_matrix_cds_get_row(in, i, (1 + lbw + ubw), values_in, indices);
        const uint_fast32_t row_len = jmtxd_matrix_brm_length_of_row(mtx, i);
        const uint_fast32_t first = jmtxd_matrix_brm_first_pos_in_row(mtx, i);
        for (uint32_t k = 0, l = 0; k < row_len; ++k)
        {
            if (l < count && indices[l] == k + first)
            {
                values[k] = values_in[l];
                l += 1;
            }
            else
            {
                values[k] = 0;
            }
        }
        jmtxd_matrix_brm_set_row(mtx, i, values);
    }

    allocator_callbacks->free(allocator_callbacks->state, indices);
    allocator_callbacks->free(allocator_callbacks->state, values_in);
    allocator_callbacks->free(allocator_callbacks->state, values);
    *p_out = mtx;
    return JMTX_RESULT_SUCCESS;
}

/**
 * Converts a CDS matrix into the BRM format. Input matrix remains untouched.
 * @param in CDS matrix to convert
 * @param p_out pointer which receives the converted matrix
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result jmtxds_convert_cds_to_brm(const jmtxd_matrix_cds* in, jmtxd_matrix_brm** p_out,
                                     const jmtx_allocator_callbacks* allocator_callbacks)
{
    if (!in)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!p_out)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (in->base.type != JMTXD_TYPE_CDS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }
    else if (allocator_callbacks->alloc == NULL || allocator_callbacks->free == NULL)
    {
        return JMTX_RESULT_BAD_PARAM;
    }
    return jmtxd_convert_cds_to_brm(in, p_out, allocator_callbacks);
}
