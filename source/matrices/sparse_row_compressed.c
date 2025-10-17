// Automatically generated from source/float/matrices/sparse_row_compressed.c on Fri Dec  1 06:43:01 2023
//
// Created by jan on 13.6.2022.
//

#include "sparse_row_compressed.h"
#include <assert.h>
#include <math.h>

static JMTX_INDEX_T crs_get_row_entries(const JMTX_NAME_TYPED(matrix_crs) * mtx, JMTX_INDEX_T row,
                                        JMTX_INDEX_T *pp_indices[1], JMTX_SCALAR_T *pp_values[1])
{
    JMTX_INDEX_T offset, len;
    if (row == 0)
    {
        offset = 0;
        len = mtx->end_of_row_offsets[0];
    }
    else
    {
        offset = mtx->end_of_row_offsets[row - 1];
        len = mtx->end_of_row_offsets[row] - offset;
    }
    *pp_indices = mtx->indices + offset;
    *pp_values = mtx->values + offset;

    return len;
}

/**
 * Inserts an entry (value-index pair) into the matrix at a specified row and a global position.
 * @param mtx matrix to which to add the entry
 * @param row what row the element is going to be in
 * @param position what is the position in the row
 * @param value what is the value of the entry
 * @param index the index of the entry
 * @return JMTX_RESULT_BAD_ALLOC on allocation failure
 */
static jmtx_result crs_insert_entry_at(JMTX_NAME_TYPED(matrix_crs) * mtx, JMTX_INDEX_T row, JMTX_INDEX_T position,
                                       JMTX_SCALAR_T value, JMTX_INDEX_T index)
{
    const JMTX_INDEX_T global_position = position + (row ? mtx->end_of_row_offsets[row - 1] : 0);
    assert(position <= mtx->n_entries);
    assert(!row || (mtx->end_of_row_offsets[row - 1] == global_position || mtx->indices[global_position - 1] < index));

    if (mtx->capacity == mtx->n_entries)
    {
        //  Reallocate arrays
        const JMTX_SIZE_T new_capacity = mtx->capacity + DEFAULT_RESERVED_ELEMENTS;
        JMTX_SCALAR_T *const new_values = mtx->base.allocator_callbacks.realloc(
            mtx->base.allocator_callbacks.state, mtx->values, sizeof(*mtx->values) * new_capacity);
        if (!new_values)
        {
            return JMTX_RESULT_BAD_ALLOC;
        }
        mtx->values = new_values;
        JMTX_INDEX_T *const new_indices = mtx->base.allocator_callbacks.realloc(
            mtx->base.allocator_callbacks.state, mtx->indices, sizeof(*mtx->indices) * new_capacity);
        if (!new_indices)
        {
            return JMTX_RESULT_BAD_ALLOC;
        }
        mtx->indices = new_indices;
        mtx->capacity = (JMTX_INDEX_T)new_capacity;
    }

    //  Check for number of values after the position
    const JMTX_INDEX_T elements_after = mtx->n_entries - global_position;
    if (elements_after)
    {
        //  Move other values out of the way
        memmove(mtx->values + global_position + 1, mtx->values + global_position,
                sizeof(*mtx->values) * (mtx->n_entries - global_position));
        memmove(mtx->indices + global_position + 1, mtx->indices + global_position,
                sizeof(*mtx->indices) * (mtx->n_entries - global_position));
    }
    //  Insert the values
    mtx->values[global_position] = value;
    mtx->indices[global_position] = index;

    //  Update offsets
    for (JMTX_INDEX_T i = row; i < mtx->base.rows; ++i)
    {
        mtx->end_of_row_offsets[i] += 1;
    }

    mtx->n_entries += 1;

    return JMTX_RESULT_SUCCESS;
}

jmtx_result JMTX_NAME_TYPED(matrix_crs_new)(JMTX_NAME_TYPED(matrix_crs) * *p_mtx, JMTX_INDEX_T rows, JMTX_INDEX_T cols,
                                            JMTX_INDEX_T reserved_entries,
                                            const jmtx_allocator_callbacks *allocator_callbacks)
{
    if (reserved_entries == 0)
    {
        reserved_entries = DEFAULT_RESERVED_ELEMENTS;
        reserved_entries =
            reserved_entries < ((uint64_t)cols * (uint64_t)rows) ? reserved_entries : ((uint64_t)cols * (uint64_t)rows);
    }
    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }

    ;
    JMTX_INDEX_T *offsets = NULL;
    JMTX_INDEX_T *indices = NULL;

    JMTX_NAME_TYPED(matrix_crs) *mtx = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*mtx));
    if (!mtx)
    {
        return JMTX_RESULT_BAD_ALLOC;
    }

    JMTX_SCALAR_T *values =
        allocator_callbacks->alloc(allocator_callbacks->state, (reserved_entries) * sizeof(*values));
    if (!values)
    {
        allocator_callbacks->free(allocator_callbacks->state, mtx);
        return JMTX_RESULT_BAD_ALLOC;
    }
    memset(values, 0, (reserved_entries) * sizeof(*values));

    indices = allocator_callbacks->alloc(allocator_callbacks->state, (reserved_entries) * sizeof(*indices));
    if (!indices)
    {
        allocator_callbacks->free(allocator_callbacks->state, indices);
        allocator_callbacks->free(allocator_callbacks->state, mtx);
        return JMTX_RESULT_BAD_ALLOC;
    }
    memset(indices, 0, (reserved_entries) * sizeof(*indices));

    offsets = allocator_callbacks->alloc(allocator_callbacks->state, (rows) * sizeof(*offsets));
    if (!offsets)
    {
        allocator_callbacks->free(allocator_callbacks->state, offsets);
        allocator_callbacks->free(allocator_callbacks->state, indices);
        allocator_callbacks->free(allocator_callbacks->state, mtx);
        return JMTX_RESULT_BAD_ALLOC;
    }
    memset(offsets, 0, (rows) * sizeof(*offsets));

    memset(mtx, 0, sizeof *mtx);
    mtx->base.cols = cols;
    mtx->base.type = JMTXD_TYPE_CRS;
    mtx->base.rows = rows;
    mtx->base.allocator_callbacks = *allocator_callbacks;
    mtx->indices = indices;
    mtx->values = values;
    mtx->capacity = reserved_entries;
    mtx->n_entries = 0;
    mtx->end_of_row_offsets = offsets;
    *p_mtx = mtx;

    return JMTX_RESULT_SUCCESS;
}

void JMTX_NAME_TYPED(matrix_crs_destroy)(JMTX_NAME_TYPED(matrix_crs) * mtx)
{
    jmtx_allocator_callbacks allocator = mtx->base.allocator_callbacks;
    allocator.free(allocator.state, mtx->indices);
    allocator.free(allocator.state, mtx->end_of_row_offsets);
    allocator.free(allocator.state, mtx->values);
    allocator.free(allocator.state, mtx);
}

jmtx_result JMTX_NAME_TYPED(matrix_crs_shrink)(JMTX_NAME_TYPED(matrix_crs) * mtx)
{
    if (mtx->n_entries == mtx->capacity)
    {
        return JMTX_RESULT_SUCCESS;
    }

    JMTX_SCALAR_T *element_new_ptr = mtx->base.allocator_callbacks.realloc(
        mtx->base.allocator_callbacks.state, mtx->values, sizeof *mtx->values * (mtx->n_entries));
    if (!element_new_ptr)
    {
        return JMTX_RESULT_BAD_ALLOC;
    }
    mtx->values = element_new_ptr;
    JMTX_INDEX_T *new_indices_ptr = mtx->base.allocator_callbacks.realloc(
        mtx->base.allocator_callbacks.state, mtx->indices, sizeof *mtx->indices * (mtx->n_entries));
    if (!new_indices_ptr)
    {
        return JMTX_RESULT_BAD_ALLOC;
    }
    mtx->indices = new_indices_ptr;
    mtx->capacity = mtx->n_entries;

    return JMTX_RESULT_SUCCESS;
}

jmtx_result JMTX_NAME_TYPED(matrix_crs_set_row)(JMTX_NAME_TYPED(matrix_crs) * mtx, JMTX_INDEX_T row, JMTX_INDEX_T n,
                                                const JMTX_INDEX_T indices[JMTX_ARRAY_ATTRIB(static n)],
                                                const JMTX_SCALAR_T values[JMTX_ARRAY_ATTRIB(static n)])
{
    jmtx_result res = JMTX_RESULT_SUCCESS;
    const JMTX_INDEX_T beginning_offset = row ? mtx->end_of_row_offsets[row - 1] : 0;
    const int32_t new_elements = (int32_t)n - (int32_t)(mtx->end_of_row_offsets[row] - beginning_offset);
    const JMTX_INDEX_T required_capacity = (JMTX_INDEX_T)((int32_t)mtx->n_entries + new_elements);
    if (mtx->capacity < required_capacity)
    {
        JMTX_SCALAR_T *new_element_ptr = mtx->base.allocator_callbacks.realloc(
            mtx->base.allocator_callbacks.state, mtx->values, sizeof *(mtx->values) * (required_capacity + 1));
        if (!new_element_ptr)
        {
            return JMTX_RESULT_BAD_ALLOC;
        }
        mtx->values = new_element_ptr;
        JMTX_INDEX_T *new_indices_ptr = mtx->base.allocator_callbacks.realloc(
            mtx->base.allocator_callbacks.state, mtx->indices, sizeof *(mtx->indices) * (required_capacity + 1));
        if (!new_indices_ptr)
        {
            return JMTX_RESULT_BAD_ALLOC;
        }
        mtx->indices = new_indices_ptr;
        mtx->capacity = required_capacity;
    }

    if (new_elements != 0)
    {
        const JMTX_INDEX_T elements_after = mtx->n_entries - mtx->end_of_row_offsets[row];
        if (elements_after)
        {
            memmove(mtx->values + mtx->end_of_row_offsets[row] + new_elements,
                    mtx->values + mtx->end_of_row_offsets[row], sizeof *mtx->values * (elements_after));
            memmove(mtx->indices + mtx->end_of_row_offsets[row] + new_elements,
                    mtx->indices + mtx->end_of_row_offsets[row], sizeof *mtx->indices * (elements_after));
        }
        memcpy(mtx->values + beginning_offset, values, sizeof *values * n);
        memcpy(mtx->indices + beginning_offset, indices, sizeof *indices * n);

        for (JMTX_INDEX_T i = row; i < mtx->base.rows; ++i)
        {
            mtx->end_of_row_offsets[i] += new_elements;
        }
        mtx->n_entries += new_elements;
    }
    else
    {
        memcpy(mtx->values + beginning_offset, values, sizeof *values * n);
        memcpy(mtx->indices + beginning_offset, indices, sizeof *indices * n);
    }

    return res;
}
void JMTX_NAME_TYPED(matrix_crs_vector_multiply)(const JMTX_NAME_TYPED(matrix_crs) * mtx,
                                                 const JMTX_SCALAR_T *restrict x, JMTX_SCALAR_T *restrict y)
{
    for (JMTX_INDEX_T i = 0; i < mtx->base.rows; ++i)
    {
        JMTX_INDEX_T *indices;
        JMTX_SCALAR_T *values;
        const JMTX_INDEX_T row_entries = crs_get_row_entries(mtx, i, &indices, &values);
        JMTX_SCALAR_T v = 0;
        for (JMTX_INDEX_T j = 0; j < row_entries; ++j)
        {
            const JMTX_INDEX_T k = indices[j];
            v += values[j] * x[k];
        }
        y[i] = v;
    }
}

jmtx_result JMTX_NAME_TYPED(matrix_crs_set_entry)(JMTX_NAME_TYPED(matrix_crs) * mtx, JMTX_INDEX_T i, JMTX_INDEX_T j,
                                                  JMTX_SCALAR_T value)
{
    jmtx_result res;
    JMTX_INDEX_T *row_indices;
    JMTX_SCALAR_T *row_values;
    const JMTX_INDEX_T n_row_elements = crs_get_row_entries(mtx, i, &row_indices, &row_values);
    //  Check if row has any values
    if (n_row_elements != 0)
    {
        //  Find first column entry less or equal to it
        const JMTX_INDEX_T possible = jmtx_internal_find_last_leq_value(n_row_elements, row_indices, j);
        if (row_indices[possible] == j)
        {
            row_values[possible] = value;
            return JMTX_RESULT_SUCCESS;
        }
        else
        {
            //  Figure out where to insert
            if (row_indices[possible] <= j)
            {
                //  Insert right after the *possible*
                res = crs_insert_entry_at(mtx, i, possible + 1, value, j);
            }
            else
            {
                //  Insert at the position of *possible*
                res = crs_insert_entry_at(mtx, i, possible, value, j);
            }
        }
    }
    else
    {
        //  Insert the new element in an empty row
        res = crs_insert_entry_at(mtx, i, 0, value, j);
    }

    return res;
}

JMTX_SCALAR_T JMTX_NAME_TYPED(matrix_crs_get_entry)(const JMTX_NAME_TYPED(matrix_crs) * mtx, JMTX_INDEX_T i,
                                                    JMTX_INDEX_T j)
{
    JMTX_INDEX_T *row_indices;
    JMTX_SCALAR_T *row_values;
    const JMTX_INDEX_T n_row_elements = crs_get_row_entries(mtx, i, &row_indices, &row_values);
    //  Check if row has any values
    if (n_row_elements != 0)
    {
        //  Find first column entry less or equal to it
        const JMTX_INDEX_T possible = jmtx_internal_find_last_leq_value(n_row_elements, row_indices, j);
        if (row_indices[possible] == j)
        {
            return row_values[possible];
        }
    }

    return 0;
}

JMTX_INDEX_T JMTX_NAME_TYPED(matrix_crs_get_row)(const JMTX_NAME_TYPED(matrix_crs) * mtx, JMTX_INDEX_T row,
                                                 JMTX_INDEX_T *p_indices[1], JMTX_SCALAR_T *p_elements[1])
{
    return crs_get_row_entries(mtx, row, p_indices, p_elements);
}

// Maybe remove this one. It lacks and all purpose
JMTX_INDEX_T JMTX_NAME_TYPED(matrix_crs_count_values)(const JMTX_NAME_TYPED(matrix_crs) * mtx, JMTX_SCALAR_T v)
{
    JMTX_INDEX_T r = 0;
    for (JMTX_INDEX_T i = 0; i < mtx->n_entries; ++i)
    {
        if (v == mtx->values[i])
        {
            r += 1;
        }
    }
    return r;
}

JMTX_INDEX_T JMTX_NAME_TYPED(matrix_crs_count_indices)(const JMTX_NAME_TYPED(matrix_crs) * mtx, JMTX_INDEX_T v)
{
    JMTX_INDEX_T r = 0;
    for (JMTX_INDEX_T i = 0; i < mtx->n_entries; ++i)
    {
        if (v == mtx->indices[i])
        {
            r += 1;
        }
    }
    return r;
}

jmtx_result JMTX_NAME_TYPED(matrix_crs_apply_unary_fn)(const JMTX_NAME_TYPED(matrix_crs) * mtx,
                                                       int (*unary_fn)(JMTX_INDEX_T i, JMTX_INDEX_T j,
                                                                       JMTX_SCALAR_T *p_value, void *param),
                                                       void *param)
{
    for (JMTX_INDEX_T i = 0; i < mtx->base.rows; ++i)
    {
        JMTX_SCALAR_T *p_elements;
        JMTX_INDEX_T *p_indices;
        const JMTX_INDEX_T n_in_row = crs_get_row_entries(mtx, i, &p_indices, &p_elements);
        for (JMTX_INDEX_T j = 0; j < n_in_row; ++j)
        {
            if ((unary_fn(i, p_indices[j], p_elements + j, param)))
            {

                return JMTX_RESULT_UNARY_RETURN;
            }
        }
    }

    return JMTX_RESULT_SUCCESS;
}

void JMTX_NAME_TYPED(matrix_crs_remove_zeros)(JMTX_NAME_TYPED(matrix_crs) * mtx)
{
    //  Update offsets
    JMTX_INDEX_T p, c = 0, r = 0;
    while (mtx->end_of_row_offsets[r] == 0)
    {
        r += 1;
    }
    for (p = 0; p < mtx->n_entries; ++p)
    {
        assert(r < mtx->base.rows);
        if (mtx->values[p] == 0)
        {
            c += 1;
        }
        while (r < mtx->base.cols && p + 1 == mtx->end_of_row_offsets[r])
        {
            assert(mtx->end_of_row_offsets[r] >= c);
            mtx->end_of_row_offsets[r] -= c;
            r += 1;
        }
    }
    assert(r == mtx->base.rows);
    assert(c <= mtx->n_entries);
#ifndef NDEBUG
    for (r = 1; r < mtx->base.rows; ++r)
    {
        assert(mtx->end_of_row_offsets[r - 1] <= mtx->end_of_row_offsets[r]);
    }
#endif

    //  Remove the actual entries
    JMTX_INDEX_T p0 = mtx->n_entries;
    while (p0 != 0)
    {
        //  Check if entry must be removed
        if (mtx->values[p0 - 1] == 0)
        {
            JMTX_INDEX_T p1 = p0;
            while (p0 != 0 && mtx->values[p0 - 1] == 0)
            {
                p0 -= 1;
            }
            memmove(mtx->values + p0, mtx->values + p1, sizeof(*mtx->values) * (mtx->n_entries - p1));
            memmove(mtx->indices + p0, mtx->indices + p1, sizeof(*mtx->indices) * (mtx->n_entries - p1));
        }
        else
        {
            p0 -= 1;
        }
    }

    mtx->n_entries -= c;
}

void JMTX_NAME_TYPED(matrix_crs_remove_bellow_magnitude)(JMTX_NAME_TYPED(matrix_crs) * mtx, JMTX_REAL_T v)
{
    //  Update offsets
    JMTX_INDEX_T p, c = 0, r = 0;
    while (mtx->end_of_row_offsets[r] == 0)
    {
        r += 1;
    }
    for (p = 0; p < mtx->n_entries; ++p)
    {
        assert(r < mtx->base.rows);
        if (JMTX_ABS(mtx->values[p]) < v)
        {
            c += 1;
        }
        while (r < mtx->base.cols && p + 1 == mtx->end_of_row_offsets[r])
        {
            assert(mtx->end_of_row_offsets[r] >= c);
            mtx->end_of_row_offsets[r] -= c;
            r += 1;
        }
    }
    assert(r == mtx->base.rows);
    assert(c <= mtx->n_entries);
#ifndef NDEBUG
    for (r = 1; r < mtx->base.rows; ++r)
    {
        assert(mtx->end_of_row_offsets[r - 1] <= mtx->end_of_row_offsets[r]);
    }
#endif

    //  Remove the actual entries
    JMTX_INDEX_T p0 = mtx->n_entries;
    while (p0 != 0)
    {
        //  Check if entry must be removed
        if (JMTX_ABS(mtx->values[p0 - 1]) < v)
        {
            JMTX_INDEX_T p1 = p0;
            while (p0 != 0 && JMTX_ABS(mtx->values[p0 - 1]) < v)
            {
                p0 -= 1;
            }
            memmove(mtx->values + p0, mtx->values + p1, sizeof(*mtx->values) * (mtx->n_entries - p1));
            memmove(mtx->indices + p0, mtx->indices + p1, sizeof(*mtx->indices) * (mtx->n_entries - p1));
        }

        p0 -= 1;
    }

    mtx->n_entries -= c;
}

JMTX_INDEX_T JMTX_NAME_TYPED(matrix_crs_entries_in_col)(const JMTX_NAME_TYPED(matrix_crs) * mtx, JMTX_INDEX_T col)
{
    JMTX_INDEX_T element_count = 0;
    for (JMTX_INDEX_T row = 0; row < mtx->base.rows && (!row || (mtx->end_of_row_offsets[row - 1] != mtx->n_entries));
         ++row)
    {
        JMTX_INDEX_T *row_indices;
        JMTX_SCALAR_T *unused_row_values;
        const JMTX_INDEX_T n_row_elements = crs_get_row_entries(mtx, row, &row_indices, &unused_row_values);
        if (n_row_elements && row_indices[0] <= col && row_indices[n_row_elements - 1] >= col)
        {
            JMTX_INDEX_T current = jmtx_internal_find_last_leq_value(n_row_elements, row_indices, col);
            if (row_indices[current] == col)
            {
                element_count += 1;
            }
        }
    }
    return element_count;
}

JMTX_INDEX_T JMTX_NAME_TYPED(matrix_crs_get_col)(const JMTX_NAME_TYPED(matrix_crs) * mtx, JMTX_INDEX_T col,
                                                 JMTX_INDEX_T n, JMTX_SCALAR_T p_values[JMTX_ARRAY_ATTRIB(n)],
                                                 JMTX_INDEX_T p_rows[JMTX_ARRAY_ATTRIB(n)])
{
    JMTX_INDEX_T k = 0;
    for (JMTX_INDEX_T row = 0;
         k < n && row < mtx->base.rows && (!row || (mtx->end_of_row_offsets[row - 1] != mtx->n_entries)); ++row)
    {
        JMTX_INDEX_T *row_indices;
        JMTX_SCALAR_T *unused_row_values;
        const JMTX_INDEX_T n_row_elements = crs_get_row_entries(mtx, row, &row_indices, &unused_row_values);
        if (n_row_elements && row_indices[0] <= col && row_indices[n_row_elements - 1] >= col)
        {
            JMTX_INDEX_T current = jmtx_internal_find_last_leq_value(n_row_elements, row_indices, col);
            if (row_indices[current] == col)
            {
                p_values[k] = mtx->values[(row ? mtx->end_of_row_offsets[row - 1] : 0) + current];
                p_rows[k] = row;
                k += 1;
            }
        }
    }
    return k;
}

jmtx_result JMTX_NAME_TYPED(matrix_crs_transpose)(const JMTX_NAME_TYPED(matrix_crs) * mtx,
                                                  JMTX_NAME_TYPED(matrix_crs) * *p_out,
                                                  const jmtx_allocator_callbacks *allocator_callbacks)
{
    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }

    const JMTX_INDEX_T cols = mtx->base.cols;
    JMTX_NAME_TYPED(matrix_crs) * out;
    jmtx_result res =
        JMTX_NAME_TYPED(matrix_crs_new)(&out, mtx->base.cols, mtx->base.rows, mtx->n_entries, allocator_callbacks);
    if (res != JMTX_RESULT_SUCCESS)
    {
        return res;
    }
    if (!allocator_callbacks)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }

    JMTX_INDEX_T *col_counts = allocator_callbacks->alloc(allocator_callbacks->state, sizeof *col_counts * cols);
    if (col_counts == NULL)
    {
        JMTX_NAME_TYPED(matrix_crs_destroy)(out);
        return JMTX_RESULT_BAD_ALLOC;
    }
    memset(col_counts, 0, sizeof *col_counts * cols);

    JMTX_INDEX_T *row_ends = out->end_of_row_offsets;
    for (JMTX_INDEX_T i = 0; i < mtx->n_entries; ++i)
    {
        col_counts[mtx->indices[i]] += 1;
    }
    row_ends[0] = col_counts[0];
    //  Compute cumsums for offsets
    for (JMTX_INDEX_T i = 1; i < cols; ++i)
    {
        row_ends[i] = col_counts[i] + row_ends[i - 1];
        col_counts[i] = 0; //   Zero the row counts so that they can be reused later for counting bucket sizes
    }
    col_counts[0] = 0;
    col_counts[cols - 1] = 0;

    for (JMTX_INDEX_T row = 0; row < mtx->base.rows; ++row)
    {
        JMTX_INDEX_T *in_cols;
        JMTX_SCALAR_T *in_vals;
        JMTX_INDEX_T n_col = crs_get_row_entries(mtx, row, &in_cols, &in_vals);

        for (JMTX_INDEX_T idx = 0; idx < n_col; ++idx)
        {
            const JMTX_INDEX_T col = in_cols[idx];
            const JMTX_INDEX_T ip = col > 0 ? row_ends[col - 1] : 0;
            const JMTX_INDEX_T n_cols = col_counts[col];

            out->values[ip + n_cols] = in_vals[idx];
            out->indices[ip + n_cols] = row;
            col_counts[col] += 1;
        }
    }
    out->n_entries = mtx->n_entries;

    allocator_callbacks->free(allocator_callbacks->state, col_counts);
    *p_out = out;

    return JMTX_RESULT_SUCCESS;
}

jmtx_result JMTX_NAME_TYPED(matrix_crs_copy)(const JMTX_NAME_TYPED(matrix_crs) * mtx,
                                             JMTX_NAME_TYPED(matrix_crs) * *p_out,
                                             const jmtx_allocator_callbacks *allocator_callbacks)
{
    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }
    JMTX_NAME_TYPED(matrix_crs) *const out = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*out));
    if (!out)
    {
        return JMTX_RESULT_BAD_ALLOC;
    }
    JMTX_SCALAR_T *const elements =
        allocator_callbacks->alloc(allocator_callbacks->state, (mtx->n_entries) * sizeof(*elements));
    if (!elements)
    {
        allocator_callbacks->free(allocator_callbacks->state, out);
        return JMTX_RESULT_BAD_ALLOC;
    }
    JMTX_INDEX_T *const indices =
        allocator_callbacks->alloc(allocator_callbacks->state, (mtx->n_entries) * sizeof *indices);
    if (!indices)
    {
        allocator_callbacks->free(allocator_callbacks->state, out);
        allocator_callbacks->free(allocator_callbacks->state, elements);
        return JMTX_RESULT_BAD_ALLOC;
    }
    JMTX_INDEX_T *const cum_sum =
        allocator_callbacks->alloc(allocator_callbacks->state, (mtx->base.rows) * sizeof *cum_sum);
    if (!cum_sum)
    {
        allocator_callbacks->free(allocator_callbacks->state, out);
        allocator_callbacks->free(allocator_callbacks->state, indices);
        allocator_callbacks->free(allocator_callbacks->state, elements);
        return JMTX_RESULT_BAD_ALLOC;
    }

    memcpy(elements, mtx->values, sizeof *elements * mtx->n_entries);
    memcpy(indices, mtx->indices, sizeof *indices * mtx->n_entries);
    memcpy(cum_sum, mtx->end_of_row_offsets, sizeof *cum_sum * (mtx->base.rows));
    memcpy(out, mtx, sizeof *out);
    out->values = elements;
    out->indices = indices;
    out->end_of_row_offsets = cum_sum;
    out->base = mtx->base;
    out->capacity = out->n_entries;
    out->base.allocator_callbacks = *allocator_callbacks;
    *p_out = out;
    return JMTX_RESULT_SUCCESS;
}

jmtx_result JMTX_NAME_TYPED(matrix_crs_build_row)(JMTX_NAME_TYPED(matrix_crs) * mtx, JMTX_INDEX_T row, JMTX_INDEX_T n,
                                                  const JMTX_INDEX_T indices[JMTX_ARRAY_ATTRIB(static n)],
                                                  const JMTX_SCALAR_T values[JMTX_ARRAY_ATTRIB(static n)])
{
    jmtx_result res = JMTX_RESULT_SUCCESS;
    const JMTX_INDEX_T required_capacity = (JMTX_INDEX_T)((int32_t)mtx->n_entries + (int32_t)n);
    if (mtx->capacity < required_capacity)
    {
        JMTX_SCALAR_T *new_element_ptr = mtx->base.allocator_callbacks.realloc(
            mtx->base.allocator_callbacks.state, mtx->values, sizeof *mtx->values * (required_capacity + 1));
        if (!new_element_ptr)
        {
            return JMTX_RESULT_BAD_ALLOC;
        }
        mtx->values = new_element_ptr;
        JMTX_INDEX_T *new_indices_ptr = mtx->base.allocator_callbacks.realloc(
            mtx->base.allocator_callbacks.state, mtx->indices, sizeof *mtx->indices * (required_capacity + 1));
        if (!new_indices_ptr)
        {
            return JMTX_RESULT_BAD_ALLOC;
        }
        mtx->indices = new_indices_ptr;
        mtx->capacity = required_capacity;
    }

    const JMTX_INDEX_T offset = row ? mtx->end_of_row_offsets[row - 1] : 0;
    memcpy(mtx->values + offset, values, sizeof *values * n);
    memcpy(mtx->indices + offset, indices, sizeof *indices * n);

    mtx->end_of_row_offsets[row] = n + offset;
    mtx->n_entries += n;

    return res;
}

JMTX_SCALAR_T JMTX_NAME_TYPED(matrix_crs_vector_multiply_row)(const JMTX_NAME_TYPED(matrix_crs) * mtx,
                                                              const JMTX_SCALAR_T *x, JMTX_INDEX_T i)
{
    JMTX_INDEX_T *indices;
    JMTX_SCALAR_T *values;
    const JMTX_INDEX_T n_row = crs_get_row_entries(mtx, i, &indices, &values);
    JMTX_SCALAR_T v = 0;
    for (JMTX_INDEX_T j = 0; j < n_row; ++j)
    {
        v += values[j] * x[indices[j]];
    }
    return v;
}

jmtx_result JMTX_NAME_TYPED(matrix_crs_add_to_entry)(JMTX_NAME_TYPED(matrix_crs) * mtx, JMTX_INDEX_T i, JMTX_INDEX_T j,
                                                     JMTX_SCALAR_T value)
{
    jmtx_result res;
    JMTX_INDEX_T *row_indices;
    JMTX_SCALAR_T *row_values;
    const JMTX_INDEX_T n_row_elements = crs_get_row_entries(mtx, i, &row_indices, &row_values);
    //  Check if row has any values
    if (n_row_elements != 0)
    {
        //  Find first column entry less or equal to it
        const JMTX_INDEX_T possible = jmtx_internal_find_last_leq_value(n_row_elements, row_indices, j);
        if (row_indices[possible] == j)
        {
            row_values[possible] += value;
            return JMTX_RESULT_SUCCESS;
        }
        else
        {
            //  Figure out where to insert
            if (row_indices[possible] <= j)
            {
                //  Insert right after the *possible*
                res = crs_insert_entry_at(mtx, i, possible + 1, value, j);
            }
            else
            {
                //  Insert at the position of *possible*
                res = crs_insert_entry_at(mtx, i, possible, value, j);
            }
        }
    }
    else
    {
        //  Insert the new element in an empty row
        res = crs_insert_entry_at(mtx, i, 0, value, j);
    }

    return res;
}

void JMTX_NAME_TYPED(matrix_crs_zero_all_entries)(const JMTX_NAME_TYPED(matrix_crs) * mtx)
{
    memset(mtx->values, 0, sizeof(*mtx->values) * mtx->n_entries);
}

void JMTX_NAME_TYPED(matrix_crs_set_all_entries)(const JMTX_NAME_TYPED(matrix_crs) * mtx, JMTX_SCALAR_T x)
{
    for (JMTX_SCALAR_T *ptr = mtx->values; ptr != mtx->values + mtx->n_entries; ++ptr)
    {
        *ptr = x;
    }
}

void JMTX_NAME_TYPED(matrix_crs_remove_row)(JMTX_NAME_TYPED(matrix_crs) * mtx, JMTX_INDEX_T row)
{
    JMTX_INDEX_T *row_indices;
    JMTX_SCALAR_T *row_values;
    const JMTX_INDEX_T removed_entry_count = crs_get_row_entries(mtx, row, &row_indices, &row_values);
    const JMTX_INDEX_T end_offset = mtx->end_of_row_offsets[row];
    //  Check if row is not empty
    if (removed_entry_count != 0)
    {
        //  Check if there are any entries following the ones being removed
        const JMTX_INDEX_T following_entries = mtx->n_entries - end_offset;
        if (following_entries != 0)
        {
            //  Move other values to their new position
            memmove(mtx->values + end_offset - removed_entry_count, mtx->values + end_offset,
                    sizeof(*mtx->values) * following_entries);
            memmove(mtx->indices + end_offset - removed_entry_count, mtx->indices + end_offset,
                    sizeof(*mtx->indices) * following_entries);
        }
        for (JMTX_INDEX_T i = row; i < mtx->base.rows - 1; ++i)
        {
            mtx->end_of_row_offsets[i] = mtx->end_of_row_offsets[i + 1] - removed_entry_count;
        }
        mtx->n_entries -= removed_entry_count;
    }
    mtx->base.rows -= 1;
}

void JMTX_NAME_TYPED(matrix_crs_remove_column)(JMTX_NAME_TYPED(matrix_crs) * mtx, JMTX_INDEX_T col)
{
    //  i: tracks the current element to check
    //  j: how many values to be removed were found
    //  r: what row we are currently in
    JMTX_INDEX_T i, j, r = 0;
    while (mtx->end_of_row_offsets[r] == 0)
    {
        r += 1;
    }
    for (i = 0, j = 0; i < mtx->n_entries; ++i)
    {
        if (mtx->indices[i] == col)
        {
            j += 1;
        }
        else
        {
            if (mtx->indices[i] > col)
            {
                mtx->indices[i] -= 1;
            }
            if (j != 0)
            {
                mtx->indices[i - j] = mtx->indices[i];
                mtx->values[i - j] = mtx->values[i];
            }
        }
        while (r < mtx->base.cols && i + 1 == mtx->end_of_row_offsets[r])
        {
            assert(mtx->end_of_row_offsets[r] >= j);
            mtx->end_of_row_offsets[r] -= j;
            r += 1;
        }
    }
    assert(mtx->n_entries >= j);
    assert(r == mtx->base.rows);
    mtx->n_entries -= j;
    mtx->base.cols -= 1;
}

void JMTX_NAME_TYPED(matrix_crs_clear)(JMTX_NAME_TYPED(matrix_crs) * mtx)
{
    mtx->n_entries = 0;
    memset(mtx->end_of_row_offsets, 0, sizeof(*mtx->end_of_row_offsets) * mtx->base.rows);
}

jmtx_result JMTX_NAME_TYPED(matrix_crs_join_vertically)(JMTX_NAME_TYPED(matrix_crs) * *output,
                                                        const jmtx_allocator_callbacks *allocator_callbacks, unsigned k,
                                                        const JMTX_NAME_TYPED(matrix_crs) *
                                                            matrix_list[JMTX_ARRAY_ATTRIB(static k)])
{
    const JMTX_INDEX_T col_count = matrix_list[0]->base.cols;
    JMTX_INDEX_T n_rows = matrix_list[0]->base.rows;

    JMTX_INDEX_T element_count = matrix_list[0]->n_entries;
    //  First one is already accounted for, so no need to start at 0
    for (unsigned i = 1; i < k; ++i)
    {
        const JMTX_NAME_TYPED(matrix_crs) *const e = matrix_list[i];
        if (e->base.cols != col_count)
        {
            return JMTX_RESULT_DIMS_MISMATCH;
        }
        n_rows += e->base.rows;
        element_count += e->n_entries;
    }

    JMTX_NAME_TYPED(matrix_crs) * out;
    jmtx_result res = JMTX_NAME_TYPED(matrix_crs_new)(&out, n_rows, col_count, element_count, allocator_callbacks);
    if (res != JMTX_RESULT_SUCCESS)
    {
        return res;
    }

    JMTX_INDEX_T pos = 0;
    for (unsigned i = 0; i < k; ++i)
    {
        const JMTX_NAME_TYPED(matrix_crs) *const mtx = matrix_list[i];
        unsigned j;
        for (j = 0; j < mtx->base.rows; ++j)
        {
            JMTX_INDEX_T *indices;
            JMTX_SCALAR_T *values;
            const JMTX_INDEX_T n = crs_get_row_entries(mtx, j, &indices, &values);
            //  All memory for this should be allocated in advance, so no need to check the return value
            const jmtx_result r1 = JMTX_NAME_TYPED(matrix_crs_build_row)(out, pos + j, n, indices, values);
            assert(r1 == JMTX_RESULT_SUCCESS);
            (void)r1;
        }
        pos += j;
    }
    assert(pos == n_rows);

    *output = out;
    return JMTX_RESULT_SUCCESS;
}

jmtx_result JMTX_NAME_TYPED(matrix_crs_new_like)(const JMTX_NAME_TYPED(matrix_crs) * mtx,
                                                 JMTX_NAME_TYPED(matrix_crs) * *p_out,
                                                 const jmtx_allocator_callbacks *allocator_callbacks,
                                                 const JMTX_SCALAR_T *p_val)
{
    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }
    JMTX_NAME_TYPED(matrix_crs) *const out = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*out));
    if (!out)
    {
        return JMTX_RESULT_BAD_ALLOC;
    }
    JMTX_SCALAR_T *const elements =
        allocator_callbacks->alloc(allocator_callbacks->state, (mtx->n_entries) * sizeof(*elements));
    if (!elements)
    {
        allocator_callbacks->free(allocator_callbacks->state, out);
        return JMTX_RESULT_BAD_ALLOC;
    }
    JMTX_INDEX_T *const indices =
        allocator_callbacks->alloc(allocator_callbacks->state, (mtx->n_entries) * sizeof *indices);
    if (!indices)
    {
        allocator_callbacks->free(allocator_callbacks->state, out);
        allocator_callbacks->free(allocator_callbacks->state, elements);
        return JMTX_RESULT_BAD_ALLOC;
    }
    JMTX_INDEX_T *const cum_sum =
        allocator_callbacks->alloc(allocator_callbacks->state, (mtx->base.rows) * sizeof *cum_sum);
    if (!cum_sum)
    {
        allocator_callbacks->free(allocator_callbacks->state, out);
        allocator_callbacks->free(allocator_callbacks->state, indices);
        allocator_callbacks->free(allocator_callbacks->state, elements);
        return JMTX_RESULT_BAD_ALLOC;
    }

    if (p_val)
    {
        const JMTX_SCALAR_T v = *p_val;
        if (v == 0)
        {
            memset(elements, 0, sizeof *elements * mtx->n_entries);
        }
        else
        {
            for (JMTX_INDEX_T i = 0; i < mtx->n_entries; ++i)
            {
                elements[i] = v;
            }
        }
    }

    memcpy(indices, mtx->indices, sizeof *indices * mtx->n_entries);
    memcpy(cum_sum, mtx->end_of_row_offsets, sizeof *cum_sum * (mtx->base.rows));
    memcpy(out, mtx, sizeof *out);
    out->values = elements;
    out->indices = indices;
    out->end_of_row_offsets = cum_sum;
    out->base = mtx->base;
    out->base.allocator_callbacks = *allocator_callbacks;
    *p_out = out;
    return JMTX_RESULT_SUCCESS;
}

JMTX_INDEX_T JMTX_NAME_TYPED(matrix_crs_find_upper_bandwidth)(const JMTX_NAME_TYPED(matrix_crs) * mtx)
{
    //  Find the greatest distance above the main diagonal
    JMTX_FAST_INT_T v_max = 0;
    for (JMTX_FAST_INT_T i = 0, p = 0; i < mtx->base.rows; ++i)
    {
        for (p = 0; p < mtx->end_of_row_offsets[i]; ++p)
        {
            if (mtx->indices[p] > i)
            {
                const JMTX_FAST_INT_T dif = mtx->indices[p] - i;
                if (dif > v_max)
                {
                    v_max = dif;
                }
            }
        }
    }
    return v_max;
}

/**
 * Finds the lower bandwidth of the matrix; what is the furthest distance of and entry below the main diagonal
 * @param mtx matrx to find the lower bandwidth of
 * @return lower bandwidth of the matrix
 */
JMTX_INDEX_T JMTX_NAME_TYPED(matrix_crs_find_lower_bandwidth)(const JMTX_NAME_TYPED(matrix_crs) * mtx)
{
    //  Find the greatest distance above the main diagonal
    JMTX_FAST_INT_T v_max = 0;
    for (JMTX_FAST_INT_T i = 0, p = 0; i < mtx->base.rows; ++i)
    {
        for (; p < mtx->end_of_row_offsets[i]; ++p)
        {
            if (mtx->indices[p] < i)
            {
                const JMTX_FAST_INT_T dif = i - mtx->indices[p];
                if (dif > v_max)
                {
                    v_max = dif;
                }
            }
        }
    }
    return v_max;
}
