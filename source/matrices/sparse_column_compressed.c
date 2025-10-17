
#include "sparse_column_compressed.h"
#include <assert.h>
#include <math.h>

static JMTX_INDEX_T ccs_get_column_entries(const JMTX_NAME_TYPED(matrix_ccs) * mtx, JMTX_INDEX_T col,
                                           JMTX_INDEX_T **pp_indices, JMTX_SCALAR_T **pp_values)
{
    JMTX_INDEX_T offset, len;
    if (col == 0)
    {
        offset = 0;
        len = mtx->end_of_column_offsets[0];
    }
    else
    {
        offset = mtx->end_of_column_offsets[col - 1];
        len = mtx->end_of_column_offsets[col] - offset;
    }
    *pp_indices = mtx->indices + offset;
    *pp_values = mtx->values + offset;

    return len;
}

/**
 * Inserts an entry (value-index pair) into the matrix at a specified column and a global position.
 * @param mtx matrix to which to add the entry
 * @param col what column the element is going to be in
 * @param position what is the position in the row
 * @param value what is the value of the entry
 * @param index the index of the entry
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC if memory allocation failed
 */
static jmtx_result ccs_insert_entry_at(JMTX_NAME_TYPED(matrix_ccs) * mtx, JMTX_INDEX_T col, JMTX_INDEX_T position,
                                       JMTX_SCALAR_T value, JMTX_INDEX_T index)
{
    const JMTX_INDEX_T global_position = position + (col ? mtx->end_of_column_offsets[col - 1] : 0);
    assert(position <= mtx->n_entries);
    assert(!col ||
           (mtx->end_of_column_offsets[col - 1] == global_position || mtx->indices[global_position - 1] < index));

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
    for (JMTX_INDEX_T i = col; i < mtx->base.cols; ++i)
    {
        mtx->end_of_column_offsets[i] += 1;
    }

    mtx->n_entries += 1;

    return JMTX_RESULT_SUCCESS;
}

jmtx_result JMTX_NAME_TYPED(matrix_ccs_new)(JMTX_NAME_TYPED(matrix_ccs) * *p_mtx, JMTX_INDEX_T rows, JMTX_INDEX_T cols,
                                            JMTX_INDEX_T reserved_entries,
                                            const jmtx_allocator_callbacks *allocator_callbacks)
{
    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }

    if (reserved_entries == 0)
    {
        reserved_entries = DEFAULT_RESERVED_ELEMENTS;
        reserved_entries =
            reserved_entries < ((uint64_t)cols * (uint64_t)rows) ? reserved_entries : ((uint64_t)cols * (uint64_t)rows);
    }

    JMTX_INDEX_T *offsets = NULL;
    JMTX_INDEX_T *indices = NULL;
    JMTX_SCALAR_T *values =
        allocator_callbacks->alloc(allocator_callbacks->state, (reserved_entries) * sizeof(*values));
    if (!values)
    {
        return JMTX_RESULT_BAD_ALLOC;
    }
    memset(values, 0, (reserved_entries) * sizeof(*values));

    indices = allocator_callbacks->alloc(allocator_callbacks->state, (reserved_entries) * sizeof *indices);
    if (!indices)
    {
        allocator_callbacks->free(allocator_callbacks->state, values);
        return JMTX_RESULT_BAD_ALLOC;
    }
    memset(indices, 0, (reserved_entries) * sizeof(*indices));

    offsets = allocator_callbacks->alloc(allocator_callbacks->state, (cols) * sizeof *offsets);
    if (!offsets)
    {
        allocator_callbacks->free(allocator_callbacks->state, indices);
        allocator_callbacks->free(allocator_callbacks->state, values);
        return JMTX_RESULT_BAD_ALLOC;
    }
    memset(offsets, 0, (cols) * sizeof(*offsets));

    JMTX_NAME_TYPED(matrix_ccs) *const this = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*this));
    if (!this)
    {
        allocator_callbacks->free(allocator_callbacks->state, offsets);
        allocator_callbacks->free(allocator_callbacks->state, indices);
        allocator_callbacks->free(allocator_callbacks->state, values);
        return JMTX_RESULT_BAD_ALLOC;
    }

    this->base.rows = rows;
    this->base.type = JMTXD_TYPE_CCS;
    this->base.cols = cols;
    this->base.allocator_callbacks = *allocator_callbacks;
    this->indices = indices;
    this->values = values;
    this->capacity = reserved_entries;
    this->n_entries = 0;
    this->end_of_column_offsets = offsets;

    *p_mtx = this;
    return JMTX_RESULT_SUCCESS;
}

void JMTX_NAME_TYPED(matrix_ccs_destroy)(JMTX_NAME_TYPED(matrix_ccs) * mtx)
{
    mtx->base.allocator_callbacks.free(mtx->base.allocator_callbacks.state, mtx->indices);
    mtx->base.allocator_callbacks.free(mtx->base.allocator_callbacks.state, mtx->end_of_column_offsets);
    mtx->base.allocator_callbacks.free(mtx->base.allocator_callbacks.state, mtx->values);
    mtx->base.allocator_callbacks.free(mtx->base.allocator_callbacks.state, mtx);
}

jmtx_result JMTX_NAME_TYPED(matrix_ccs_shrink)(JMTX_NAME_TYPED(matrix_ccs) * mtx)
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

jmtx_result JMTX_NAME_TYPED(matrix_ccs_set_col)(JMTX_NAME_TYPED(matrix_ccs) * mtx, JMTX_INDEX_T col, JMTX_INDEX_T n,
                                                const JMTX_INDEX_T *indices, const JMTX_SCALAR_T *values)
{
    const JMTX_INDEX_T beginning_offset = col ? mtx->end_of_column_offsets[col - 1] : 0;
    const int32_t new_elements = (int32_t)n - (int32_t)(mtx->end_of_column_offsets[col] - beginning_offset);
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
        const JMTX_INDEX_T elements_after = mtx->n_entries - mtx->end_of_column_offsets[col];
        if (elements_after)
        {
            memmove(mtx->values + mtx->end_of_column_offsets[col] + new_elements,
                    mtx->values + mtx->end_of_column_offsets[col], sizeof *mtx->values * (elements_after));
            memmove(mtx->indices + mtx->end_of_column_offsets[col] + new_elements,
                    mtx->indices + mtx->end_of_column_offsets[col], sizeof *mtx->indices * (elements_after));
        }
        memcpy(mtx->values + beginning_offset, values, sizeof *values * n);
        memcpy(mtx->indices + beginning_offset, indices, sizeof *indices * n);

        for (JMTX_INDEX_T i = col; i < mtx->base.cols; ++i)
        {
            mtx->end_of_column_offsets[i] += new_elements;
        }
        mtx->n_entries += new_elements;
    }
    else
    {
        memcpy(mtx->values + beginning_offset, values, sizeof *values * n);
        memcpy(mtx->indices + beginning_offset, indices, sizeof *indices * n);
    }
    return JMTX_RESULT_SUCCESS;
}

void JMTX_NAME_TYPED(matrix_ccs_vector_multiply)(const JMTX_NAME_TYPED(matrix_ccs) * mtx,
                                                 const JMTX_SCALAR_T *restrict x, JMTX_SCALAR_T *restrict y)
{
    for (JMTX_INDEX_T i = 0; i < mtx->base.cols; ++i)
    {
        JMTX_INDEX_T *indices;
        JMTX_SCALAR_T *values;
        const JMTX_INDEX_T n_elements = ccs_get_column_entries(mtx, i, &indices, &values);
        JMTX_SCALAR_T v = 0;
        for (JMTX_INDEX_T j = 0; j < n_elements; ++j)
        {
            const JMTX_INDEX_T k = indices[j];
            v += values[j] * x[k];
        }
        y[i] = v;
    }
}

jmtx_result JMTX_NAME_TYPED(matrix_ccs_set_entry)(JMTX_NAME_TYPED(matrix_ccs) * mtx, JMTX_INDEX_T i, JMTX_INDEX_T j,
                                                  JMTX_SCALAR_T value)
{
    jmtx_result res;
    JMTX_INDEX_T *col_indices;
    JMTX_SCALAR_T *col_values;
    const JMTX_INDEX_T n_col_elements = ccs_get_column_entries(mtx, j, &col_indices, &col_values);
    //  Check if row has any values
    if (n_col_elements != 0)
    {
        //  Find first column entry less or equal to it
        const JMTX_INDEX_T possible = jmtx_internal_find_last_leq_value(n_col_elements, col_indices, i);
        if (col_indices[possible] == i)
        {
            col_values[possible] = value;
            return JMTX_RESULT_SUCCESS;
        }
        else
        {
            //  Figure out where to insert
            if (col_indices[possible] <= i)
            {
                //  Insert right after the *possible*
                res = ccs_insert_entry_at(mtx, j, possible + 1, value, i);
            }
            else
            {
                //  Insert at the position of *possible*
                res = ccs_insert_entry_at(mtx, j, possible, value, i);
            }
        }
    }
    else
    {
        //  Insert the new element in an empty row
        res = ccs_insert_entry_at(mtx, j, 0, value, i);
    }

    return res;
}

JMTX_SCALAR_T JMTX_NAME_TYPED(matrix_ccs_get_entry)(const JMTX_NAME_TYPED(matrix_ccs) * mtx, JMTX_INDEX_T i,
                                                    JMTX_INDEX_T j)
{
    JMTX_INDEX_T *col_indices;
    JMTX_SCALAR_T *col_values;
    const JMTX_INDEX_T n_col_elements = ccs_get_column_entries(mtx, j, &col_indices, &col_values);
    //  Check if row has any values
    if (n_col_elements != 0)
    {
        //  Find first column entry less or equal to it
        const JMTX_INDEX_T possible = jmtx_internal_find_last_leq_value(n_col_elements, col_indices, i);
        if (col_indices[possible] == i)
        {
            return col_values[possible];
        }
    }
    return 0.0f;
}

JMTX_INDEX_T JMTX_NAME_TYPED(matrix_ccs_get_col)(const JMTX_NAME_TYPED(matrix_ccs) * mtx, JMTX_INDEX_T col,
                                                 JMTX_INDEX_T **p_indices, JMTX_SCALAR_T **p_elements)
{
    return ccs_get_column_entries(mtx, col, p_indices, p_elements);
}

JMTX_INDEX_T JMTX_NAME_TYPED(matrix_ccs_count_values)(const JMTX_NAME_TYPED(matrix_ccs) * mtx, JMTX_SCALAR_T v)
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

JMTX_INDEX_T JMTX_NAME_TYPED(matrix_ccs_count_indices)(const JMTX_NAME_TYPED(matrix_ccs) * mtx, JMTX_INDEX_T v)
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

jmtx_result JMTX_NAME_TYPED(matrix_ccs_apply_unary_fn)(const JMTX_NAME_TYPED(matrix_ccs) * mtx,
                                                       int (*unary_fn)(JMTX_INDEX_T i, JMTX_INDEX_T j,
                                                                       JMTX_SCALAR_T *p_element, void *param),
                                                       void *param)
{
    for (JMTX_INDEX_T j = 0; j < mtx->base.cols; ++j)
    {
        JMTX_SCALAR_T *p_elements;
        JMTX_INDEX_T *p_indices;
        const JMTX_INDEX_T n_in_row = ccs_get_column_entries(mtx, j, &p_indices, &p_elements);
        for (JMTX_INDEX_T i = 0; i < n_in_row; ++i)
        {
            if ((unary_fn(j, p_indices[i], p_elements + i, param)))
            {
                return JMTX_RESULT_UNARY_RETURN;
            }
        }
    }
    return JMTX_RESULT_SUCCESS;
}

void JMTX_NAME_TYPED(matrix_ccs_remove_zeros)(JMTX_NAME_TYPED(matrix_ccs) * mtx)
{
    //  Update offsets
    JMTX_INDEX_T p, c = 0, r = 0;
    while (mtx->end_of_column_offsets[r] == 0)
    {
        r += 1;
    }
    for (p = 0; p < mtx->n_entries; ++p)
    {
        assert(r < mtx->base.cols);
        if (mtx->values[p] == 0)
        {
            c += 1;
        }
        while (r < mtx->base.cols && p + 1 == mtx->end_of_column_offsets[r])
        {
            assert(mtx->end_of_column_offsets[r] >= c);
            mtx->end_of_column_offsets[r] -= c;
            r += 1;
        }
    }
    assert(r == mtx->base.cols);
    assert(c <= mtx->n_entries);
#ifndef NDEBUG
    for (r = 1; r < mtx->base.cols; ++r)
    {
        assert(mtx->end_of_column_offsets[r - 1] <= mtx->end_of_column_offsets[r]);
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

void JMTX_NAME_TYPED(matrix_ccs_remove_bellow)(JMTX_NAME_TYPED(matrix_ccs) * mtx, JMTX_REAL_T v)
{
    //  Update offsets
    JMTX_INDEX_T p, c = 0, r = 0;
    while (mtx->end_of_column_offsets[r] == 0)
    {
        r += 1;
    }
    for (p = 0; p < mtx->n_entries; ++p)
    {
        assert(r < mtx->base.cols);
        if (JMTX_ABS(mtx->values[p]) < v)
        {
            c += 1;
        }
        while (r < mtx->base.cols && p + 1 == mtx->end_of_column_offsets[r])
        {
            assert(mtx->end_of_column_offsets[r] >= c);
            mtx->end_of_column_offsets[r] -= c;
            r += 1;
        }
    }
    assert(r == mtx->base.rows);
    assert(c <= mtx->n_entries);
#ifndef NDEBUG
    for (r = 1; r < mtx->base.cols; ++r)
    {
        assert(mtx->end_of_column_offsets[r - 1] <= mtx->end_of_column_offsets[r]);
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

    //  Beef
    mtx->n_entries -= c;
}

JMTX_INDEX_T JMTX_NAME_TYPED(matrix_ccs_elements_in_row)(const JMTX_NAME_TYPED(matrix_ccs) * mtx, JMTX_INDEX_T row)
{
    JMTX_INDEX_T element_count = 0;
    for (JMTX_INDEX_T col = 0;
         col < mtx->base.cols && (!col || (mtx->end_of_column_offsets[col - 1] != mtx->n_entries)); ++col)
    {
        JMTX_INDEX_T *col_indices;
        JMTX_SCALAR_T *unused_col_values;
        const JMTX_INDEX_T n_col_elements = ccs_get_column_entries(mtx, col, &col_indices, &unused_col_values);
        if (n_col_elements && col_indices[0] <= row && col_indices[n_col_elements - 1] >= row)
        {
            JMTX_INDEX_T current = jmtx_internal_find_last_leq_value(n_col_elements, col_indices, row);
            if (col_indices[current] == row)
            {
                element_count += 1;
            }
        }
    }
    return element_count;
}

JMTX_INDEX_T JMTX_NAME_TYPED(matrix_ccs_get_row)(const JMTX_NAME_TYPED(matrix_ccs) * mtx, JMTX_INDEX_T row,
                                                 JMTX_INDEX_T n, JMTX_SCALAR_T *p_values, JMTX_INDEX_T *p_columns)
{
    JMTX_INDEX_T k = 0;
    for (JMTX_INDEX_T col = 0;
         k < n && col < mtx->base.cols && (!col || (mtx->end_of_column_offsets[col - 1] != mtx->n_entries)); ++col)
    {
        JMTX_INDEX_T *col_indices;
        JMTX_SCALAR_T *unused_col_values;
        const JMTX_INDEX_T n_col_elements = ccs_get_column_entries(mtx, col, &col_indices, &unused_col_values);
        if (n_col_elements && col_indices[0] <= row && col_indices[n_col_elements - 1] >= row)
        {
            const JMTX_INDEX_T current = jmtx_internal_find_last_leq_value(n_col_elements, col_indices, row);
            if (col_indices[current] == row)
            {
                p_values[k] = mtx->values[(col ? mtx->end_of_column_offsets[col - 1] : 0) + current];
                p_columns[k] = col;
                k += 1;
            }
        }
    }

    return k;
}

jmtx_result JMTX_NAME_TYPED(matrix_ccs_transpose)(const JMTX_NAME_TYPED(matrix_ccs) * mtx,
                                                  JMTX_NAME_TYPED(matrix_ccs) * *p_out,
                                                  const jmtx_allocator_callbacks *allocator_callbacks)
{
    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }

    const JMTX_INDEX_T rows = mtx->base.rows;
    JMTX_NAME_TYPED(matrix_ccs) * out;
    jmtx_result res =
        JMTX_NAME_TYPED(matrix_ccs_new)(&out, mtx->base.cols, mtx->base.rows, mtx->n_entries, allocator_callbacks);
    if (res != JMTX_RESULT_SUCCESS)
    {
        return res;
    }
    if (!allocator_callbacks)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }

    JMTX_INDEX_T *row_counts = allocator_callbacks->alloc(allocator_callbacks->state, sizeof *row_counts * rows);
    if (row_counts == NULL)
    {
        JMTX_NAME_TYPED(matrix_ccs_destroy)(out);
        return JMTX_RESULT_BAD_ALLOC;
    }
    memset(row_counts, 0, sizeof *row_counts * rows);

    JMTX_INDEX_T *col_ends = out->end_of_column_offsets;
    for (JMTX_INDEX_T i = 0; i < mtx->n_entries; ++i)
    {
        row_counts[mtx->indices[i]] += 1;
    }
    col_ends[0] = row_counts[0];
    //  Compute cumsums for offsets
    for (JMTX_INDEX_T i = 1; i < rows; ++i)
    {
        col_ends[i] = row_counts[i] + col_ends[i - 1];
        row_counts[i] = 0; //   Zero the row counts so that they can be reused later for counting bucket sizes
    }
    row_counts[0] = 0;
    row_counts[rows - 1] = 0;

    for (JMTX_INDEX_T col = 0; col < mtx->base.cols; ++col)
    {
        JMTX_INDEX_T *in_rows;
        JMTX_SCALAR_T *in_vals;
        JMTX_INDEX_T n_col = ccs_get_column_entries(mtx, col, &in_rows, &in_vals);

        for (JMTX_INDEX_T idx = 0; idx < n_col; ++idx)
        {
            const JMTX_INDEX_T row = in_rows[idx];
            const JMTX_INDEX_T ip = row > 0 ? col_ends[row - 1] : 0;
            const JMTX_INDEX_T n_rows = row_counts[row];

            out->values[ip + n_rows] = in_vals[idx];
            out->indices[ip + n_rows] = col;
            row_counts[row] += 1;
        }
    }
    out->n_entries = mtx->n_entries;

    allocator_callbacks->free(allocator_callbacks->state, row_counts);
    *p_out = out;

    return JMTX_RESULT_SUCCESS;
}

jmtx_result JMTX_NAME_TYPED(matrix_ccs_copy)(const JMTX_NAME_TYPED(matrix_ccs) * mtx,
                                             JMTX_NAME_TYPED(matrix_ccs) * *p_out,
                                             const jmtx_allocator_callbacks *allocator_callbacks)
{
    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }

    JMTX_NAME_TYPED(matrix_ccs) *const this = allocator_callbacks->alloc(allocator_callbacks->state, sizeof(*this));
    if (!this)
    {
        return JMTX_RESULT_BAD_ALLOC;
    }

    JMTX_SCALAR_T *const elements =
        mtx->base.allocator_callbacks.alloc(mtx->base.allocator_callbacks.state, (mtx->n_entries) * sizeof *elements);
    if (!elements)
    {
        allocator_callbacks->free(allocator_callbacks->state, this);
        return JMTX_RESULT_BAD_ALLOC;
    }
    JMTX_INDEX_T *const indices =
        mtx->base.allocator_callbacks.alloc(mtx->base.allocator_callbacks.state, (mtx->n_entries) * sizeof *indices);
    if (!indices)
    {
        mtx->base.allocator_callbacks.free(mtx->base.allocator_callbacks.state, elements);
        allocator_callbacks->free(allocator_callbacks->state, this);
        return JMTX_RESULT_BAD_ALLOC;
    }
    JMTX_INDEX_T *const cum_sum =
        mtx->base.allocator_callbacks.alloc(mtx->base.allocator_callbacks.state, (mtx->base.cols) * sizeof *cum_sum);
    if (!cum_sum)
    {
        mtx->base.allocator_callbacks.free(mtx->base.allocator_callbacks.state, indices);
        mtx->base.allocator_callbacks.free(mtx->base.allocator_callbacks.state, elements);
        allocator_callbacks->free(allocator_callbacks->state, this);
        return JMTX_RESULT_BAD_ALLOC;
    }

    memcpy(elements, mtx->values, sizeof *elements * mtx->n_entries);
    memcpy(indices, mtx->indices, sizeof *indices * mtx->n_entries);
    memcpy(cum_sum, mtx->end_of_column_offsets, sizeof *cum_sum * (mtx->base.cols));
    this->base = mtx->base;
    this->values = elements;
    this->indices = indices;
    this->end_of_column_offsets = cum_sum;
    this->capacity = mtx->n_entries;
    this->n_entries = this->n_entries;
    *p_out = this;
    return JMTX_RESULT_SUCCESS;
}

void JMTX_NAME_TYPED(matrix_zero_all_entries)(const JMTX_NAME_TYPED(matrix_ccs) * mtx)
{
    memset(mtx->values, 0, sizeof(*mtx->values) * mtx->n_entries);
}

void JMTX_NAME_TYPED(matrix_ccs_set_all_entries)(const JMTX_NAME_TYPED(matrix_ccs) * mtx, JMTX_SCALAR_T x)
{
    for (JMTX_SCALAR_T *ptr = mtx->values; ptr != mtx->values + mtx->n_entries; ++ptr)
    {
        *ptr = x;
    }
}

jmtx_result JMTX_NAME_TYPED(matrix_ccs_build_col)(JMTX_NAME_TYPED(matrix_ccs) * mtx, JMTX_INDEX_T col, JMTX_INDEX_T n,
                                                  const JMTX_INDEX_T *indices, const JMTX_SCALAR_T *values)
{
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

    const JMTX_INDEX_T offset = col ? mtx->end_of_column_offsets[col - 1] : 0;
    memcpy(mtx->values + offset, values, sizeof *values * n);
    memcpy(mtx->indices + offset, indices, sizeof *indices * n);

    mtx->end_of_column_offsets[col] = n + offset;
    mtx->n_entries += n;

    return JMTX_RESULT_SUCCESS;
}

/**
 * Finds the upper bandwidth of the matrix; what is the furthest distance of and entry above the main diagonal
 * @param mtx matrx to find the upper bandwidth of
 * @return upper bandwidth of the matrix
 */
JMTX_INDEX_T JMTX_NAME_TYPED(matrix_ccs_find_upper_bandwidth)(const JMTX_NAME_TYPED(matrix_ccs) * mtx)
{
    //  Find the greatest distance above the main diagonal
    JMTX_FAST_INT_T v_max = 0;
    for (JMTX_FAST_INT_T i = 0, p = 0; i < mtx->base.cols; ++i)
    {
        for (; p < mtx->end_of_column_offsets[i]; ++p)
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

/**
 * Finds the lower bandwidth of the matrix; what is the furthest distance of and entry below the main diagonal
 * @param mtx matrx to find the lower bandwidth of
 * @return lower bandwidth of the matrix
 */
JMTX_INDEX_T JMTX_NAME_TYPED(matrix_ccs_find_lower_bandwidth)(const JMTX_NAME_TYPED(matrix_ccs) * mtx)
{
    //  Find the greatest distance above the main diagonal
    JMTX_FAST_INT_T v_max = 0;
    for (JMTX_FAST_INT_T i = 0, p = 0; i < mtx->base.cols; ++i)
    {
        for (; p < mtx->end_of_column_offsets[i]; ++p)
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
