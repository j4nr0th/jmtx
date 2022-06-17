//
// Created by jan on 13.6.2022.
//

#include "sparse_row_compressed.h"
#define DEFAULT_RESERVED_ELEMENTS 64

#define FILL_VALUE 0xDEADBEEF

static void beef_it_up(scalar_t* ptr, size_t elements)
{
    if (elements < 2)
    {
        const uint32_t BeefBuffer2[2] = {0xDEADBEEF, 0xDEADBEEF};
        memcpy(ptr, BeefBuffer2, elements * sizeof(*ptr));
    }

    while (elements & 0x3)
    {
        *(uint32_t*)ptr = 0xDEADBEEF;
        --elements;
        ++ptr;
    }

    for (uint i = 0; i < elements / 4; ++i)
    {
        *((uint32_t*)(ptr + i)) = 0xDEADBEEF;
    }
    memcpy(ptr + elements / 4, ptr, elements / 4 * sizeof(*ptr));
    memcpy(ptr + elements / 2, ptr, elements / 2 * sizeof(*ptr));
    //  Beefed
}

static int beef_check(const scalar_t* ptr, size_t elements)
{
    const uint32_t* const buffer = (const uint32_t*)ptr;
    int beef_count = 0;
    for (uint i = 0; i < elements; ++i)
    {
        beef_count += (buffer[i] == 0xDEADBEEF);
    }
    return (beef_count << 16) | 0x0000BEEF;
}

mtx_res_t matrix_crs_new(CrsMatrix* mtx, uint columns, uint rows, uint reserved_elements)
{
    if (reserved_elements == 0)
    {
        reserved_elements = DEFAULT_RESERVED_ELEMENTS;
        reserved_elements = reserved_elements < columns * rows ? reserved_elements : columns * rows;
    }
    mtx_res_t res = 0;
    uint* elements_per_row = NULL;
    uint* indices = NULL;
    scalar_t* p_elements = calloc(1 + reserved_elements, sizeof*p_elements);
    if (p_elements == NULL)
    {
        res = -1;
        goto fail1;
    }

    indices = calloc(1 + reserved_elements, sizeof*indices);
    if (indices == NULL)
    {
        res = -1;
        goto fail2;
    }

    elements_per_row = calloc(rows + 1, sizeof*elements_per_row);
    if (elements_per_row == NULL)
    {
        res = -1;
        goto fail3;
    }
    beef_it_up(p_elements + 1, reserved_elements);
    _Static_assert(sizeof(scalar_t) == sizeof(uint), "element and index sizes must be the same");
    beef_it_up((scalar_t*)indices + 1, reserved_elements);
    for (uint i = 0; i < rows + 1; ++i)
    {
        elements_per_row[i] = 1;
    }
    memset(mtx, 0, sizeof*mtx);
    mtx->indices = indices;
    mtx->columns = columns;
    mtx->type = mtx_type_crs;
    mtx->elements = p_elements;
    mtx->rows = rows;
    mtx->capacity = reserved_elements;
    mtx->n_elements = 0;
    mtx->elements_before = elements_per_row;

    return res;
    fail3: free(indices);
    fail2: free(elements_per_row);
    fail1: free(p_elements);
    return res;
}

mtx_res_t matrix_crs_destroy(CrsMatrix* mtx)
{
    mtx_res_t res = 0;
    free(mtx->indices);
    free(mtx->elements_before);
    free(mtx->elements);
    return res;
}

mtx_res_t matrix_crs_shrink(CrsMatrix* mtx)
{
    mtx_res_t res = 0;

    if (mtx->n_elements == mtx->capacity) return 0;
    mtx->capacity = mtx->n_elements;
    scalar_t* element_new_ptr = realloc(mtx->elements, sizeof*mtx->elements * (mtx->n_elements + 1));
    if (!element_new_ptr)
    {
        res = -1;
        goto end;
    }
    mtx->elements = element_new_ptr;
    uint* new_indices_ptr = realloc(mtx->indices, sizeof*mtx->indices * (mtx->n_elements + 1));
    if (!new_indices_ptr)
    {
        res = -1;
        goto end;
    }
    mtx->indices = new_indices_ptr;

end:
    return res;
}

mtx_res_t matrix_crs_set_row(CrsMatrix* mtx, uint row, uint n, const uint* indices, const scalar_t* elements)
{
    mtx_res_t res = 0;
    const int new_elements = (int)n - (int)(mtx->elements_before[row + 1] - mtx->elements_before[row]);
    const int required_capacity = (int)mtx->n_elements + new_elements;
    if (mtx->capacity < required_capacity)
    {
        scalar_t* new_element_ptr = realloc(mtx->elements, sizeof*mtx->elements * (required_capacity + 1));
        if (!new_element_ptr)
        {
            res = -1;
            goto end;
        }
        mtx->elements = new_element_ptr;
        uint* new_indices_ptr = realloc(mtx->indices, sizeof*mtx->indices * (required_capacity + 1));
        if (!new_indices_ptr)
        {
            res = -1;
            goto end;
        }
        mtx->indices = new_indices_ptr;
    }

    const uint elements_after = mtx->elements_before[mtx->rows] - mtx->elements_before[row + 1];
    if (elements_after)
    {
        memmove(mtx->elements + mtx->elements_before[row + 1] + new_elements, mtx->elements + mtx->elements_before[row + 1],
                sizeof*mtx->elements * (elements_after));
        memmove(mtx->indices + mtx->elements_before[row + 1] + new_elements, mtx->indices + mtx->elements_before[row + 1],
            sizeof*mtx->indices * (elements_after));
    }
    memcpy(mtx->elements + mtx->elements_before[row], elements, sizeof*elements * n);
    memcpy(mtx->indices + mtx->elements_before[row], indices, sizeof*indices * n);

    for (uint i = row; i < mtx->rows; ++i)
    {
        mtx->elements_before[i + 1] += new_elements;
    }
    mtx->n_elements += new_elements;
end:
    return res;
}

mtx_res_t matrix_crs_vector_multiply(const CrsMatrix* mtx, const scalar_t* restrict x, scalar_t* restrict y)
{
    mtx_res_t res = 0;

    for (uint i = 0; i < mtx->rows; ++i)
    {
        const uint* indices = mtx->indices + mtx->elements_before[i];
        const scalar_t* row_ptr = mtx->elements + mtx->elements_before[i];
        scalar_t v = 0;
        const uint n_elements = mtx->elements_before[i + 1] - mtx->elements_before[i];
        for (uint j = 0; j < n_elements; ++j)
        {
            const uint k = indices[j];
            v += row_ptr[j] * x[k];
        }
        y[i] = v;
    }

    return res;
}

mtx_res_t matrix_crs_set_element(CrsMatrix* mtx, uint i, uint j, scalar_t x)
{
    mtx_res_t res = 0;

    const uint n_row_elements = mtx->elements_before[i + 1] - mtx->elements_before[i];
    const uint* const row_indices = mtx->indices + mtx->elements_before[i];
    uint current = 0;
    for (
            uint size = n_row_elements, step = n_row_elements / 2;
            step != 0;
            step = size / 2
            )
    {
        if (row_indices[current + step] < j)
        {
            current += step;
            size -= step;
        }
        else if (row_indices[current + step] == j)
        {
            current += step;
            break;
        }
        else
        {
            size = step;
        }
    }

    if (n_row_elements && row_indices[current] == j)
    {
        *(mtx->elements + mtx->elements_before[i] + current) = x;
    }
    else
    {
        if (mtx->capacity == mtx->n_elements)
        {
            const uint new_capacity = mtx->capacity + DEFAULT_RESERVED_ELEMENTS;
            scalar_t* const new_element_ptr = realloc(mtx->elements, sizeof*mtx->elements * (new_capacity + 1));
            if (!new_element_ptr)
            {
                res = -1;
                goto end;
            }
            mtx->elements = new_element_ptr;
            beef_it_up(new_element_ptr + 1 + mtx->n_elements, new_capacity - mtx->capacity);
            uint* const new_index_ptr = realloc(mtx->indices, sizeof*mtx->indices * (new_capacity + 1));
            if (!new_index_ptr)
            {
                res = -1;
                goto end;
            }
            mtx->indices = new_index_ptr;
            beef_it_up((scalar_t *)(new_index_ptr + 1 + mtx->n_elements), new_capacity - mtx->capacity);
            mtx->capacity = new_capacity;
        }
        if (mtx->elements_before[mtx->rows] - mtx->elements_before[i + 1] != 0)
        {
            current += (n_row_elements != 0);
            memmove(mtx->elements + mtx->elements_before[i] + current + 1, mtx->elements + mtx->elements_before[i] + current, (mtx->elements_before[mtx->rows] - mtx->elements_before[i + 1] + n_row_elements - current) * sizeof(*mtx->elements));
            memmove(mtx->indices + mtx->elements_before[i] + current + 1, mtx->indices + mtx->elements_before[i] + current, (mtx->elements_before[mtx->rows] - mtx->elements_before[i + 1] + n_row_elements - current) * sizeof(*mtx->indices));
        }
        *(mtx->elements + mtx->elements_before[i] + current) = x;
        *(mtx->indices + mtx->elements_before[i] + current) = j;
        mtx->n_elements += 1;
        for (uint k = i; k < mtx->rows; ++k)
        {
            mtx->elements_before[k + 1] += 1;
        }
    }
end:
    return res;
}

mtx_res_t matrix_crs_get_element(const CrsMatrix* mtx, uint i, uint j, scalar_t* x)
{
    mtx_res_t res = 0;
    const uint n_row_elements = mtx->elements_before[i + 1] - mtx->elements_before[i];
    const uint* const row_indices = mtx->indices + mtx->elements_before[i];
    uint index = 0;
    uint current = 0;
    for (
            uint size = n_row_elements, step = n_row_elements / 2;
            step != 0;
            step = size / 2
        )
    {
        if (row_indices[current + step] < j)
        {
            current += step;
            size -= step;
        }
        else if (row_indices[current + step] == j)
        {
            current += step;
            break;
        }
        else
        {
            size = step;
        }
    }

    if (n_row_elements && row_indices[current] == j)
    {
        *x = *(mtx->elements + mtx->elements_before[i] + current);
    }
    else
    {
        *x = (scalar_t)0.0;
    }

    return res;
}

mtx_res_t matrix_crs_get_row(const CrsMatrix* mtx, uint row, uint* n, uint** p_indices, scalar_t** p_elements)
{
    mtx_res_t res = 0;
    const uint n_row = mtx->elements_before[row + 1] - mtx->elements_before[row];
    *n = n_row;
    if (n_row)
    {
        *p_indices = mtx->indices + mtx->elements_before[row];
        *p_elements = mtx->elements + mtx->elements_before[row];
    }

    return res;
}

mtx_res_t matrix_crs_beef_check(const CrsMatrix* mtx, int* p_beef_status)
{
    const int beef_status = beef_check(mtx->elements + 1, mtx->n_elements);
    *p_beef_status = beef_status;
    return 0;
}

mtx_res_t matrix_crs_apply_unary_fn(const CrsMatrix* mtx, int (*unary_fn)(uint i, uint j, scalar_t* p_element, void* param), void* param)
{
    for (uint i = 0; i < mtx->rows; ++i)
    {
        const uint n_in_row = mtx->elements_before[i + 1] - mtx->elements_before[i];
        scalar_t* const p_elements = mtx->elements + mtx->elements_before[i];
        const uint* const p_indices = mtx->indices + mtx->elements_before[i];
        for (uint j = 0; j < n_in_row; ++j)
        {
            int res;
            if ((res = unary_fn(i, p_indices[j], p_elements + j, param)))
            {
                return res;
            }
        }
    }
    return 0;
}

mtx_res_t matrix_crs_remove_zeros(CrsMatrix* mtx)
{
    uint zero_count, k, l;
    for (zero_count = 0, k = 0; k < mtx->n_elements; ++k)
    {
        zero_count += (mtx->elements[k + 1] == (scalar_t)0.0);
    }
    if (!zero_count) return 0;

    uint* const zero_indices = calloc(zero_count, sizeof*zero_indices);
    if (!zero_indices) return -1;

    for (k = 0, l = 0; l < zero_count && k < mtx->n_elements; ++k)
    {
        if (mtx->elements[k + 1] == (scalar_t)0.0)
        {
            zero_indices[l++] = k;
        }
    }

    const uint original_zero_count = zero_count;
    while (zero_count)
    {
        uint consecutive_zeros;
        for (consecutive_zeros = 1; consecutive_zeros < zero_count; ++consecutive_zeros)
        {
            if (zero_indices[zero_count - consecutive_zeros] != zero_indices[zero_count - consecutive_zeros - 1])
                break;
        }

        memmove(mtx->elements + 1 + zero_indices[zero_count - consecutive_zeros],
                mtx->elements + 1 + zero_indices[zero_count - 1] + 1,
                (mtx->n_elements - zero_indices[zero_count - 1] - 1) * sizeof*mtx->elements);
        memmove(mtx->indices + 1 + zero_indices[zero_count - consecutive_zeros],
                mtx->indices + 1 + zero_indices[zero_count - 1] + 1,
                (mtx->n_elements - zero_indices[zero_count - 1] - 1) * sizeof*mtx->indices);

        mtx->n_elements -= consecutive_zeros;
        zero_count -= consecutive_zeros;
    }

    for (k = 0, l = 0, zero_count = 0; k < mtx->rows + 1; ++k)
    {
        while (l < original_zero_count && zero_indices[l] + 1 < mtx->elements_before[k])
        {
            zero_count += 1;
            l += 1;
        }
        mtx->elements_before[k] -= zero_count;
    }

    //  Beef
    beef_it_up(mtx->elements + 1 + mtx->n_elements, original_zero_count);
    _Static_assert(sizeof(scalar_t) == sizeof(uint), "Size of index and scalar must be the same for beef");
    beef_it_up((scalar_t*)(mtx->indices + 1 + mtx->n_elements), original_zero_count);


    free(zero_indices);
    return 0;
}

mtx_res_t matrix_crs_remove_bellow(CrsMatrix* mtx, scalar_t v)
{
    //  AND-ing with 0x7FFFFFFF should remove the sign bit, allowing to use fast comparison of float's absolute values
    const uint int_abs_v = ((*(uint*)&v) & 0x7FFFFFFF);
    const scalar_t abs_v = *(scalar_t*)&int_abs_v;
    _Static_assert(sizeof(scalar_t) == sizeof(uint), "Size of scalar and uint must be the same for this to work");
    uint zero_count, k, l;
    for (zero_count = 0, k = 0; k < mtx->n_elements; ++k)
    {
        const uint element_abs = ((uint*)mtx->elements)[k + 1] & 0x7FFFFFFF;
        zero_count += (*(scalar_t*)&element_abs < abs_v);
    }
    if (!zero_count) return 0;

    uint* const zero_indices = calloc(zero_count, sizeof*zero_indices);
    if (!zero_indices) return -1;

    for (k = 0, l = 0; l < zero_count && k < mtx->n_elements; ++k)
    {
        const uint element_abs = ((uint*)mtx->elements)[k + 1] & 0x7FFFFFFF;
        if (*(scalar_t*)&element_abs < abs_v)
        {
            zero_indices[l++] = k;
        }
    }

    const uint original_zero_count = zero_count;
    while (zero_count)
    {
        uint consecutive_zeros;
        for (consecutive_zeros = 1; consecutive_zeros < zero_count; ++consecutive_zeros)
        {
            if (zero_indices[zero_count - consecutive_zeros] != zero_indices[zero_count - consecutive_zeros - 1])
                break;
        }

        memmove(mtx->elements + 1 + zero_indices[zero_count - consecutive_zeros],
                mtx->elements + 1 + zero_indices[zero_count - 1] + 1,
                (mtx->n_elements - zero_indices[zero_count - 1] - 1) * sizeof*mtx->elements);
        memmove(mtx->indices + 1 + zero_indices[zero_count - consecutive_zeros],
                mtx->indices + 1 + zero_indices[zero_count - 1] + 1,
                (mtx->n_elements - zero_indices[zero_count - 1] - 1) * sizeof*mtx->indices);

        mtx->n_elements -= consecutive_zeros;
        zero_count -= consecutive_zeros;
    }

    for (k = 0, l = 0, zero_count = 0; k < mtx->rows + 1; ++k)
    {
        while (l < original_zero_count && zero_indices[l] + 1 < mtx->elements_before[k])
        {
            zero_count += 1;
            l += 1;
        }
        mtx->elements_before[k] -= zero_count;
    }

    //  Beef
    beef_it_up(mtx->elements + 1 + mtx->n_elements, original_zero_count);
    _Static_assert(sizeof(scalar_t) == sizeof(uint), "Size of index and scalar must be the same for beef");
    beef_it_up((scalar_t*)(mtx->indices + 1 + mtx->n_elements), original_zero_count);


    free(zero_indices);
    return 0;
}

mtx_res_t matrix_crs_elements_in_column(const CrsMatrix* mtx, uint col, uint* p_n)
{
    uint element_count = 0;
    for (uint row = 0; row < mtx->rows && mtx->elements_before[row] != mtx->elements_before[mtx->rows]; ++row)
    {
        uint current = 0;
        const uint n_row_elements = mtx->elements_before[row + 1] -mtx->elements_before[row];
        const uint* row_indices = mtx->indices + mtx->elements_before[row];
        for (
                uint size = n_row_elements, step = n_row_elements / 2;
                step != 0;
                step = size / 2
                )
        {
            if (row_indices[current + step] < col)
            {
                current += step;
                size -= step;
            }
            else if (row_indices[current + step] == col)
            {
                current += step;
                break;
            }
            else
            {
                size = step;
            }
        }

        if (n_row_elements && row_indices[current] == col)
        {
            element_count += 1;
        }
    }
    *p_n = element_count;

    return 0;
}

mtx_res_t matrix_crs_get_column(const CrsMatrix* mtx, uint col, uint n, scalar_t* p_elements, uint* p_rows)
{
    uint k = 0;
    for (uint row = 0; row < mtx->rows && mtx->elements_before[row] != mtx->elements_before[mtx->rows] && k != n; ++row)
    {
        uint current = 0;
        const uint n_row_elements = mtx->elements_before[row + 1] -mtx->elements_before[row];
        const uint* row_indices = mtx->indices + mtx->elements_before[row];
        for (
                uint size = n_row_elements, step = n_row_elements / 2;
                step != 0;
                step = size / 2
                )
        {
            if (row_indices[current + step] < col)
            {
                current += step;
                size -= step;
            }
            else if (row_indices[current + step] == col)
            {
                current += step;
                break;
            }
            else
            {
                size = step;
            }
        }

        if (n_row_elements && row_indices[current] == col)
        {
            p_elements[k] = mtx->elements[mtx->elements_before[row] + current];
            p_rows[k] = row;
            k += 1;
        }
    }

    return 0;
}

mtx_res_t matrix_crs_transpose(const CrsMatrix* restrict mtx, CrsMatrix* restrict out)
{
    const uint n_elements = mtx->n_elements;
    const uint new_rows = mtx->columns;
    const uint new_cols = mtx->rows;

    uint* const column_cum_counts = calloc(new_rows + 1, sizeof*column_cum_counts);
    if (!column_cum_counts) return -1;
    uint* const new_indices = calloc(n_elements + 1, sizeof*new_indices);
    if (!new_indices)
    {
        free(column_cum_counts);
        return -1;
    }
    scalar_t* const new_elements = calloc(n_elements + 1, sizeof*new_elements);
    if (!new_elements)
    {
        free(column_cum_counts);
        free(new_indices);
        return -1;
    }

    *column_cum_counts = 1;

    for (uint j = 0, n, p = 1; j < mtx->columns; ++j)
    {
        matrix_crs_elements_in_column(mtx, j, &n);

        matrix_crs_get_column(mtx, j, n, new_elements + p, new_indices + p);
        p += n;
        column_cum_counts[j + 1] = column_cum_counts[j] + n;
    }

    memcpy(out, mtx, sizeof*out);
    out->elements_before = column_cum_counts;
    out->elements = new_elements;
    out->indices = new_indices;
    out->capacity = n_elements;
    out->rows = new_rows;
    out->columns = new_cols;

    return 0;
}

mtx_res_t matrix_crs_copy(const CrsMatrix* mtx, CrsMatrix* out)
{
    scalar_t* const elements = calloc(1 + mtx->n_elements, sizeof *elements);
    if (!elements) return -1;
    uint* const indices = calloc(1 + mtx->n_elements, sizeof *indices);
    if (!indices)
    {
        free(indices);
        return -1;
    }
    uint* const cum_sum = calloc(1 + mtx->rows, sizeof *cum_sum);
    if (!cum_sum)
    {
        free(indices);
        free(cum_sum);
        return -1;
    }

    memcpy(elements + 1, mtx->elements + 1, sizeof* elements * mtx->n_elements);
    memcpy(indices + 1, mtx->indices + 1, sizeof* indices * mtx->n_elements);
    memcpy(cum_sum, mtx->elements_before, sizeof* cum_sum * (mtx->rows + 1));
    memcpy(out, mtx, sizeof *out);
    out->elements = elements;
    out->indices = indices;
    out->elements_before = cum_sum;

    return 0;
}

mtx_res_t matrix_crs_build_row(CrsMatrix* mtx, uint row, uint n, const uint* indices, const scalar_t* elements)
{
    mtx_res_t res = 0;
    const int required_capacity = (int)mtx->n_elements + (int)n;
    if (mtx->capacity < required_capacity)
    {
        scalar_t* new_element_ptr = realloc(mtx->elements, sizeof*mtx->elements * (required_capacity + 1));
        if (!new_element_ptr)
        {
            res = -1;
            goto end;
        }
        mtx->elements = new_element_ptr;
        uint* new_indices_ptr = realloc(mtx->indices, sizeof*mtx->indices * (required_capacity + 1));
        if (!new_indices_ptr)
        {
            res = -1;
            goto end;
        }
        mtx->indices = new_indices_ptr;
    }

    memcpy(mtx->elements + mtx->elements_before[row], elements, sizeof*elements * n);
    memcpy(mtx->indices + mtx->elements_before[row], indices, sizeof*indices * n);

    mtx->elements_before[row + 1] = n + mtx->elements_before[row];
    mtx->n_elements += n;
    end:
    return res;
}
