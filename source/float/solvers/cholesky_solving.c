//
// Created by jan on 6.11.2023.
//

#include <assert.h>
#include <math.h>
#include "../../../include/jmtx/float/solvers/cholesky_solving.h"
#include "../matrices/sparse_row_compressed_internal.h"
#include "../../../include/jmtx/float/solvers/incomplete_cholesky_decomposition.h"
#include "../../../include/jmtx/float/matrices/sparse_conversion.h"

/**
 * Solves a problem A x = C C^T x = y, where C is a lower triangular matrix.
 * @param c lower triangular matrix C in the CCS format
 * @param ct transpose of the matrix C in the CRS format
 * @param y memory containing forcing vector
 * @param x memory which receives the solution
 */
void jmtx_cholesky_solve(const jmtx_matrix_crs* c, const jmtx_matrix_crs* ct, const float* restrict y, float* restrict x)
{
    const uint32_t n = c->base.cols;
    x[0] = y[0];
    //  First is the forward substitution for C v = y
    for (uint32_t i = 1; i < n; ++i)
    {
        uint32_t* indices;
        float* values;
        uint32_t count = jmtx_matrix_crs_get_row(c, i, &indices, &values);
        assert(indices[count - 1] == (uint32_t)i);

        float v = 0;
        for (uint32_t j = 0; j < count - 1; ++j)
        {
            assert(indices[j] < i);
            v += values[j] * x[indices[j]];
        }
        x[i] = (y[i] - v) / values[count - 1];
    }
    //  Then the backward substitution for C^T x = v
    for (int32_t i = (int32_t)n - 1; i >= 0; --i)
    {
        uint32_t* indices;
        float* values;
        uint32_t count = jmtx_matrix_crs_get_row(ct, i, &indices, &values);
        assert(indices[0] == (uint32_t)i);

        float v = 0;
        for (uint32_t j = 1; j < count; ++j)
        {
            assert(indices[j] > (uint32_t)i);
            v += values[j] * x[indices[j]];
        }
        x[i] = (x[i] - v) / values[0];
    }
}

/**
 * Solves a problem A x = C C^T x = y, where C is a lower triangular matrix. This version of the function stores the
 * solution vector x back into the same memory where the forcing vector was.
 * @param c lower triangular matrix C in the CCS format
 * @param ct transpose of the matrix C in the CRS format
 * @param x memory which contains the forcing vector and receives the solution
 */
void jmtx_cholesky_solve_inplace(const jmtx_matrix_crs* c, const jmtx_matrix_crs* ct, float* restrict x)
{
    const uint32_t n = c->base.cols;
    //  First is the forward substitution for C v = y
    for (uint32_t i = 1; i < n; ++i)
    {
        uint32_t* indices;
        float* values;
        uint32_t count = jmtx_matrix_crs_get_row(c, i, &indices, &values);
        assert(indices[count - 1] == (uint32_t)i);

        float v = 0;
        for (uint32_t j = 0; j < count - 1; ++j)
        {
            assert(indices[j] < i);
            v += values[j] * x[indices[j]];
        }
        x[i] = (x[i] - v) / values[count - 1];
    }
    //  Then the backward substitution for C^T x = v
    for (int32_t i = (int32_t)n - 1; i >= 0; --i)
    {
        uint32_t* indices;
        float* values;
        uint32_t count = jmtx_matrix_crs_get_row(ct, i, &indices, &values);
        assert(indices[0] == (uint32_t)i);

        float v = 0;
        for (uint32_t j = 1; j < count; ++j)
        {
            assert(indices[j] > (uint32_t)i);
            v += values[j] * x[indices[j]];
        }
        x[i] = (x[i] - v) / values[0];
    }
}

static inline int check_vector_overlaps(const unsigned n, const size_t size, const void* ptrs[static const n])
{
    for (unsigned i = 0; i < n; ++i)
    {
        const uintptr_t p1 = (uintptr_t)ptrs[i];
        for (unsigned j = i + 1; j < n; ++j)
        {
            const uintptr_t p2 = (uintptr_t)ptrs[j];
            if (p1 > p2)
            {
                if (p2 + size > p1)
                {
                    return 1;
                }
            }
            else if (p1 + size > p2)
            {
                return 1;
            }
        }
    }

    return 0;
}
/**
 * Solves a problem A x = C C^T x = y, where C is a lower triangular matrix.
 * @param c lower triangular matrix C in the CCS format
 * @param ct transpose of the matrix C in the CRS format
 * @param y memory containing forcing vector
 * @param x memory which receives the solution
 * @returns JMTX_RESULT_SUCCESS if successful, otherwise an error code indicating error in the input parameters
 */
jmtx_result jmtxs_cholesky_solve(const jmtx_matrix_crs* c, const jmtx_matrix_crs* ct, uint32_t n,
                                 const float y[static restrict n], float x[restrict n])
{
    if (!c)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (c->base.type != JMTX_TYPE_CRS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (c->base.rows != n || c->base.cols != n)
    {
        return JMTX_RESULT_BAD_MATRIX;
    }

    if (!ct)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (ct->base.type != JMTX_TYPE_CRS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (ct->base.rows != n || ct->base.cols != n)
    {
        return JMTX_RESULT_BAD_MATRIX;
    }

    if (!x)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!y)
    {
        return JMTX_RESULT_NULL_PARAM;
    }

    const void* ptrs[] = {x, y};
    if (check_vector_overlaps(sizeof(ptrs) / sizeof(*ptrs), sizeof(*x) * n, ptrs))
    {
        return JMTX_RESULT_BAD_PARAM;
    }
    jmtx_cholesky_solve(c, ct, y, x);
    return JMTX_RESULT_SUCCESS;
}

/**
 * Solves a problem A x = C C^T x = y, where C is a lower triangular matrix. This version of the function stores the
 * solution vector x back into the same memory where the forcing vector was.
 * @param c lower triangular matrix C in the CCS format
 * @param ct transpose of the matrix C in the CRS format
 * @param x memory which contains the forcing vector and receives the solution
 * @returns JMTX_RESULT_SUCCESS if successful, otherwise an error code indicating error in the input parameters
 */
jmtx_result jmtxs_cholesky_solve_inplace(const jmtx_matrix_crs* c, const jmtx_matrix_crs* ct, uint32_t n,
                                         float x[static n])
{
    if (!c)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (c->base.type != JMTX_TYPE_CRS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (c->base.rows != n || c->base.cols != n)
    {
        return JMTX_RESULT_BAD_MATRIX;
    }

    if (!ct)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (ct->base.type != JMTX_TYPE_CRS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (ct->base.rows != n || ct->base.cols != n)
    {
        return JMTX_RESULT_BAD_MATRIX;
    }

    if (!x)
    {
        return JMTX_RESULT_NULL_PARAM;
    }

    jmtx_cholesky_solve_inplace(c, ct, x);
    return JMTX_RESULT_SUCCESS;
}
