// Automatically generated from source/float/solvers/cholesky_solving.c on Sun Dec 17 20:12:49 2023
//
// Created by jan on 6.11.2023.
//

#include <assert.h>
#include <math.h>
#include "../matrices/sparse_row_compressed_internal.h"
#include "../../../include/jmtx/double/decompositions/incomplete_cholesky_decomposition.h"
#include "../../../include/jmtx/double/matrices/sparse_conversion.h"
#include "../../../include/jmtx/double/solvers/cholesky_solving.h"

/**
 * Solves a problem A x = C C^T x = y, where C is a lower triangular matrix.
 * @param c lower triangular matrix C in the CCS format
 * @param ct transpose of the matrix C in the CRS format
 * @param y memory containing forcing vector
 * @param x memory which receives the solution
 */
void jmtxd_solve_direct_cholesky_crs(const jmtxd_matrix_crs* c, const jmtxd_matrix_crs* ct, const double* restrict y, double* restrict x)
{
    const uint32_t n = c->base.cols;
    x[0] = y[0];
    //  First is the forward substitution for C v = y
    for (uint32_t i = 1; i < n; ++i)
    {
        uint32_t* indices;
        double* values;
        uint32_t count = jmtxd_matrix_crs_get_row(c, i, &indices, &values);
        assert(indices[count - 1] == (uint32_t)i);

        double v = 0;
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
        double* values;
        uint32_t count = jmtxd_matrix_crs_get_row(ct, i, &indices, &values);
        assert(indices[0] == (uint32_t)i);

        double v = 0;
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
void jmtxd_solve_direct_cholesky_crs_inplace(const jmtxd_matrix_crs* c, const jmtxd_matrix_crs* ct, double* restrict x)
{
    const uint32_t n = c->base.cols;
    //  First is the forward substitution for C v = y
    for (uint32_t i = 1; i < n; ++i)
    {
        uint32_t* indices;
        double* values;
        uint32_t count = jmtxd_matrix_crs_get_row(c, i, &indices, &values);
        assert(indices[count - 1] == (uint32_t)i);

        double v = 0;
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
        double* values;
        uint32_t count = jmtxd_matrix_crs_get_row(ct, i, &indices, &values);
        assert(indices[0] == (uint32_t)i);

        double v = 0;
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
jmtx_result jmtxds_solve_direct_cholesky_crs(const jmtxd_matrix_crs* c, const jmtxd_matrix_crs* ct, uint32_t n,
                                 const double y[static restrict n], double x[restrict n])
{
    if (!c)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (c->base.type != JMTXD_TYPE_CRS)
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
    if (ct->base.type != JMTXD_TYPE_CRS)
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
    jmtxd_solve_direct_cholesky_crs(c, ct, y, x);
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
jmtx_result jmtxds_solve_direct_cholesky_crs_inplace(const jmtxd_matrix_crs* c, const jmtxd_matrix_crs* ct, uint32_t n,
                                         double x[static n])
{
    if (!c)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (c->base.type != JMTXD_TYPE_CRS)
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
    if (ct->base.type != JMTXD_TYPE_CRS)
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

    jmtxd_solve_direct_cholesky_crs_inplace(c, ct, x);
    return JMTX_RESULT_SUCCESS;
}
