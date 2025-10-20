#include "cholesky_solving.h"
#include "../decompositions/incomplete_cholesky_decomposition.h"
#include "../matrices/sparse_conversion.h"
#include "../matrices/sparse_row_compressed.h"
#include <assert.h>
#include <math.h>

/**
 * Solves a problem A x = C C^T x = y, where C is a lower triangular matrix.
 * @param c lower triangular matrix C in the CCS format
 * @param ct transpose of the matrix C in the CRS format
 * @param y memory containing forcing vector
 * @param x memory which receives the solution
 */
void JMTX_NAME_TYPED(solve_direct_cholesky_crs)(const JMTX_NAME_TYPED(matrix_crs) * c,
                                                const JMTX_NAME_TYPED(matrix_crs) * ct, const JMTX_SCALAR_T *restrict y,
                                                JMTX_SCALAR_T *restrict x)
{
    const JMTX_INDEX_T n = c->base.cols;
    x[0] = y[0];
    //  First is the forward substitution for C v = y
    for (JMTX_INDEX_T i = 1; i < n; ++i)
    {
        JMTX_INDEX_T *indices;
        JMTX_SCALAR_T *values;
        JMTX_INDEX_T count = JMTX_NAME_TYPED(matrix_crs_get_row)(c, i, &indices, &values);
        assert(indices[count - 1] == (JMTX_INDEX_T)i);

        JMTX_SCALAR_T v = 0;
        for (JMTX_INDEX_T j = 0; j < count - 1; ++j)
        {
            assert(indices[j] < i);
            v += values[j] * x[indices[j]];
        }
        x[i] = (y[i] - v) / values[count - 1];
    }
    //  Then the backward substitution for C^T x = v
    for (int32_t i = (int32_t)n - 1; i >= 0; --i)
    {
        JMTX_INDEX_T *indices;
        JMTX_SCALAR_T *values;
        const JMTX_INDEX_T count = JMTX_NAME_TYPED(matrix_crs_get_row)(ct, i, &indices, &values);
        assert(indices[0] == (JMTX_INDEX_T)i);

        JMTX_SCALAR_T v = 0;
        for (JMTX_INDEX_T j = 1; j < count; ++j)
        {
            assert(indices[j] > (JMTX_INDEX_T)i);
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
void JMTX_NAME_TYPED(solve_direct_cholesky_crs_inplace)(const JMTX_NAME_TYPED(matrix_crs) * c,
                                                        const JMTX_NAME_TYPED(matrix_crs) * ct,
                                                        JMTX_SCALAR_T *restrict x)
{
    const JMTX_INDEX_T n = c->base.cols;
    //  First is the forward substitution for C v = y
    for (JMTX_INDEX_T i = 1; i < n; ++i)
    {
        JMTX_INDEX_T *indices;
        JMTX_SCALAR_T *values;
        const JMTX_INDEX_T count = JMTX_NAME_TYPED(matrix_crs_get_row)(c, i, &indices, &values);
        assert(indices[count - 1] == (JMTX_INDEX_T)i);

        JMTX_SCALAR_T v = 0;
        for (JMTX_INDEX_T j = 0; j < count - 1; ++j)
        {
            assert(indices[j] < i);
            v += values[j] * x[indices[j]];
        }
        x[i] = (x[i] - v) / values[count - 1];
    }
    //  Then the backward substitution for C^T x = v
    for (int32_t i = (int32_t)n - 1; i >= 0; --i)
    {
        JMTX_INDEX_T *indices;
        JMTX_SCALAR_T *values;
        JMTX_INDEX_T count = JMTX_NAME_TYPED(matrix_crs_get_row)(ct, i, &indices, &values);
        assert(indices[0] == (JMTX_INDEX_T)i);

        JMTX_SCALAR_T v = 0;
        for (JMTX_INDEX_T j = 1; j < count; ++j)
        {
            assert(indices[j] > (JMTX_INDEX_T)i);
            v += values[j] * x[indices[j]];
        }
        x[i] = (x[i] - v) / values[0];
    }
}
