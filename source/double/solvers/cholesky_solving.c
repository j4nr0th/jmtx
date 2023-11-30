// Automatically generated from source/float/solvers/cholesky_solving.c on Thu Nov 30 20:33:08 2023
//
// Created by jan on 6.11.2023.
//

#include <assert.h>
#include <math.h>
#include "../../../include/jmtx/double/solvers/cholesky_solving.h"
#include "../matrices/sparse_row_compressed_internal.h"
#include "../../../include/jmtx/double/solvers/incomplete_cholesky_decomposition.h"
#include "../../../include/jmtx/double/matrices/sparse_conversion.h"

void jmtxd_cholesky_solve(const jmtxd_matrix_crs* c, const jmtxd_matrix_crs* ct, const double* restrict y, double* restrict x)
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

void jmtxd_cholesky_solve_inplace(const jmtxd_matrix_crs* c, const jmtxd_matrix_crs* ct, double* restrict x)
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
