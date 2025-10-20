//
// Created by jan on 1.1.2024.
//

#ifndef JMTX_GMRES_INTERNAL_H
#define JMTX_GMRES_INTERNAL_H

#include "../matrices/band_row_major.h"
#include "../matrices/sparse_diagonal_compressed.h"
#include "../matrices/sparse_row_compressed.h"

JMTX_INDEX_T JMTX_NAME_TYPED(gmresm_round_cds)(
    const JMTX_NAME_TYPED(matrix_cds) * mtx, const JMTX_INDEX_T n, const JMTX_INDEX_T m, const JMTX_REAL_T y_mag,
    const JMTX_REAL_T tol, const JMTX_SCALAR_T residual[JMTX_ARRAY_ATTRIB(const restrict static n)],
    JMTX_SCALAR_T x[JMTX_ARRAY_ATTRIB(const restrict static n)], JMTX_NAME_TYPED(matrix_brm) * r,
    JMTX_SCALAR_T ck[JMTX_ARRAY_ATTRIB(const restrict m)], JMTX_SCALAR_T sk[JMTX_ARRAY_ATTRIB(const restrict m)],
    JMTX_SCALAR_T g[JMTX_ARRAY_ATTRIB(const restrict m)], JMTX_SCALAR_T alpha[JMTX_ARRAY_ATTRIB(const restrict m)],
    JMTX_SCALAR_T h[JMTX_ARRAY_ATTRIB(const restrict m)], JMTX_SCALAR_T p_mat[JMTX_ARRAY_ATTRIB(const restrict m *n)]);

JMTX_INDEX_T JMTX_NAME_TYPED(gmresm_round_crs)(
    const JMTX_NAME_TYPED(matrix_crs) * mtx, const JMTX_INDEX_T n, const JMTX_INDEX_T m, const JMTX_REAL_T y_mag,
    const JMTX_REAL_T tol, const JMTX_SCALAR_T residual[JMTX_ARRAY_ATTRIB(const restrict static n)],
    JMTX_SCALAR_T x[JMTX_ARRAY_ATTRIB(const restrict static n)], JMTX_NAME_TYPED(matrix_brm) * r,
    JMTX_SCALAR_T ck[JMTX_ARRAY_ATTRIB(const restrict m)], JMTX_SCALAR_T sk[JMTX_ARRAY_ATTRIB(const restrict m)],
    JMTX_SCALAR_T g[JMTX_ARRAY_ATTRIB(const restrict m)], JMTX_SCALAR_T alpha[JMTX_ARRAY_ATTRIB(const restrict m)],
    JMTX_SCALAR_T h[JMTX_ARRAY_ATTRIB(const restrict m)], JMTX_SCALAR_T p_mat[JMTX_ARRAY_ATTRIB(const restrict m *n)]);

#endif // JMTX_GMRES_INTERNAL_H
