//
// Created by jan on 1.1.2024.
//

#ifndef JMTX_GMRES_INTERNAL_H
#define JMTX_GMRES_INTERNAL_H

#include "../../../include/jmtx/float/matrices/band_row_major.h"
#include "../matrices/sparse_diagonal_compressed_internal.h"

uint32_t jmtx_gmresm_round_cds(const jmtx_matrix_cds *mtx, const uint32_t n, const uint32_t m, const float y_mag,
                               const float tol, const float residual[JMTX_ARRAY_ATTRIB(const restrict static n)],
                               float x[JMTX_ARRAY_ATTRIB(const restrict static n)], jmtx_matrix_brm *r,
                               float ck[JMTX_ARRAY_ATTRIB(const restrict m)],
                               float sk[JMTX_ARRAY_ATTRIB(const restrict m)],
                               float g[JMTX_ARRAY_ATTRIB(const restrict m)],
                               float alpha[JMTX_ARRAY_ATTRIB(const restrict m)],
                               float h[JMTX_ARRAY_ATTRIB(const restrict m)],
                               float p_mat[JMTX_ARRAY_ATTRIB(const restrict m *n)]);

uint32_t jmtx_gmresm_round_crs(const jmtx_matrix_crs *mtx, const uint32_t n, const uint32_t m, const float y_mag,
                               const float tol, const float residual[JMTX_ARRAY_ATTRIB(const restrict static n)],
                               float x[JMTX_ARRAY_ATTRIB(const restrict static n)], jmtx_matrix_brm *r,
                               float ck[JMTX_ARRAY_ATTRIB(const restrict m)],
                               float sk[JMTX_ARRAY_ATTRIB(const restrict m)],
                               float g[JMTX_ARRAY_ATTRIB(const restrict m)],
                               float alpha[JMTX_ARRAY_ATTRIB(const restrict m)],
                               float h[JMTX_ARRAY_ATTRIB(const restrict m)],
                               float p_mat[JMTX_ARRAY_ATTRIB(const restrict m *n)]);

#endif // JMTX_GMRES_INTERNAL_H
