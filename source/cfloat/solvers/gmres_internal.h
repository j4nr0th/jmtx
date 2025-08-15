//
// Created by jan on 1.1.2024.
//

#ifndef JMTXC_GMRES_INTERNAL_H
#define JMTXC_GMRES_INTERNAL_H

#include "../../../include/jmtx/cfloat/matrices/band_row_major.h"
#include "../matrices/sparse_diagonal_compressed_internal.h"
#include "../matrices/sparse_row_compressed_internal.h"

uint32_t jmtxc_gmresm_round_cds(const jmtxc_matrix_cds *mtx, const uint32_t n, const uint32_t m,
                                const _Complex float y_mag, const float tol,
                                const _Complex float residual[JMTX_ARRAY_ATTRIB(const restrict static n)],
                                _Complex float x[JMTX_ARRAY_ATTRIB(const restrict static n)], jmtxc_matrix_brm *r,
                                _Complex float ck[JMTX_ARRAY_ATTRIB(const restrict m)],
                                _Complex float sk[JMTX_ARRAY_ATTRIB(const restrict m)],
                                _Complex float g[JMTX_ARRAY_ATTRIB(const restrict m)],
                                _Complex float alpha[JMTX_ARRAY_ATTRIB(const restrict m)],
                                _Complex float h[JMTX_ARRAY_ATTRIB(const restrict m)],
                                _Complex float p_mat[JMTX_ARRAY_ATTRIB(const restrict m *n)]);

uint32_t jmtxc_gmresm_round_crs(const jmtxc_matrix_crs *mtx, const uint32_t n, const uint32_t m,
                                const _Complex float y_mag, const float tol,
                                const _Complex float residual[JMTX_ARRAY_ATTRIB(const restrict static n)],
                                _Complex float x[JMTX_ARRAY_ATTRIB(const restrict static n)], jmtxc_matrix_brm *r,
                                _Complex float ck[JMTX_ARRAY_ATTRIB(const restrict m)],
                                _Complex float sk[JMTX_ARRAY_ATTRIB(const restrict m)],
                                _Complex float g[JMTX_ARRAY_ATTRIB(const restrict m)],
                                _Complex float alpha[JMTX_ARRAY_ATTRIB(const restrict m)],
                                _Complex float h[JMTX_ARRAY_ATTRIB(const restrict m)],
                                _Complex float p_mat[JMTX_ARRAY_ATTRIB(const restrict m *n)]);

#endif // JMTXC_GMRES_INTERNAL_H
