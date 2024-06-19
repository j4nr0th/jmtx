//
// Created by jan on 1.1.2024.
//

#ifndef JMTXC_GMRES_INTERNAL_H
#define JMTXC_GMRES_INTERNAL_H

#include "../matrices/sparse_row_compressed_internal.h"
#include "../matrices/sparse_diagonal_compressed_internal.h"
#include "../../../include/jmtx/cfloat/matrices/band_row_major.h"

uint32_t jmtxc_gmresm_round_cds(const jmtxc_matrix_cds* mtx, const uint32_t n, const uint32_t m, const _Complex float y_mag,
                                const float tol, const _Complex float residual[const restrict static n],
                                _Complex float x[const restrict static n], jmtxc_matrix_brm* r, _Complex float ck[const restrict m],
                                _Complex float sk[const restrict m], _Complex float g[const restrict m], _Complex float alpha[const restrict m],
                                _Complex float h[const restrict m], _Complex float p_mat[const restrict m * n]);

uint32_t jmtxc_gmresm_round_crs(const jmtxc_matrix_crs* mtx, const uint32_t n, const uint32_t m, const _Complex float y_mag,
                                const float tol, const _Complex float residual[const restrict static n],
                                _Complex float x[const restrict static n], jmtxc_matrix_brm* r, _Complex float ck[const restrict m],
                                _Complex float sk[const restrict m], _Complex float g[const restrict m], _Complex float alpha[const restrict m],
                                _Complex float h[const restrict m], _Complex float p_mat[const restrict m * n]);

#endif //JMTXC_GMRES_INTERNAL_H
