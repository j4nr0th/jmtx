//
// Created by jan on 1.1.2024.
//

#ifndef JMTX_GMRES_INTERNAL_H
#define JMTX_GMRES_INTERNAL_H

#include "../matrices/sparse_diagonal_compressed_internal.h"
#include "../../../include/jmtx/float/matrices/band_row_major.h"

uint32_t jmtx_gmresm_round_cds(const jmtx_matrix_cds* mtx, const uint32_t n, const uint32_t m, const float y_mag,
                                const float tol, const float residual[const restrict static n],
                                float x[const restrict static n], jmtx_matrix_brm* r, float ck[const restrict m],
                                float sk[const restrict m], float g[const restrict m], float alpha[const restrict m],
                                float h[const restrict m], float p_mat[const restrict m * n]);

uint32_t jmtx_gmresm_round_crs(const jmtx_matrix_crs* mtx, const uint32_t n, const uint32_t m, const float y_mag,
                                const float tol, const float residual[const restrict static n],
                                float x[const restrict static n], jmtx_matrix_brm* r, float ck[const restrict m],
                                float sk[const restrict m], float g[const restrict m], float alpha[const restrict m],
                                float h[const restrict m], float p_mat[const restrict m * n]);

#endif //JMTX_GMRES_INTERNAL_H
