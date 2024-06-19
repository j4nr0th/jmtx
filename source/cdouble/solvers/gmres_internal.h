//
// Created by jan on 1.1.2024.
//

#ifndef JMTXZ_GMRES_INTERNAL_H
#define JMTXZ_GMRES_INTERNAL_H

#include "../matrices/sparse_diagonal_compressed_internal.h"
#include "../matrices/sparse_row_compressed_internal.h"
#include "../../../include/jmtx/cdouble/matrices/band_row_major.h"

uint32_t jmtxz_gmresm_round_cds(const jmtxz_matrix_cds* mtx, const uint32_t n, const uint32_t m, const _Complex double y_mag,
                                const double tol, const _Complex double residual[const restrict static n],
                                _Complex double x[const restrict static n], jmtxz_matrix_brm* r, _Complex double ck[const restrict m],
                                _Complex double sk[const restrict m], _Complex double g[const restrict m], _Complex double alpha[const restrict m],
                                _Complex double h[const restrict m], _Complex double p_mat[const restrict m * n]);

uint32_t jmtxz_gmresm_round_crs(const jmtxz_matrix_crs* mtx, const uint32_t n, const uint32_t m, const _Complex double y_mag,
                                const double tol, const _Complex double residual[const restrict static n],
                                _Complex double x[const restrict static n], jmtxz_matrix_brm* r, _Complex double ck[const restrict m],
                                _Complex double sk[const restrict m], _Complex double g[const restrict m], _Complex double alpha[const restrict m],
                                _Complex double h[const restrict m], _Complex double p_mat[const restrict m * n]);

#endif //JMTXZ_GMRES_INTERNAL_H
