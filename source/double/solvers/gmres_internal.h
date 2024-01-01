//
// Created by jan on 1.1.2024.
//

#ifndef JMTXD_GMRES_INTERNAL_H
#define JMTXD_GMRES_INTERNAL_H

#include "../matrices/sparse_diagonal_compressed_internal.h"
#include "../../../include/jmtx/double/matrices/band_row_major.h"

uint32_t jmtxd_gmresm_round_cds(const jmtxd_matrix_cds* mtx, const uint32_t n, const uint32_t m, const double y_mag,
                                const double tol, const double residual[const restrict static n],
                                double x[const restrict static n], jmtxd_matrix_brm* r, double ck[const restrict m],
                                double sk[const restrict m], double g[const restrict m], double alpha[const restrict m],
                                double h[const restrict m], double p_mat[const restrict m * n]);

#endif //JMTXD_GMRES_INTERNAL_H
