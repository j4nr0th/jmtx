//
// Created by jan on 1.1.2024.
//

#ifndef JMTXD_GMRES_INTERNAL_H
#define JMTXD_GMRES_INTERNAL_H

#include "../../../include/jmtx/double/matrices/band_row_major.h"
#include "../matrices/sparse_diagonal_compressed_internal.h"

uint32_t jmtxd_gmresm_round_cds(const jmtxd_matrix_cds *mtx, const uint32_t n, const uint32_t m, const double y_mag,
                                const double tol, const double residual[JMTX_ARRAY_ATTRIB(const restrict static n)],
                                double x[JMTX_ARRAY_ATTRIB(const restrict static n)], jmtxd_matrix_brm *r,
                                double ck[JMTX_ARRAY_ATTRIB(const restrict m)],
                                double sk[JMTX_ARRAY_ATTRIB(const restrict m)],
                                double g[JMTX_ARRAY_ATTRIB(const restrict m)],
                                double alpha[JMTX_ARRAY_ATTRIB(const restrict m)],
                                double h[JMTX_ARRAY_ATTRIB(const restrict m)],
                                double p_mat[JMTX_ARRAY_ATTRIB(const restrict m *n)]);

uint32_t jmtxd_gmresm_round_crs(const jmtxd_matrix_crs *mtx, const uint32_t n, const uint32_t m, const double y_mag,
                                const double tol, const double residual[JMTX_ARRAY_ATTRIB(const restrict static n)],
                                double x[JMTX_ARRAY_ATTRIB(const restrict static n)], jmtxd_matrix_brm *r,
                                double ck[JMTX_ARRAY_ATTRIB(const restrict m)],
                                double sk[JMTX_ARRAY_ATTRIB(const restrict m)],
                                double g[JMTX_ARRAY_ATTRIB(const restrict m)],
                                double alpha[JMTX_ARRAY_ATTRIB(const restrict m)],
                                double h[JMTX_ARRAY_ATTRIB(const restrict m)],
                                double p_mat[JMTX_ARRAY_ATTRIB(const restrict m *n)]);

#endif // JMTXD_GMRES_INTERNAL_H
