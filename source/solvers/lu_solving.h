//
// Created by jan on 6.11.2023.
//

#ifndef JMTX_LU_SOLVING_H
#define JMTX_LU_SOLVING_H
#include "../matrices/sparse_row_compressed.h"

void jmtx_lu_solve(const jmtx_matrix_crs* l, const jmtx_matrix_crs* u, const float* restrict y, float* restrict x);

void jmtx_lu_solve_inplace(const jmtx_matrix_crs* l, const jmtx_matrix_crs* u, float* restrict x);

#endif //JMTX_LU_SOLVING_H
