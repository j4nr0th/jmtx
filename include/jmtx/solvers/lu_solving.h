//
// Created by jan on 6.11.2023.
//

#ifndef JMTX_LU_SOLVING_H
#define JMTX_LU_SOLVING_H

#include "../matrices/sparse_row_compressed.h"
#include "solver_base.h"

void jmtx_lu_solve(const jmtx_matrix_crs* l, const jmtx_matrix_crs* u, const float* restrict y, float* restrict x);

void jmtx_lu_solve_inplace(const jmtx_matrix_crs* l, const jmtx_matrix_crs* u, float* restrict x);

jmtx_result jmtx_incomplete_lu_decomposition_solve(
        const jmtx_matrix_crs* mtx, const float* y, float* x, float* aux_vec, jmtx_solver_arguments* args,
        const jmtx_allocator_callbacks* allocator_callbacks);

jmtx_result jmtx_incomplete_lu_decomposition_solve_precomputed(
        const jmtx_matrix_crs* mtx, const jmtx_matrix_crs* l, const jmtx_matrix_crs* u, const float* y, float* x,
        float* aux_vec, jmtx_solver_arguments* args);

jmtx_result jmtx_incomplete_lu_decomposition_solve_parallel(
        const jmtx_matrix_crs* mtx, const float* y, float* x, float* aux_vec, jmtx_solver_arguments* args,
        const jmtx_allocator_callbacks* allocator_callbacks);

jmtx_result jmtx_incomplete_lu_decomposition_solve_precomputed_parallel(
        const jmtx_matrix_crs* mtx, const jmtx_matrix_crs* l, const jmtx_matrix_crs* u, const float* y, float* x,
        float* aux_vec, jmtx_solver_arguments* args);

#endif //JMTX_LU_SOLVING_H