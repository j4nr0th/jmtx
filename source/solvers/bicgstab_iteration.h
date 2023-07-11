//
// Created by jan on 17.6.2022.
//

#ifndef MTXLIB_BICGSTAB_ITERATION_H
#define MTXLIB_BICGSTAB_ITERATION_H
#include "../matrices/sparse_row_compressed.h"
#ifdef MTX_ERROR_MESSAGES
#include "errors.h"
#endif


/*
 * Biconjugate gradient stabilized method is an iterative method by H. A. van der Vorst (https://en.wikipedia.org/wiki/Biconjugate_gradient_stabilized_method)
 */

jmtx_result bicgstab_crs(
        const jmtx_matrix_crs* mtx, const jmtx_scalar_t* y, jmtx_scalar_t* x, jmtx_scalar_t convergence_dif,
        uint32_t n_max_iter, uint32_t* p_iter, jmtx_scalar_t* p_final_error, jmtx_scalar_t* p_error,
        const jmtx_allocator_callbacks* allocator_callbacks);

jmtx_result bicgstab_crs_mt(
        const jmtx_matrix_crs* mtx, const jmtx_scalar_t* y, jmtx_scalar_t* x, jmtx_scalar_t convergence_dif,
        uint32_t n_max_iter, uint32_t* p_iter, jmtx_scalar_t* p_final_error, jmtx_scalar_t* p_error,
        const jmtx_allocator_callbacks* allocator_callbacks, uint32_t n_thrds);

#endif //MTXLIB_BICGSTAB_ITERATION_H
