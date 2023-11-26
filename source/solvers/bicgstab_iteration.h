//
// Created by jan on 17.6.2022.
//

#ifndef JMTX_BICGSTAB_ITERATION_H
#define JMTX_BICGSTAB_ITERATION_H
#ifndef JMTX_SPARSE_ROW_COMPRESSED_H
    #include "../../include/jmtx/matrices/sparse_row_compressed.h"
#endif


/*
 * Biconjugate gradient stabilized method is an iterative method by H. A. van der Vorst (https://en.wikipedia.org/wiki/Biconjugate_gradient_stabilized_method)
 */

jmtx_result jmtx_bicgstab_crs(
        const jmtx_matrix_crs* mtx, const float* y, float* x, float convergence_dif,
        uint32_t n_max_iter, uint32_t* p_iter, float* p_final_error, float* p_error,
        const jmtx_allocator_callbacks* allocator_callbacks);

jmtx_result jmtx_bicgstab_crs_mt(
        const jmtx_matrix_crs* mtx, const float* y, float* x, float convergence_dif,
        uint32_t n_max_iter, uint32_t* p_iter, float* p_final_error, float* p_error,
        const jmtx_allocator_callbacks* allocator_callbacks, uint32_t n_thrds);

#endif //JMTX_BICGSTAB_ITERATION_H
