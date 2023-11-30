// Automatically generated from source/float/solvers/bicgstab_iteration.h on Thu Nov 30 20:33:08 2023
//
// Created by jan on 17.6.2022.
//

#ifndef JMTX_BICGSTAB_ITERATION_H
#define JMTX_BICGSTAB_ITERATION_H
#ifndef JMTX_SPARSE_ROW_COMPRESSED_H
    #include "../../../include/jmtx/double/matrices/sparse_row_compressed.h"
#endif


/*
 * Biconjugate gradient stabilized method is an iterative method by H. A. van der Vorst (https://en.wikipedia.org/wiki/Biconjugate_gradient_stabilized_method)
 */

jmtx_result jmtxd_bicgstab_crs(
        const jmtxd_matrix_crs* mtx, const double* y, double* x, double convergence_dif,
        uint32_t n_max_iter, uint32_t* p_iter, double* p_final_error, double* p_error,
        const jmtx_allocator_callbacks* allocator_callbacks);

jmtx_result jmtxd_bicgstab_crs_mt(
        const jmtxd_matrix_crs* mtx, const double* y, double* x, double convergence_dif,
        uint32_t n_max_iter, uint32_t* p_iter, double* p_final_error, double* p_error,
        const jmtx_allocator_callbacks* allocator_callbacks, uint32_t n_thrds);

#endif //JMTX_BICGSTAB_ITERATION_H
