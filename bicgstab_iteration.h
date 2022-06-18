//
// Created by jan on 17.6.2022.
//

#ifndef MTXLIB_BICGSTAB_ITERATION_H
#define MTXLIB_BICGSTAB_ITERATION_H
#include "sparse_row_compressed.h"

/*
 * Biconjugate gradient stabilized method is an iterative method by H. A. van der Vorst (https://en.wikipedia.org/wiki/Biconjugate_gradient_stabilized_method)
 */

mtx_res_t bicgstab_crs(
        const CrsMatrix* mtx, const scalar_t* y, scalar_t* x, scalar_t convergence_dif, uint n_max_iter, uint* p_iter,
        scalar_t* p_error);

mtx_res_t bicgstab_crs_mt(
        const CrsMatrix* mtx, const scalar_t* y, scalar_t* x, scalar_t convergence_dif, uint n_max_iter, uint* p_iter,
        scalar_t* p_error, uint n_thrds);

#endif //MTXLIB_BICGSTAB_ITERATION_H
