// Automatically generated from source/float/solvers/gauss_seidel_iteration.c on Sun Dec 17 16:44:56 2023
//
// Created by jan on 16.6.2022.
//

#include "../matrices/sparse_row_compressed_internal.h"
#include "../../../include/jmtx/cfloat/solvers/gauss_seidel_iteration.h"

#include <math.h>
#include <complex.h>



/*
 * Gauss-Seidel is an iterative method for solving the system Ax = y. It works by splitting the matrix A into
 * matrices D (diagonal), L (lower triangular with zero diagonal), and U (upper triangular with zero diagonal), such
 * that A = D + L + U. The equation is then expressed as x_(n+1) = (D + L)^{-1} (y - U x_(n)). This means that all the
 * functions in this file require that the diagonals of the matrices are non-zero.
 */

/**
 * Uses Gauss-Seidel (https://en.wikipedia.org/wiki/Gauss%E2%80%93Seidel_method)
 * to solve the linear system Ax = y
 *
 * @param mtx pointer to the memory where matrix A is stored as a band row major matrix
 * @param y pointer to the memory where the vector y is stored
 * @param x pointer to the memory where the solution vector x will be stored
 * @param aux_vec auxiliary memory for a vector of the same size as x and y
 * @param args::in_convergence_criterion tolerance to determine if the solution is close enough
 * @param args::in_max_iterations number of iterations to stop at
 * @param args::out_last_error receives the value of the error criterion at the final iteration
 * @param args::out_last_iteration receives the number of the final iteration
 * @param args::opt_error_evolution (optional) pointer to an array of length max_iterations, that receives the error value of each
 * iteration
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_NOT_CONVERGED if it hasn't reached given stopping criterion,
 * in case of failure it returns the associated error code
 */
jmtx_result jmtxc_solve_iterative_gauss_seidel_crs(const jmtxc_matrix_crs* mtx, const _Complex float* restrict y, _Complex float* restrict x,
                                  _Complex float* restrict aux_vec1, jmtx_solver_arguments* args)
{
    //  Length of x and y
    const uint32_t n = mtx->base.cols;
    _Complex float* const div_factors = aux_vec1;
    float mag_y = 0;
    for (uint32_t i = 0; i < n; ++i)
    {
        const _Complex float d = jmtxc_matrix_crs_get_entry(mtx, i, i);
        div_factors[i] = 1 / d;
        mag_y += conjf(y[i]) * y[i];
    }
    mag_y = sqrtf(mag_y);


    float err;
    uint32_t n_iterations = 0;
    for(;;)
    {

        for (uint32_t i = 0; i < n; ++i)
        {
            _Complex float res = jmtxc_matrix_crs_vector_multiply_row(mtx, x, i);
            const _Complex float new_x = (y[i] - res) * div_factors[i] + x[i];
            x[i] = new_x;
        }

        err = 0;
        for (uint32_t i = 0; i < n; ++i)
        {
            const _Complex float val = y[i] - jmtxc_matrix_crs_vector_multiply_row(mtx, x, i);
            err += conjf(val) * val;
        }
        err = sqrtf(err) / mag_y;
        if (args->opt_error_evolution)
        {
            args->opt_error_evolution[n_iterations] = err;
        }

        if (n_iterations == args->in_max_iterations)
        {
            break;
        }
        if (args->opt_error_evolution)
        {
            args->opt_error_evolution[n_iterations] = err;
        }
        n_iterations += 1;
        if (err < args->in_convergence_criterion)
        {
            break;
        }
    }

    args->out_last_iteration = n_iterations;
    args->out_last_error = err;

    if (!isfinite(err) || err >= args->in_convergence_criterion)
    {
        return JMTX_RESULT_NOT_CONVERGED;
    }

    return JMTX_RESULT_SUCCESS;
}


static inline int check_vector_overlaps(const unsigned n, const size_t size, const void* ptrs[static const n])
{
    for (unsigned i = 0; i < n; ++i)
    {
        const uintptr_t p1 = (uintptr_t)ptrs[i];
        for (unsigned j = i + 1; j < n; ++j)
        {
            const uintptr_t p2 = (uintptr_t)ptrs[j];
            if (p1 > p2)
            {
                if (p2 + size > p1)
                {
                    return 1;
                }
            }
            else if (p1 + size > p2)
            {
                return 1;
            }
        }
    }

    return 0;
}

/**
 * Uses Gauss-Seidel (https://en.wikipedia.org/wiki/Gauss%E2%80%93Seidel_method)
 * to solve the linear system Ax = y
 *
 * @param mtx pointer to the memory where matrix A is stored as a band row major matrix
 * @param y pointer to the memory where the vector y is stored
 * @param x pointer to the memory where the solution vector x will be stored
 * @param aux_vec auxiliary memory for a vector of the same size as x and y
 * @param args::in_convergence_criterion tolerance to determine if the solution is close enough
 * @param args::in_max_iterations number of iterations to stop at
 * @param args::out_last_error receives the value of the error criterion at the final iteration
 * @param args::out_last_iteration receives the number of the final iteration
 * @param args::opt_error_evolution (optional) pointer to an array of length max_iterations, that receives the error value of each
 * iteration
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_NOT_CONVERGED if it hasn't reached given stopping criterion,
 * in case of failure it returns the associated error code
 */
jmtx_result jmtxcs_solve_iterative_gauss_seidel_crs(const jmtxc_matrix_crs* mtx, uint32_t n, const _Complex float y[static restrict n],
                                   _Complex float x[restrict n], _Complex float aux_vec1[restrict n], jmtx_solver_arguments* args)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTXC_TYPE_CRS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (mtx->base.rows != n || mtx->base.cols != n)
    {
        return JMTX_RESULT_BAD_MATRIX;
    }
    if (!x)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!aux_vec1)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!args)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    const void* ptrs[] = {y, x, aux_vec1};
    if (check_vector_overlaps(sizeof(ptrs)/sizeof(*ptrs), n * sizeof(*x), ptrs))
    {
        return JMTX_RESULT_BAD_PARAM;
    }
    return jmtxc_solve_iterative_gauss_seidel_crs(mtx, y, x, aux_vec1, args);
}
