#include "gauss_seidel_iteration.h"
#include "../matrices/sparse_row_compressed.h"

#include <math.h>

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
 * @param args::opt_error_evolution (optional) pointer to an array of length max_iterations, that receives the error
 * value of each iteration
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_NOT_CONVERGED if it hasn't reached given stopping criterion,
 * in case of failure it returns the associated error code
 */
jmtx_result JMTX_NAME_TYPED(solve_iterative_gauss_seidel_crs)(const JMTX_NAME_TYPED(matrix_crs) * mtx,
                                                              const JMTX_SCALAR_T *restrict y,
                                                              JMTX_SCALAR_T *restrict x,
                                                              JMTX_SCALAR_T *restrict aux_vec1,
                                                              JMTX_NAME_TYPED(solver_arguments) * args)
{
    //  Length of x and y
    const JMTX_INDEX_T n = mtx->base.cols;
    JMTX_SCALAR_T *const div_factors = aux_vec1;
    JMTX_REAL_T mag_y = 0;
    for (JMTX_INDEX_T i = 0; i < n; ++i)
    {
        const JMTX_SCALAR_T d = JMTX_NAME_TYPED(matrix_crs_get_entry)(mtx, i, i);
        div_factors[i] = 1 / d;
        mag_y += JMTX_DOT(y[i], y[i]);
    }
    mag_y = JMTX_REAL_ROOT(mag_y);

    JMTX_REAL_T err;
    JMTX_INDEX_T n_iterations = 0;
    for (;;)
    {

        for (JMTX_INDEX_T i = 0; i < n; ++i)
        {
            const JMTX_SCALAR_T res = JMTX_NAME_TYPED(matrix_crs_vector_multiply_row)(mtx, x, i);
            const JMTX_SCALAR_T new_x = (y[i] - res) * div_factors[i] + x[i];
            x[i] = new_x;
        }

        err = 0;
        for (JMTX_INDEX_T i = 0; i < n; ++i)
        {
            const JMTX_SCALAR_T val = y[i] - JMTX_NAME_TYPED(matrix_crs_vector_multiply_row)(mtx, x, i);
            err += JMTX_DOT(val, val);
        }
        err = JMTX_REAL_ROOT(err) / mag_y;
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
