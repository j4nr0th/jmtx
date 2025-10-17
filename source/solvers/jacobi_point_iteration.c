#include "jacobi_point_iteration.h"
#include "../matrices/band_row_major.h"
#include "../matrices/sparse_diagonal_compressed.h"
#include "../matrices/sparse_row_compressed.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>

#include <omp.h>

/**
 * Uses Jacobi point iteration (also known as Jacobi method: https://en.wikipedia.org/wiki/Jacobi_method)
 * to solve the linear system Ax = y
 *
 * @param mtx pointer to the memory where matrix A is stored as a compressed row sparse matrix
 * @param y pointer to the memory where the vector y is stored
 * @param x pointer to the memory where the solution vector x will be written
 * @param aux_vec1 auxiliary memory for a vector of the same size as x and y
 * @param aux_vec2 auxiliary memory for a vector of the same size as x and y
 * @param args::in_convergence_criterion tolerance to determine if the solution is close enough
 * @param args::in_max_iterations number of iterations to stop at
 * @param args::out_last_error receives the value of the error criterion at the final iteration
 * @param args::out_last_iteration receives the number of the final iteration
 * @param args::opt_error_evolution (optional) pointer to an array of length max_iterations, that receives the error
 * value of each iteration
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_NOT_CONVERGED if it hasn't reached given stopping criterion,
 * in case of failure it returns the associated error code
 */
jmtx_result JMTX_NAME_TYPED(solve_iterative_jacobi_crs)(const JMTX_NAME_TYPED(matrix_crs) * mtx,
                                                        const JMTX_SCALAR_T *restrict y, JMTX_SCALAR_T *restrict x,
                                                        JMTX_SCALAR_T *restrict aux_vec1,
                                                        JMTX_SCALAR_T *restrict aux_vec2,
                                                        JMTX_NAME_TYPED(solver_arguments) * args)
{
    //  Length of x and y
    const JMTX_INDEX_T n = mtx->base.cols;
    JMTX_SCALAR_T *const div_factor = aux_vec1;
    JMTX_REAL_T y_mag = 0;
    //  Initial guess by assuming that mtx is a diagonal matrix
    for (JMTX_INDEX_T i = 0; i < n; ++i)
    {
        JMTX_SCALAR_T d = JMTX_NAME_TYPED(matrix_crs_get_entry)(mtx, i, i);
        x[i] = y[i] / d;
        div_factor[i] = 1.0f / d;
        y_mag += JMTX_DOT(y[i], y[i]);
    }
    y_mag = JMTX_REAL_ROOT(y_mag);

    //  Memory used to store result of the current iteration
    JMTX_SCALAR_T *const auxiliary_x = aux_vec2;

    JMTX_SCALAR_T *x0 = auxiliary_x;
    JMTX_SCALAR_T *x1 = x;
    JMTX_REAL_T err;
    JMTX_INDEX_T n_iterations = 0;
    do
    {
        err = 0.0f;
        {
            JMTX_SCALAR_T *tmp = x1;
            x1 = x0;
            x0 = tmp;
        }

        //  For each entry, find the corresponding row in matrix A - D and compute the dot product between x and that
        //  row
        for (JMTX_INDEX_T i = 0; i < n; ++i)
        {
            JMTX_SCALAR_T *row_ptr;
            JMTX_INDEX_T *index_ptr;
            JMTX_INDEX_T n_elements = JMTX_NAME_TYPED(matrix_crs_get_row)(mtx, i, &index_ptr, &row_ptr);
            JMTX_SCALAR_T res = 0;
            for (JMTX_INDEX_T j = 0; j < n_elements; ++j)
            {
                if (i != index_ptr[j])
                {
                    res += row_ptr[j] * x0[index_ptr[j]];
                }
            }
            //  Multiplication of vector x by D⁻¹
            x1[i] = (y[i] - res) * div_factor[i];
        }

        for (JMTX_INDEX_T i = 0; i < n; ++i)
        {
            JMTX_SCALAR_T val = JMTX_NAME_TYPED(matrix_crs_vector_multiply_row)(mtx, x1, i);
            val -= y[i];
            err += JMTX_DOT(val, val);
        }
        //  Have "err" as ratio between magnitude of y vector and magnitude of residual
        err = JMTX_REAL_ROOT(err) / (JMTX_REAL_T)y_mag;
        if (args->opt_error_evolution)
        {
            args->opt_error_evolution[n_iterations] = err;
        }
        n_iterations += 1;
    } while (err > args->in_convergence_criterion && n_iterations < args->in_max_iterations);

    if (x1 == auxiliary_x)
    {
        memcpy(x, auxiliary_x, sizeof *x * n);
    }

    args->out_last_error = err;
    args->out_last_iteration = n_iterations;
    //    LEAVE_FUNCTION();
    return n_iterations == args->in_max_iterations ? JMTX_RESULT_NOT_CONVERGED : JMTX_RESULT_SUCCESS;
}

/**
 * Uses Jacobi point iteration (also known as Jacobi method: https://en.wikipedia.org/wiki/Jacobi_method)
 * to solve the linear system Ax = y. Uses a relaxation factor ω for the following relation:
 *
 *  x(n + 1) = ω D⁻¹ (y - A x(n)) + x(n)
 *
 *  This may allow for better convergence
 *
 * @param mtx pointer to the memory where matrix A is stored as a compressed row sparse matrix
 * @param y pointer to the memory where the vector y is stored
 * @param x pointer to the memory where the solution vector x will be stored
 * @param relaxation_factor relaxation factor used for the iterations, which must be greater than zero
 * @param aux_vec1 auxiliary memory for a vector of the same size as x and y
 * @param aux_vec2 auxiliary memory for a vector of the same size as x and y
 * @param args::in_convergence_criterion tolerance to determine if the solution is close enough
 * @param args::in_max_iterations number of iterations to stop at
 * @param args::out_last_error receives the value of the error criterion at the final iteration
 * @param args::out_last_iteration receives the number of the final iteration
 * @param args::opt_error_evolution (optional) pointer to an array of length max_iterations, that receives the error
 * value of each iteration
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_NOT_CONVERGED if it hasn't reached given stopping criterion,
 * in case of failure it returns the associated error code
 */
jmtx_result JMTX_NAME_TYPED(solve_iterative_jacobi_relaxed_crs)(
    const JMTX_NAME_TYPED(matrix_crs) * mtx, const JMTX_SCALAR_T *restrict y, JMTX_SCALAR_T *restrict x,
    JMTX_SCALAR_T relaxation_factor, JMTX_SCALAR_T *restrict aux_vec1, JMTX_SCALAR_T *restrict aux_vec2,
    JMTX_NAME_TYPED(solver_arguments) * args)
{
    //  Length of x and y
    const JMTX_INDEX_T n = mtx->base.cols;
    JMTX_SCALAR_T *div_factor = aux_vec1;
    JMTX_REAL_T y_mag = 0;
    //  Initial guess by assuming that mtx is a diagonal matrix
    for (JMTX_INDEX_T i = 0; i < n; ++i)
    {
        const JMTX_SCALAR_T d = JMTX_NAME_TYPED(matrix_crs_get_entry)(mtx, i, i);
        x[i] = y[i] / d;
        div_factor[i] = relaxation_factor / d;
        y_mag += JMTX_DOT(y[i], y[i]);
    }
    y_mag = JMTX_REAL_ROOT(y_mag);

    //  Memory used to store result of the current iteration
    JMTX_SCALAR_T *const auxiliary_x = aux_vec2;

    JMTX_SCALAR_T *x0 = auxiliary_x;
    JMTX_SCALAR_T *x1 = x;
    JMTX_REAL_T err;
    JMTX_INDEX_T n_iterations = 0;
    do
    {
        err = 0.0f;
        {
            JMTX_SCALAR_T *tmp = x1;
            x1 = x0;
            x0 = tmp;
        }

        //  For each entry, find the corresponding row in matrix A - D and compute the dot product between x and that
        //  row
        for (JMTX_INDEX_T i = 0; i < n; ++i)
        {
            const JMTX_SCALAR_T res = JMTX_NAME_TYPED(matrix_crs_vector_multiply_row)(mtx, x0, i);
            x1[i] = x0[i] + (y[i] - res) * div_factor[i];
        }

        JMTX_NAME_TYPED(matrix_crs_vector_multiply)(mtx, x1, x0);
        for (JMTX_INDEX_T i = 0; i < n; ++i)
        {
            const JMTX_SCALAR_T val = y[i] - x0[i];
            err += JMTX_DOT(val, val);
        }
        err = JMTX_REAL_ROOT(err) / (JMTX_REAL_T)y_mag;
        if (args->opt_error_evolution)
        {
            args->opt_error_evolution[n_iterations] = err;
        }
        n_iterations += 1;
    } while (err > args->in_convergence_criterion && n_iterations < args->in_max_iterations);

    if (x1 == auxiliary_x)
    {
        memcpy(x, auxiliary_x, sizeof *x * n);
    }

    args->out_last_iteration = n_iterations;
    args->out_last_error = err;

    return n_iterations == args->in_max_iterations ? JMTX_RESULT_NOT_CONVERGED : JMTX_RESULT_SUCCESS;
}

/**
 * Uses Jacobi point iteration (also known as Jacobi method: https://en.wikipedia.org/wiki/Jacobi_method)
 * to solve the linear system Ax = y.
 *
 * This version of the function uses OpenMP to solve the system in parallel.
 *
 * @param mtx pointer to the memory where matrix A is stored as a compressed row sparse matrix
 * @param y pointer to the memory where the vector y is stored
 * @param x pointer to the memory where the solution vector x will be stored
 * @param aux_vec1 auxiliary memory for a vector of the same size as x and y
 * @param aux_vec2 auxiliary memory for a vector of the same size as x and y
 * @param args::in_convergence_criterion tolerance to determine if the solution is close enough
 * @param args::in_max_iterations number of iterations to stop at
 * @param args::out_last_error receives the value of the error criterion at the final iteration
 * @param args::out_last_iteration receives the number of the final iteration
 * @param args::opt_error_evolution (optional) pointer to an array of length max_iterations, that receives the error
 * value of each iteration
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_NOT_CONVERGED if it hasn't reached given stopping criterion,
 * in case of failure it returns the associated error code
 */
jmtx_result JMTX_NAME_TYPED(solve_iterative_jacobi_crs_parallel)(
    const JMTX_NAME_TYPED(matrix_crs) * mtx, const JMTX_SCALAR_T *restrict y, JMTX_SCALAR_T *restrict x,
    JMTX_SCALAR_T *restrict aux_vector1, JMTX_SCALAR_T *restrict aux_vector2, JMTX_NAME_TYPED(solver_arguments) * args)
{
    //  Length of x and y
    const JMTX_INDEX_T n = mtx->base.cols;
    JMTX_REAL_T y_mag = 0;
    //  Memory used to store result of the current iteration
    JMTX_SCALAR_T *x0 = aux_vector2;
    JMTX_SCALAR_T *x1 = x;
    JMTX_REAL_T err = 0.0f;
    JMTX_INDEX_T n_iterations = 0;
    JMTX_REAL_T *const p_error_evolution = args->opt_error_evolution;
    const JMTX_REAL_T convergence_dif = args->in_convergence_criterion;
    const JMTX_INDEX_T max_iterations = args->in_max_iterations;
#pragma omp parallel default(none) shared(err, x0, x1, y, mtx, aux_vector1, y_mag, p_error_evolution, n_iterations,    \
                                              convergence_dif, max_iterations, n)
    {
#pragma omp for schedule(static)
        for (JMTX_INDEX_T i = 0; i < n; ++i)
        {
            const JMTX_SCALAR_T d = JMTX_NAME_TYPED(matrix_crs_get_entry)(mtx, i, i);
            aux_vector1[i] = 1.0f / d;
            const JMTX_SCALAR_T mag = JMTX_DOT(y[i], y[i]);
            y_mag += mag;
        }

#pragma omp master
        {
            y_mag = JMTX_REAL_ROOT(y_mag);
        }

        do
        {
#pragma omp barrier
#pragma omp master
            {
                err = 0.0f;
                JMTX_SCALAR_T *tmp = x1;
                x1 = x0;
                x0 = tmp;
            }
#pragma omp barrier

            //  For each entry, find the corresponding row in matrix A - D and compute the dot product between x and
            //  that row
#pragma omp for schedule(static)
            for (JMTX_INDEX_T i = 0; i < n; ++i)
            {
                x1[i] = x0[i] + (y[i] - JMTX_NAME_TYPED(matrix_crs_vector_multiply_row)(mtx, x0, i)) * aux_vector1[i];
            }

#pragma omp for reduction(+ : err) schedule(static)
            for (JMTX_INDEX_T i = 0; i < n; ++i)
            {
                const JMTX_SCALAR_T val = JMTX_NAME_TYPED(matrix_crs_vector_multiply_row)(mtx, x1, i) - y[i];
                err += JMTX_DOT(val, val);
            }

#pragma omp master
            {
                err = JMTX_REAL_ROOT(err) / (JMTX_REAL_T)y_mag;
                if (p_error_evolution != NULL)
                {
                    p_error_evolution[n_iterations] = err;
                }
                n_iterations += 1;
            }
#pragma omp barrier
        } while (err > convergence_dif && n_iterations < max_iterations);
    }

    if (x1 == aux_vector2)
    {
        memcpy(x, aux_vector2, sizeof *x * n);
    }
    args->out_last_iteration = n_iterations;
    args->out_last_error = err;

    return n_iterations == max_iterations || !isfinite(err) ? JMTX_RESULT_NOT_CONVERGED : JMTX_RESULT_SUCCESS;
}

/**
 * Uses Jacobi point iteration (also known as Jacobi method: https://en.wikipedia.org/wiki/Jacobi_method)
 * to solve the linear system Ax = y
 *
 * @param mtx pointer to the memory where matrix A is stored as a compressed diagonal sparse matrix
 * @param y pointer to the memory where the vector y is stored
 * @param x pointer to the memory where the solution vector x will be stored
 * @param aux_vec1 auxiliary memory for a vector of the same size as x and y
 * @param aux_vec2 auxiliary memory for a vector of the same size as x and y
 * @param args::in_convergence_criterion tolerance to determine if the solution is close enough
 * @param args::in_max_iterations number of iterations to stop at
 * @param args::out_last_error receives the value of the error criterion at the final iteration
 * @param args::out_last_iteration receives the number of the final iteration
 * @param args::opt_error_evolution (optional) pointer to an array of length max_iterations, that receives the error
 * value of each iteration
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_NOT_CONVERGED if it hasn't reached given stopping criterion,
 * in case of failure it returns the associated error code
 */
jmtx_result JMTX_NAME_TYPED(solve_iterative_jacobi_cds)(const JMTX_NAME_TYPED(matrix_cds) * mtx,
                                                        const JMTX_SCALAR_T *restrict y, JMTX_SCALAR_T *restrict x,
                                                        JMTX_SCALAR_T *restrict aux_vec1,
                                                        JMTX_SCALAR_T *restrict aux_vec2,
                                                        JMTX_NAME_TYPED(solver_arguments) * args)
{
    //  Length of x and y
    const JMTX_INDEX_T n = mtx->base.cols;
    JMTX_SCALAR_T *const div_factor = aux_vec1;
    JMTX_REAL_T y_mag = 0;
    //  Initial guess by assuming that mtx is a diagonal matrix
    for (JMTX_INDEX_T i = 0; i < n; ++i)
    {
        JMTX_SCALAR_T d = mtx->main_diagonal[i];
        x[i] = y[i] / d;
        div_factor[i] = 1.0f / d;
        y_mag += JMTX_DOT(y[i], y[i]);
    }
    y_mag = JMTX_REAL_ROOT(y_mag);

    //  Memory used to store result of the current iteration
    JMTX_SCALAR_T *const auxiliary_x = aux_vec2;

    JMTX_SCALAR_T *x0 = auxiliary_x;
    JMTX_SCALAR_T *x1 = x;
    JMTX_REAL_T err;
    JMTX_INDEX_T n_iterations = 0;
    do
    {
        err = 0.0f;
        {
            JMTX_SCALAR_T *tmp = x1;
            x1 = x0;
            x0 = tmp;
        }

        //  For each entry, find the corresponding row in matrix A - D and compute the dot product between x and that
        //  row
        JMTX_NAME_TYPED(matrix_cds_vector_multiply)(mtx, x0, x1);
        for (JMTX_INDEX_T i = 0; i < n; ++i)
        {
            //  Multiplication of vector x by D⁻¹
            x1[i] = (y[i] - x1[i]) * div_factor[i];
        }

        JMTX_NAME_TYPED(matrix_cds_vector_multiply)(mtx, x1, x0);
        for (JMTX_INDEX_T i = 0; i < n; ++i)
        {
            const JMTX_SCALAR_T val = y[i] - x0[i];
            err += JMTX_DOT(val, val);
        }
        //  Have "err" as ratio between magnitude of y vector and magnitude of residual
        err = JMTX_REAL_ROOT(err) / (JMTX_REAL_T)y_mag;
        if (args->opt_error_evolution)
        {
            args->opt_error_evolution[n_iterations] = err;
        }
        n_iterations += 1;
    } while (err > args->in_convergence_criterion && n_iterations < args->in_max_iterations);

    if (x1 == auxiliary_x)
    {
        memcpy(x, auxiliary_x, sizeof *x * n);
    }

    args->out_last_error = err;
    args->out_last_iteration = n_iterations;

    return n_iterations == args->in_max_iterations ? JMTX_RESULT_NOT_CONVERGED : JMTX_RESULT_SUCCESS;
}

/**
 * Uses Jacobi point iteration (also known as Jacobi method: https://en.wikipedia.org/wiki/Jacobi_method)
 * to solve the linear system Ax = y. Uses a relaxation factor ω for the following relation:
 *
 *  x(n + 1) = ω D⁻¹ (y - A x(n)) + x(n)
 *
 *  This may allow for better convergence
 *
 * @param mtx pointer to the memory where matrix A is stored as a compressed diagonal sparse matrix
 * @param y pointer to the memory where the vector y is stored
 * @param x pointer to the memory where the solution vector x will be stored
 * @param relaxation_factor relaxation factor used for the iterations, which must be greater than zero
 * @param aux_vec1 auxiliary memory for a vector of the same size as x and y
 * @param aux_vec2 auxiliary memory for a vector of the same size as x and y
 * @param args::in_convergence_criterion tolerance to determine if the solution is close enough
 * @param args::in_max_iterations number of iterations to stop at
 * @param args::out_last_error receives the value of the error criterion at the final iteration
 * @param args::out_last_iteration receives the number of the final iteration
 * @param args::opt_error_evolution (optional) pointer to an array of length max_iterations, that receives the error
 * value of each iteration
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_NOT_CONVERGED if it hasn't reached given stopping criterion,
 * in case of failure it returns the associated error code
 */
jmtx_result JMTX_NAME_TYPED(solve_iterative_jacobi_relaxed_cds)(
    const JMTX_NAME_TYPED(matrix_cds) * mtx, const JMTX_SCALAR_T *restrict y, JMTX_SCALAR_T *restrict x,
    JMTX_SCALAR_T relaxation_factor, JMTX_SCALAR_T *restrict aux_vec1, JMTX_SCALAR_T *restrict aux_vec2,
    JMTX_NAME_TYPED(solver_arguments) * args)
{
    //  Length of x and y
    const JMTX_INDEX_T n = mtx->base.cols;
    JMTX_SCALAR_T *div_factor = aux_vec1;
    JMTX_REAL_T y_mag = 0;
    //  Initial guess by assuming that mtx is a diagonal matrix
    for (JMTX_INDEX_T i = 0; i < n; ++i)
    {
        const JMTX_SCALAR_T d = mtx->main_diagonal[i];
        x[i] = y[i] / d;
        div_factor[i] = relaxation_factor / d;
        y_mag += JMTX_DOT(y[i], y[i]);
    }
    y_mag = JMTX_REAL_ROOT(y_mag);

    //  Memory used to store result of the current iteration
    JMTX_SCALAR_T *const auxiliary_x = aux_vec2;

    JMTX_SCALAR_T *x0 = auxiliary_x;
    JMTX_SCALAR_T *x1 = x;
    JMTX_REAL_T err;
    JMTX_INDEX_T n_iterations = 0;
    do
    {
        err = 0.0f;
        {
            JMTX_SCALAR_T *tmp = x1;
            x1 = x0;
            x0 = tmp;
        }

        //  For each entry, find the corresponding row in matrix A - D and compute the dot product between x and that
        //  row
        JMTX_NAME_TYPED(matrix_cds_vector_multiply)(mtx, x0, x1);
        for (JMTX_INDEX_T i = 0; i < n; ++i)
        {
            x1[i] = x0[i] + (y[i] - x1[i]) * div_factor[i];
        }

        JMTX_NAME_TYPED(matrix_cds_vector_multiply)(mtx, x1, x0);
        for (JMTX_INDEX_T i = 0; i < n; ++i)
        {
            const JMTX_SCALAR_T val = y[i] - x0[i];
            err += JMTX_DOT(val, val);
        }
        err = JMTX_REAL_ROOT(err) / (JMTX_REAL_T)y_mag;
        if (args->opt_error_evolution)
        {
            args->opt_error_evolution[n_iterations] = err;
        }
        n_iterations += 1;
    } while (err > args->in_convergence_criterion && n_iterations < args->in_max_iterations);

    if (x1 == auxiliary_x)
    {
        memcpy(x, auxiliary_x, sizeof *x * n);
    }

    args->out_last_iteration = n_iterations;
    args->out_last_error = err;

    return n_iterations == args->in_max_iterations ? JMTX_RESULT_NOT_CONVERGED : JMTX_RESULT_SUCCESS;
}

/**
 * Uses Jacobi point iteration (also known as Jacobi method: https://en.wikipedia.org/wiki/Jacobi_method)
 * to solve the linear system Ax = y
 *
 * @param mtx pointer to the memory where matrix A is stored as a band row major matrix
 * @param y pointer to the memory where the vector y is stored
 * @param x pointer to the memory where the solution vector x will be stored
 * @param aux_vec1 auxiliary memory for a vector of the same size as x and y
 * @param aux_vec2 auxiliary memory for a vector of the same size as x and y
 * @param args::in_convergence_criterion tolerance to determine if the solution is close enough
 * @param args::in_max_iterations number of iterations to stop at
 * @param args::out_last_error receives the value of the error criterion at the final iteration
 * @param args::out_last_iteration receives the number of the final iteration
 * @param args::opt_error_evolution (optional) pointer to an array of length max_iterations, that receives the error
 * value of each iteration
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_NOT_CONVERGED if it hasn't reached given stopping criterion,
 * in case of failure it returns the associated error code
 */
jmtx_result JMTX_NAME_TYPED(solve_iterative_jacobi_brm)(const JMTX_NAME_TYPED(matrix_brm) * mtx,
                                                        const JMTX_SCALAR_T *restrict y, JMTX_SCALAR_T *restrict x,
                                                        JMTX_SCALAR_T *restrict aux_vec1,
                                                        JMTX_SCALAR_T *restrict aux_vec2,
                                                        JMTX_NAME_TYPED(solver_arguments) * args)
{
    //  Length of x and y
    const JMTX_INDEX_T n = mtx->base.cols;
    JMTX_SCALAR_T *const div_factor = aux_vec1;
    JMTX_REAL_T y_mag = 0;
    //  Initial guess by assuming that mtx is a diagonal matrix
    for (JMTX_INDEX_T i = 0; i < n; ++i)
    {
        const JMTX_SCALAR_T d = JMTX_NAME_TYPED(matrix_brm_get_entry)(mtx, i, i);
        x[i] = y[i] / d;
        div_factor[i] = 1.0f / d;
        y_mag += JMTX_DOT(y[i], y[i]);
    }
    y_mag = JMTX_REAL_ROOT(y_mag);

    //  Memory used to store result of the current iteration
    JMTX_SCALAR_T *const auxiliary_x = aux_vec2;

    JMTX_SCALAR_T *x0 = auxiliary_x;
    JMTX_SCALAR_T *x1 = x;
    JMTX_REAL_T err;
    JMTX_INDEX_T n_iterations = 0;
    do
    {
        err = 0.0f;
        {
            JMTX_SCALAR_T *tmp = x1;
            x1 = x0;
            x0 = tmp;
        }

        //  For each entry, find the corresponding row in matrix A - D and compute the dot product between x and that
        //  row
        JMTX_NAME_TYPED(matrix_brm_vector_multiply)(mtx, x0, x1);
        for (JMTX_INDEX_T i = 0; i < n; ++i)
        {
            //  Multiplication of vector x by D⁻¹
            x1[i] = (y[i] - x1[i]) * div_factor[i];
        }

        JMTX_NAME_TYPED(matrix_brm_vector_multiply)(mtx, x1, x0);
        for (JMTX_INDEX_T i = 0; i < n; ++i)
        {
            const JMTX_SCALAR_T val = y[i] - x0[i];
            err += JMTX_DOT(val, val);
        }
        //  Have "err" as ratio between magnitude of y vector and magnitude of residual
        err = JMTX_REAL_ROOT(err) / (JMTX_REAL_T)y_mag;
        if (args->opt_error_evolution)
        {
            args->opt_error_evolution[n_iterations] = err;
        }
        n_iterations += 1;
    } while (err > args->in_convergence_criterion && n_iterations < args->in_max_iterations);

    if (x1 == auxiliary_x)
    {
        memcpy(x, auxiliary_x, sizeof *x * n);
    }

    args->out_last_error = err;
    args->out_last_iteration = n_iterations;

    return n_iterations == args->in_max_iterations ? JMTX_RESULT_NOT_CONVERGED : JMTX_RESULT_SUCCESS;
}

/**
 * Uses Jacobi point iteration (also known as Jacobi method: https://en.wikipedia.org/wiki/Jacobi_method)
 * to solve the linear system Ax = y. Uses a relaxation factor ω for the following relation:
 *
 *  x(n + 1) = ω D⁻¹ (y - A x(n)) + x(n)
 *
 *  This may allow for better convergence
 *
 * @param mtx pointer to the memory where matrix A is stored as a band row major matrix
 * @param y pointer to the memory where the vector y is stored
 * @param x pointer to the memory where the solution vector x will be stored
 * @param aux_vec1 auxiliary memory for a vector of the same size as x and y
 * @param aux_vec2 auxiliary memory for a vector of the same size as x and y
 * @param relaxation_factor relaxation factor used for the iterations, which must be greater than zero
 * @param args::in_convergence_criterion tolerance to determine if the solution is close enough
 * @param args::in_max_iterations number of iterations to stop at
 * @param args::out_last_error receives the value of the error criterion at the final iteration
 * @param args::out_last_iteration receives the number of the final iteration
 * @param args::opt_error_evolution (optional) pointer to an array of length max_iterations, that receives the error
 * value of each iteration
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_NOT_CONVERGED if it hasn't reached given stopping criterion,
 * in case of failure it returns the associated error code
 */
jmtx_result JMTX_NAME_TYPED(solve_iterative_jacobi_relaxed_brm)(
    const JMTX_NAME_TYPED(matrix_brm) * mtx, const JMTX_SCALAR_T *restrict y, JMTX_SCALAR_T *restrict x,
    JMTX_SCALAR_T relaxation_factor, JMTX_SCALAR_T *restrict aux_vec1, JMTX_SCALAR_T *restrict aux_vec2,
    JMTX_NAME_TYPED(solver_arguments) * args)
{
    //  Length of x and y
    const JMTX_INDEX_T n = mtx->base.cols;
    JMTX_SCALAR_T *div_factor = aux_vec1;
    JMTX_REAL_T y_mag = 0;
    //  Initial guess by assuming that mtx is a diagonal matrix
    for (JMTX_INDEX_T i = 0; i < n; ++i)
    {
        const JMTX_SCALAR_T d = JMTX_NAME_TYPED(matrix_brm_get_entry)(mtx, i, i);
        x[i] = y[i] / d;
        div_factor[i] = relaxation_factor / d;
        y_mag += JMTX_DOT(y[i], y[i]);
    }
    y_mag = JMTX_REAL_ROOT(y_mag);

    //  Memory used to store result of the current iteration
    JMTX_SCALAR_T *const auxiliary_x = aux_vec2;

    JMTX_SCALAR_T *x0 = auxiliary_x;
    JMTX_SCALAR_T *x1 = x;
    JMTX_REAL_T err;
    JMTX_INDEX_T n_iterations = 0;
    do
    {
        err = 0.0f;
        {
            JMTX_SCALAR_T *tmp = x1;
            x1 = x0;
            x0 = tmp;
        }

        //  For each entry, find the corresponding row in matrix A - D and compute the dot product between x and that
        //  row
        JMTX_NAME_TYPED(matrix_brm_vector_multiply)(mtx, x0, x1);
        for (JMTX_INDEX_T i = 0; i < n; ++i)
        {
            x1[i] = x0[i] + (y[i] - x1[i]) * div_factor[i];
        }

        JMTX_NAME_TYPED(matrix_brm_vector_multiply)(mtx, x1, x0);
        for (JMTX_INDEX_T i = 0; i < n; ++i)
        {
            const JMTX_SCALAR_T val = y[i] - x0[i];
            err += JMTX_DOT(val, val);
        }
        err = JMTX_REAL_ROOT(err) / (JMTX_REAL_T)y_mag;
        if (args->opt_error_evolution)
        {
            args->opt_error_evolution[n_iterations] = err;
        }
        n_iterations += 1;
    } while (err > args->in_convergence_criterion && n_iterations < args->in_max_iterations);

    if (x1 == auxiliary_x)
    {
        memcpy(x, auxiliary_x, sizeof *x * n);
    }

    args->out_last_iteration = n_iterations;
    args->out_last_error = err;

    return n_iterations == args->in_max_iterations ? JMTX_RESULT_NOT_CONVERGED : JMTX_RESULT_SUCCESS;
}
/**
 * Uses Jacobi point iteration (also known as Jacobi method: https://en.wikipedia.org/wiki/Jacobi_method)
 * to solve the linear system Ax = y.
 *
 * This version of the function uses OpenMP to solve the system in parallel.
 *
 * @param mtx pointer to the memory where matrix A is stored as a band row major matrix
 * @param y pointer to the memory where the vector y is stored
 * @param x pointer to the memory where the solution vector x will be stored
 * @param aux_vec1 auxiliary memory for a vector of the same size as x and y
 * @param aux_vec2 auxiliary memory for a vector of the same size as x and y
 * @param args::in_convergence_criterion tolerance to determine if the solution is close enough
 * @param args::in_max_iterations number of iterations to stop at
 * @param args::out_last_error receives the value of the error criterion at the final iteration
 * @param args::out_last_iteration receives the number of the final iteration
 * @param args::opt_error_evolution (optional) pointer to an array of length max_iterations, that receives the error
 * value of each iteration
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_NOT_CONVERGED if it hasn't reached given stopping criterion,
 * in case of failure it returns the associated error code
 */
jmtx_result JMTX_NAME_TYPED(solve_iterative_jacobi_brm_parallel)(
    const JMTX_NAME_TYPED(matrix_brm) * mtx, const JMTX_SCALAR_T *restrict y, JMTX_SCALAR_T *restrict x,
    JMTX_SCALAR_T *restrict aux_vector1, JMTX_SCALAR_T *restrict aux_vector2, JMTX_NAME_TYPED(solver_arguments) * args)
{
    //  Length of x and y
    const JMTX_INDEX_T n = mtx->base.cols;
    JMTX_REAL_T y_mag = 0;
    //  Memory used to store result of the current iteration
    JMTX_SCALAR_T *x0 = aux_vector2;
    JMTX_SCALAR_T *x1 = x;
    JMTX_REAL_T err = 0;
    JMTX_INDEX_T n_iterations = 0;
    JMTX_REAL_T *const p_error_evolution = args->opt_error_evolution;
    const JMTX_REAL_T convergence_dif = args->in_convergence_criterion;
    const JMTX_INDEX_T max_iterations = args->in_max_iterations;

#pragma omp parallel default(none) shared(err, x0, x1, y, mtx, aux_vector1, y_mag, p_error_evolution, n_iterations,    \
                                              convergence_dif, max_iterations, n)
    {
#pragma omp for schedule(static)
        for (JMTX_INDEX_T i = 0; i < n; ++i)
        {
            const JMTX_SCALAR_T d = JMTX_NAME_TYPED(matrix_brm_get_entry)(mtx, i, i);
            aux_vector1[i] = 1.0f / d;
            const JMTX_SCALAR_T mag = JMTX_DOT(y[i], y[i]);
            y_mag += mag;
        }

#pragma omp master
        {
            y_mag = JMTX_REAL_ROOT(y_mag);
        }

        do
        {
#pragma omp barrier
#pragma omp master
            {
                err = 0.0f;
                JMTX_SCALAR_T *tmp = x1;
                x1 = x0;
                x0 = tmp;
            }
#pragma omp barrier

            //  For each entry, find the corresponding row in matrix A - D and compute the dot product between x and
            //  that row
#pragma omp for schedule(static)
            for (JMTX_INDEX_T i = 0; i < n; ++i)
            {
                x1[i] = x0[i] + (y[i] - JMTX_NAME_TYPED(matrix_brm_vector_multiply_row)(mtx, x0, i)) * aux_vector1[i];
            }

#pragma omp for reduction(+ : err) schedule(static)
            for (JMTX_INDEX_T i = 0; i < n; ++i)
            {
                const JMTX_SCALAR_T val = JMTX_NAME_TYPED(matrix_brm_vector_multiply_row)(mtx, x1, i) - y[i];
                err += JMTX_DOT(val, val);
            }

#pragma omp master
            {
                err = JMTX_REAL_ROOT(err) / (JMTX_REAL_T)y_mag;
                if (p_error_evolution != NULL)
                {
                    p_error_evolution[n_iterations] = err;
                }
                n_iterations += 1;
            }
#pragma omp barrier
        } while (err > convergence_dif && n_iterations < max_iterations);
    }

    if (x1 == aux_vector2)
    {
        memcpy(x, aux_vector2, sizeof *x * n);
    }
    args->out_last_iteration = n_iterations;
    args->out_last_error = err;

    return n_iterations == max_iterations || !isfinite(err) ? JMTX_RESULT_NOT_CONVERGED : JMTX_RESULT_SUCCESS;
}
