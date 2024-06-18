// Automatically generated from source/float/solvers/jacobi_point_iteration.c on Thu Dec 14 17:26:53 2023
//
// Created by jan on 15.6.2022.
//

#include "../matrices/sparse_row_compressed_internal.h"
#include "../matrices/sparse_diagonal_compressed_internal.h"
#include "../matrices/band_row_major_internal.h"
#include "../../../include/jmtx/double/solvers/jacobi_point_iteration.h"
#include <math.h>
#include <stdio.h>
#include <assert.h>

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
 * @param args::opt_error_evolution (optional) pointer to an array of length max_iterations, that receives the error value of each
 * iteration
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_NOT_CONVERGED if it hasn't reached given stopping criterion,
 * in case of failure it returns the associated error code
 */
jmtx_result jmtxd_solve_iterative_jacobi_crs(
        const jmtxd_matrix_crs* mtx, const double* restrict y, double* restrict x, double* restrict aux_vec1, double* restrict aux_vec2,
        jmtxd_solver_arguments* args)
{
    //  Length of x and y
    const uint32_t n = mtx->base.cols;
    double* const div_factor = aux_vec1;
    double y_mag = 0;
    //  Initial guess by assuming that mtx is a diagonal matrix
    for (uint32_t i = 0; i < n; ++i)
    {
        double d = jmtxd_matrix_crs_get_entry(mtx, i, i);
        x[i] = y[i] / d;
        div_factor[i] = 1.0f / d;
        y_mag += y[i] * y[i];
    }
    y_mag = sqrt(y_mag);

    //  Memory used to store result of the current iteration
    double* const auxiliary_x = aux_vec2;

    double* x0 = auxiliary_x;
    double* x1 = x;
    double err;
    uint32_t n_iterations = 0;
    do
    {
        err = 0.0f;
        {
            double* tmp = x1;
            x1 = x0;
            x0 = tmp;
        }

        //  For each entry, find the corresponding row in matrix A - D and compute the dot product between x and that row
        for (uint32_t i = 0; i < n; ++i)
        {
            double* row_ptr;
            uint32_t* index_ptr;
            uint32_t n_elements = jmtxd_matrix_crs_get_row(mtx, i, &index_ptr, &row_ptr);
            double res = 0;
            for (uint32_t j = 0; j < n_elements; ++j)
            {
                if (i != index_ptr[j])
                {
                    res += row_ptr[j] * x0[index_ptr[j]];
                }
            }
            //  Multiplication of vector x by D⁻¹
            x1[i] = (y[i] - res) * div_factor[i];
        }

        for (uint32_t i = 0; i < n; ++i)
        {
            double val = jmtxd_matrix_crs_vector_multiply_row(mtx, x1, i);
            val -= y[i];
            err += val * val;
        }
        //  Have "err" as ratio between magnitude of y vector and magnitude of residual
        err = sqrt(err) / (double)y_mag;
        if (args->opt_error_evolution)
        {
            args->opt_error_evolution[n_iterations] = err;
        }
        n_iterations += 1;
    } while(err > args->in_convergence_criterion && n_iterations < args->in_max_iterations);


    if (x1 == auxiliary_x)
    {
        memcpy(x, auxiliary_x, sizeof*x * n);
    }

    args->out_last_error = err;
    args->out_last_iteration = n_iterations;
//    LEAVE_FUNCTION();
    return n_iterations == args->in_max_iterations ? JMTX_RESULT_NOT_CONVERGED : JMTX_RESULT_SUCCESS;
}

/**
 * Uses Jacobi point iteration (also known as Jacobi method: https://en.wikipedia.org/wiki/Jacobi_method)
 * to solve the linear system Ax = y
 *
 * @param mtx pointer to the memory where matrix A is stored as a compressed row sparse matrix
 * @param n size of the matrix and input vectors
 * @param y pointer to the memory where the vector y is stored
 * @param x pointer to the memory where the solution vector x will be stored
 * @param aux_vec1 auxiliary memory for a vector of the same size as x and y
 * @param aux_vec2 auxiliary memory for a vector of the same size as x and y
 * @param args::in_convergence_criterion tolerance to determine if the solution is close enough
 * @param args::in_max_iterations number of iterations to stop at
 * @param args::out_last_error receives the value of the error criterion at the final iteration
 * @param args::out_last_iteration receives the number of the final iteration
 * @param args::opt_error_evolution (optional) pointer to an array of length max_iterations, that receives the error value of each
 * iteration
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_NOT_CONVERGED if it hasn't reached given stopping criterion,
 * in case of failure it returns the associated error code
 */
jmtx_result jmtxds_solve_iterative_jacobi_crs(
        const jmtxd_matrix_crs* mtx, const uint32_t n, const double y[static restrict n], double x[restrict n], double aux_vec1[restrict n], double aux_vec2[restrict n],
        jmtxd_solver_arguments* args)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.rows != n || mtx->base.cols != n)
    {
        return JMTX_RESULT_BAD_MATRIX;
    }
    if (mtx->base.type != JMTXD_TYPE_CRS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (!x)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!args)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!aux_vec1)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!aux_vec2)
    {
        return JMTX_RESULT_NULL_PARAM;
    }



    //  Length of x and y
    double* const div_factor = aux_vec1;
    double y_mag = 0;
    //  Initial guess by assuming that mtx is a diagonal matrix
    for (uint32_t i = 0; i < n; ++i)
    {
        double d = jmtxd_matrix_crs_get_entry(mtx, i, i);
        if (d == 0.0f)
        {
            //  Diagonal entry is zero!
            //  Can't solve this one with Jacobi
            return JMTX_RESULT_BAD_MATRIX;
        }
        x[i] = y[i] / d;
        div_factor[i] = 1.0f / d;
        y_mag += y[i] * y[i];
    }
    y_mag = sqrt(y_mag);

    //  Memory used to store result of the current iteration
    double* const auxiliary_x = aux_vec2;

    double* x0 = auxiliary_x;
    double* x1 = x;
    double err;
    uint_fast32_t n_iterations = 0;
    do
    {
        err = 0.0f;
        {
            double* tmp = x1;
            x1 = x0;
            x0 = tmp;
        }

        //  For each entry, find the corresponding row in matrix A - D and compute the dot product between x and that row
        for (uint_fast32_t i = 0; i < n; ++i)
        {
            const double res = jmtxd_matrix_crs_vector_multiply_row(mtx, x0, i);
            //  Multiplication of vector x by D⁻¹
            x1[i] = (y[i] - res) * div_factor[i];
        }

        //  Vector x0 no longer needed
        jmtxd_matrix_crs_vector_multiply(mtx, x1, x0);
        //  Loop is easier to vectorize like this
        for (uint_fast32_t i = 0; i < n; ++i)
        {
            const double r = y[i] - x0[i];
            err += r * r;
        }
        //  Have "err" as ratio between magnitude of y vector and magnitude of residual
        err = sqrt(err) / (double)y_mag;
        if (args->opt_error_evolution)
        {
            args->opt_error_evolution[n_iterations] = err;
        }
        n_iterations += 1;
    } while(err > args->in_convergence_criterion && n_iterations < args->in_max_iterations);


    if (x1 == auxiliary_x)
    {
        memcpy(x, auxiliary_x, sizeof*x * n);
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
 * @param args::opt_error_evolution (optional) pointer to an array of length max_iterations, that receives the error value of each
 * iteration
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_NOT_CONVERGED if it hasn't reached given stopping criterion,
 * in case of failure it returns the associated error code
 */
jmtx_result jmtxd_solve_iterative_jacobi_relaxed_crs(
        const jmtxd_matrix_crs* mtx, const double* restrict y, double* restrict x, double relaxation_factor, double* restrict aux_vec1,
        double* restrict aux_vec2, jmtxd_solver_arguments* args)
{
    //  Length of x and y
    const uint32_t n = mtx->base.cols;
    double* div_factor = aux_vec1;
    double y_mag = 0;
    //  Initial guess by assuming that mtx is a diagonal matrix
    for (uint32_t i = 0; i < n; ++i)
    {
        const double d = jmtxd_matrix_crs_get_entry(mtx, i, i);
        x[i] = y[i] / d;
        div_factor[i] = relaxation_factor / d;
        y_mag += y[i] * y[i];
    }
    y_mag = sqrt(y_mag);

    //  Memory used to store result of the current iteration
    double* const auxiliary_x = aux_vec2;

    double* x0 = auxiliary_x;
    double* x1 = x;
    double err;
    uint32_t n_iterations = 0;
    do
    {
        err = 0.0f;
        {
            double* tmp = x1;
            x1 = x0;
            x0 = tmp;
        }

        //  For each entry, find the corresponding row in matrix A - D and compute the dot product between x and that row
        for (uint32_t i = 0; i < n; ++i)
        {
            const double res = jmtxd_matrix_crs_vector_multiply_row(mtx, x0, i);
            x1[i] = x0[i] + (y[i] - res) * div_factor[i];
        }

        jmtxd_matrix_crs_vector_multiply(mtx, x1, x0);
        for (uint32_t i = 0; i < n; ++i)
        {
            const double val = y[i] - x0[i];
            err += val * val;
        }
        err = sqrt(err) / (double)y_mag;
        if (args->opt_error_evolution)
        {
            args->opt_error_evolution[n_iterations] = err;
        }
        n_iterations += 1;
    } while(err > args->in_convergence_criterion && n_iterations < args->in_max_iterations);


    if (x1 == auxiliary_x)
    {
        memcpy(x, auxiliary_x, sizeof*x * n);
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
 * @param args::opt_error_evolution (optional) pointer to an array of length max_iterations, that receives the error value of each
 * iteration
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_NOT_CONVERGED if it hasn't reached given stopping criterion,
 * in case of failure it returns the associated error code
 */
jmtx_result jmtxd_solve_iterative_jacobi_crs_parallel(
        const jmtxd_matrix_crs* mtx, const double* restrict y, double* restrict x, double* restrict aux_vector1, double* restrict aux_vector2,
        jmtxd_solver_arguments* args)
{
    //  Length of x and y
    const uint32_t n = mtx->base.cols;
    double y_mag = 0;
    //  Memory used to store result of the current iteration
    double* x0 = aux_vector2;
    double* x1 = x;
    double err = 0.0f;
    uint32_t n_iterations = 0;
    double* const p_error_evolution = args->opt_error_evolution;
    const double convergence_dif = args->in_convergence_criterion;
    const uint32_t max_iterations = args->in_max_iterations;
#pragma omp parallel default(none) shared(err, x0, x1, y, mtx, aux_vector1, y_mag, p_error_evolution, n_iterations,\
    convergence_dif, max_iterations, n)
    {
#pragma omp for schedule(static)
        for (uint32_t i = 0; i < n; ++i)
        {
            const double d = jmtxd_matrix_crs_get_entry(mtx, i, i);
            aux_vector1[i] = 1.0f / d;
            const double mag = y[i] * y[i];
            y_mag += mag;
        }

#pragma omp master
        {
            y_mag = sqrt(y_mag);
        }

        do
        {
#pragma omp barrier
#pragma omp master
            {
                err = 0.0f;
                double* tmp = x1;
                x1 = x0;
                x0 = tmp;
            }
#pragma omp barrier

            //  For each entry, find the corresponding row in matrix A - D and compute the dot product between x and that row
#pragma omp for schedule(static)
            for (uint32_t i = 0; i < n; ++i)
            {
                x1[i] = x0[i] + (y[i] - jmtxd_matrix_crs_vector_multiply_row(mtx, x0, i)) * aux_vector1[i];
            }

#pragma omp for reduction(+:err) schedule(static)
            for (uint32_t i = 0; i < n; ++i)
            {
                const double val = jmtxd_matrix_crs_vector_multiply_row(mtx, x1, i) - y[i];
                err += val * val;
            }


#pragma omp master
            {
                err = sqrt(err) / (double) y_mag;
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
        memcpy(x, aux_vector2, sizeof*x * n);
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
 * @param args::opt_error_evolution (optional) pointer to an array of length max_iterations, that receives the error value of each
 * iteration
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_NOT_CONVERGED if it hasn't reached given stopping criterion,
 * in case of failure it returns the associated error code
 */
jmtx_result jmtxd_solve_iterative_jacobi_cds(
        const jmtxd_matrix_cds* mtx, const double* restrict y, double* restrict x, double* restrict aux_vec1, double* restrict aux_vec2,
        jmtxd_solver_arguments* args)
{
    //  Length of x and y
    const uint32_t n = mtx->base.cols;
    double* const div_factor = aux_vec1;
    double y_mag = 0;
    //  Initial guess by assuming that mtx is a diagonal matrix
    for (uint32_t i = 0; i < n; ++i)
    {
        double d = mtx->main_diagonal[i];
        x[i] = y[i] / d;
        div_factor[i] = 1.0f / d;
        y_mag += y[i] * y[i];
    }
    y_mag = sqrt(y_mag);

    //  Memory used to store result of the current iteration
    double* const auxiliary_x = aux_vec2;

    double* x0 = auxiliary_x;
    double* x1 = x;
    double err;
    uint32_t n_iterations = 0;
    do
    {
        err = 0.0f;
        {
            double* tmp = x1;
            x1 = x0;
            x0 = tmp;
        }

        //  For each entry, find the corresponding row in matrix A - D and compute the dot product between x and that row
        jmtxd_matrix_cds_vector_multiply(mtx, x0, x1);
        for (uint32_t i = 0; i < n; ++i)
        {
            //  Multiplication of vector x by D⁻¹
            x1[i] = (y[i] - x1[i]) * div_factor[i];
        }

        jmtxd_matrix_cds_vector_multiply(mtx, x1, x0);
        for (uint32_t i = 0; i < n; ++i)
        {
            const double val = y[i] - x0[i];
            err += val * val;
        }
        //  Have "err" as ratio between magnitude of y vector and magnitude of residual
        err = sqrt(err) / (double)y_mag;
        if (args->opt_error_evolution)
        {
            args->opt_error_evolution[n_iterations] = err;
        }
        n_iterations += 1;
    } while(err > args->in_convergence_criterion && n_iterations < args->in_max_iterations);


    if (x1 == auxiliary_x)
    {
        memcpy(x, auxiliary_x, sizeof*x * n);
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
 * @param args::opt_error_evolution (optional) pointer to an array of length max_iterations, that receives the error value of each
 * iteration
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_NOT_CONVERGED if it hasn't reached given stopping criterion,
 * in case of failure it returns the associated error code
 */
jmtx_result jmtxd_solve_iterative_jacobi_relaxed_cds(
        const jmtxd_matrix_cds* mtx, const double* restrict y, double* restrict x, double relaxation_factor, double* restrict aux_vec1,
        double* restrict aux_vec2, jmtxd_solver_arguments* args)
{
    //  Length of x and y
    const uint32_t n = mtx->base.cols;
    double* div_factor = aux_vec1;
    double y_mag = 0;
    //  Initial guess by assuming that mtx is a diagonal matrix
    for (uint32_t i = 0; i < n; ++i)
    {
        const double d = mtx->main_diagonal[i];
        x[i] = y[i] / d;
        div_factor[i] = relaxation_factor / d;
        y_mag += y[i] * y[i];
    }
    y_mag = sqrt(y_mag);

    //  Memory used to store result of the current iteration
    double* const auxiliary_x = aux_vec2;

    double* x0 = auxiliary_x;
    double* x1 = x;
    double err;
    uint32_t n_iterations = 0;
    do
    {
        err = 0.0f;
        {
            double* tmp = x1;
            x1 = x0;
            x0 = tmp;
        }

        //  For each entry, find the corresponding row in matrix A - D and compute the dot product between x and that row
        jmtxd_matrix_cds_vector_multiply(mtx, x0, x1);
        for (uint32_t i = 0; i < n; ++i)
        {
            x1[i] = x0[i] + (y[i] - x1[i]) * div_factor[i];
        }

        jmtxd_matrix_cds_vector_multiply(mtx, x1, x0);
        for (uint32_t i = 0; i < n; ++i)
        {
            const double val = y[i] - x0[i];
            err += val * val;
        }
        err = sqrt(err) / (double)y_mag;
        if (args->opt_error_evolution)
        {
            args->opt_error_evolution[n_iterations] = err;
        }
        n_iterations += 1;
    } while(err > args->in_convergence_criterion && n_iterations < args->in_max_iterations);


    if (x1 == auxiliary_x)
    {
        memcpy(x, auxiliary_x, sizeof*x * n);
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
 * @param args::opt_error_evolution (optional) pointer to an array of length max_iterations, that receives the error value of each
 * iteration
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_NOT_CONVERGED if it hasn't reached given stopping criterion,
 * in case of failure it returns the associated error code
 */
jmtx_result jmtxd_solve_iterative_jacobi_brm(
        const jmtxd_matrix_brm* mtx, const double* restrict y, double* restrict x, double* restrict aux_vec1, double* restrict aux_vec2,
        jmtxd_solver_arguments* args)
{
    //  Length of x and y
    const uint32_t n = mtx->base.cols;
    double* const div_factor = aux_vec1;
    double y_mag = 0;
    //  Initial guess by assuming that mtx is a diagonal matrix
    for (uint32_t i = 0; i < n; ++i)
    {
        const double d = jmtxd_matrix_brm_get_entry(mtx, i, i);
        x[i] = y[i] / d;
        div_factor[i] = 1.0f / d;
        y_mag += y[i] * y[i];
    }
    y_mag = sqrt(y_mag);

    //  Memory used to store result of the current iteration
    double* const auxiliary_x = aux_vec2;

    double* x0 = auxiliary_x;
    double* x1 = x;
    double err;
    uint32_t n_iterations = 0;
    do
    {
        err = 0.0f;
        {
            double* tmp = x1;
            x1 = x0;
            x0 = tmp;
        }

        //  For each entry, find the corresponding row in matrix A - D and compute the dot product between x and that row
        jmtxd_matrix_brm_vector_multiply(mtx, x0, x1);
        for (uint32_t i = 0; i < n; ++i)
        {
            //  Multiplication of vector x by D⁻¹
            x1[i] = (y[i] - x1[i]) * div_factor[i];
        }

        jmtxd_matrix_brm_vector_multiply(mtx, x1, x0);
        for (uint32_t i = 0; i < n; ++i)
        {
            const double val = y[i] - x0[i];
            err += val * val;
        }
        //  Have "err" as ratio between magnitude of y vector and magnitude of residual
        err = sqrt(err) / (double)y_mag;
        if (args->opt_error_evolution)
        {
            args->opt_error_evolution[n_iterations] = err;
        }
        n_iterations += 1;
    } while(err > args->in_convergence_criterion && n_iterations < args->in_max_iterations);


    if (x1 == auxiliary_x)
    {
        memcpy(x, auxiliary_x, sizeof*x * n);
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
 * @param args::opt_error_evolution (optional) pointer to an array of length max_iterations, that receives the error value of each
 * iteration
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_NOT_CONVERGED if it hasn't reached given stopping criterion,
 * in case of failure it returns the associated error code
 */
jmtx_result jmtxd_solve_iterative_jacobi_relaxed_brm(
        const jmtxd_matrix_brm* mtx, const double* restrict y, double* restrict x, double relaxation_factor, double* restrict aux_vec1,
        double* restrict aux_vec2, jmtxd_solver_arguments* args)
{
    //  Length of x and y
    const uint32_t n = mtx->base.cols;
    double* div_factor = aux_vec1;
    double y_mag = 0;
    //  Initial guess by assuming that mtx is a diagonal matrix
    for (uint32_t i = 0; i < n; ++i)
    {
        const double d = jmtxd_matrix_brm_get_entry(mtx, i, i);
        x[i] = y[i] / d;
        div_factor[i] = relaxation_factor / d;
        y_mag += y[i] * y[i];
    }
    y_mag = sqrt(y_mag);

    //  Memory used to store result of the current iteration
    double* const auxiliary_x = aux_vec2;

    double* x0 = auxiliary_x;
    double* x1 = x;
    double err;
    uint32_t n_iterations = 0;
    do
    {
        err = 0.0f;
        {
            double* tmp = x1;
            x1 = x0;
            x0 = tmp;
        }

        //  For each entry, find the corresponding row in matrix A - D and compute the dot product between x and that row
        jmtxd_matrix_brm_vector_multiply(mtx, x0, x1);
        for (uint32_t i = 0; i < n; ++i)
        {
            x1[i] = x0[i] + (y[i] - x1[i]) * div_factor[i];
        }

        jmtxd_matrix_brm_vector_multiply(mtx, x1, x0);
        for (uint32_t i = 0; i < n; ++i)
        {
            const double val = y[i] - x0[i];
            err += val * val;
        }
        err = sqrt(err) / (double)y_mag;
        if (args->opt_error_evolution)
        {
            args->opt_error_evolution[n_iterations] = err;
        }
        n_iterations += 1;
    } while(err > args->in_convergence_criterion && n_iterations < args->in_max_iterations);


    if (x1 == auxiliary_x)
    {
        memcpy(x, auxiliary_x, sizeof*x * n);
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
 * @param args::opt_error_evolution (optional) pointer to an array of length max_iterations, that receives the error value of each
 * iteration
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_NOT_CONVERGED if it hasn't reached given stopping criterion,
 * in case of failure it returns the associated error code
 */
jmtx_result jmtxd_solve_iterative_jacobi_brm_parallel(
        const jmtxd_matrix_brm* mtx, const double* restrict y, double* restrict x, double* restrict aux_vector1, double* restrict aux_vector2,
        jmtxd_solver_arguments* args)
{
    //  Length of x and y
    const uint32_t n = mtx->base.cols;
    double y_mag = 0;
    //  Memory used to store result of the current iteration
    double* x0 = aux_vector2;
    double* x1 = x;
    double err = 0.0f;
    uint32_t n_iterations = 0;
    double* const p_error_evolution = args->opt_error_evolution;
    const double convergence_dif = args->in_convergence_criterion;
    const uint32_t max_iterations = args->in_max_iterations;
#pragma omp parallel default(none) shared(err, x0, x1, y, mtx, aux_vector1, y_mag, p_error_evolution, n_iterations,\
    convergence_dif, max_iterations, n)
    {
#pragma omp for schedule(static)
        for (uint32_t i = 0; i < n; ++i)
        {
            const double d = jmtxd_matrix_brm_get_entry(mtx, i, i);
            aux_vector1[i] = 1.0f / d;
            const double mag = y[i] * y[i];
            y_mag += mag;
        }

#pragma omp master
        {
            y_mag = sqrt(y_mag);
        }

        do
        {
#pragma omp barrier
#pragma omp master
            {
                err = 0.0f;
                double* tmp = x1;
                x1 = x0;
                x0 = tmp;
            }
#pragma omp barrier

            //  For each entry, find the corresponding row in matrix A - D and compute the dot product between x and that row
#pragma omp for schedule(static)
            for (uint32_t i = 0; i < n; ++i)
            {
                x1[i] = x0[i] + (y[i] - jmtxd_matrix_brm_vector_multiply_row(mtx, x0, i)) * aux_vector1[i];
            }

#pragma omp for reduction(+:err) schedule(static)
            for (uint32_t i = 0; i < n; ++i)
            {
                const double val = jmtxd_matrix_brm_vector_multiply_row(mtx, x1, i) - y[i];
                err += val * val;
            }


#pragma omp master
            {
                err = sqrt(err) / (double) y_mag;
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
        memcpy(x, aux_vector2, sizeof*x * n);
    }
    args->out_last_iteration = n_iterations;
    args->out_last_error = err;

    return n_iterations == max_iterations || !isfinite(err) ? JMTX_RESULT_NOT_CONVERGED : JMTX_RESULT_SUCCESS;
}
