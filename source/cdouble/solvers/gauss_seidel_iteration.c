// Automatically generated from source/cfloat/solvers/gauss_seidel_iteration.c on Fri Dec  1 18:48:06 2023
// Automatically generated from source/cdouble/solvers/gauss_seidel_iteration.c on Fri Dec  1 17:36:03 2023
//
// Created by jan on 16.6.2022.
//

#include "../../../include/jmtx/cdouble/solvers/gauss_seidel_iteration.h"
#include "../matrices/sparse_row_compressed_internal.h"

#include <math.h>
#include <complex.h>


jmtx_result jmtxz_gauss_seidel_crs(const jmtxz_matrix_crs* mtx, const _Complex double* restrict y, _Complex double* restrict x, jmtxd_solver_arguments* args)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTXZ_TYPE_CRS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (!y)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!x)
    {
        return JMTX_RESULT_NULL_PARAM;
    }

    if (!args)
    {
        return JMTX_RESULT_NULL_PARAM;
    }

    //  Length of x and y
    const uint32_t n = mtx->base.cols;

    //  Initial guess of x is just y, which would be the case if matrix is just I
    memcpy(x, y, n * sizeof *x);
    //  Improve the guess by assuming that mtx is a diagonal matrix
    for (uint32_t i = 0; i < n; ++i)
    {
        _Complex double d = jmtxz_matrix_crs_get_entry(mtx, i, i);
        x[i] /= d;
    }


    double err;
    uint32_t n_iterations = 0;
    do
    {

        for (uint32_t i = 0; i < n; ++i)
        {
            _Complex double* row_ptr;
            uint32_t* index_ptr;
            uint32_t n_elements = jmtxz_matrix_crs_get_row(mtx, i, &index_ptr, &row_ptr);
            _Complex double res = (_Complex double)0.0;
            uint32_t k = 0;
            for (uint32_t j = 0; j < n_elements; ++j)
            {
                if (i != index_ptr[j])
                {
                    res += row_ptr[j] * x[index_ptr[j]];
                }
                else
                {
                    k = j;
                }
            }
            const _Complex double new_x = (y[i] - res) / row_ptr[k];
            x[i] = new_x;
        }

        err = 0;
        for (uint32_t i = 0; i < n; ++i)
        {
            _Complex double val;
            val = jmtxz_matrix_crs_vector_multiply_row(mtx, x, i);
            val -= y[i];
            err += val * val;
        }
        err = sqrt(err) / (double)n;
        if (args->opt_error_evolution)
        {
            args->opt_error_evolution[n_iterations] = err;
        }

        n_iterations += 1;
    } while(err > args->in_convergence_criterion && n_iterations < args->in_max_iterations);


    args->out_last_iteration = n_iterations;
    args->out_last_error = err;

    return n_iterations == args->in_max_iterations ? JMTX_RESULT_NOT_CONVERGED : JMTX_RESULT_SUCCESS;
}

jmtx_result jmtxz_gauss_seidel_crs_parallel(
        const jmtxz_matrix_crs* const mtx, const _Complex double* restrict y, _Complex double* restrict x, _Complex double* restrict aux_vector,
        jmtxd_solver_arguments* args)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTXZ_TYPE_CRS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (mtx->base.cols != mtx->base.rows)
    {
        //  Only doing square matrices
        return JMTX_RESULT_BAD_MATRIX;
    }
    if (!y)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!x)
    {
//        REPORT_ERROR_MESSAGE("Vector x pointer was null");
//        LEAVE_FUNCTION();
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!aux_vector)
    {
        return JMTX_RESULT_NULL_PARAM;
    }

    //  Length of x and y
    const uint32_t n = mtx->base.cols;
    double err = 0;
    uint32_t n_iterations = 0;

    _Complex double mag_y = 0;
    const double convergence_dif = args->in_convergence_criterion;
    double* const p_error_evolution = args->opt_error_evolution;
    const uint32_t n_max_iter = args->in_max_iterations;
#pragma omp parallel default(none) shared(n, err, n_iterations, x, y, mtx, aux_vector, convergence_dif, p_error_evolution, n_max_iter, mag_y)
    {
#pragma omp for reduction(+:mag_y) schedule(static)
        for (uint32_t i = 0; i < n; ++i)
        {
            _Complex double d = jmtxz_matrix_crs_get_entry(mtx, i, i);
            aux_vector[i] = 1.0f/d;
            mag_y += conjf(y[i]) * y[i];
        }

#pragma omp single
        {
            mag_y = sqrt(mag_y);
            err = 0;
        }

        do
        {
//#pragma omp barrier

#pragma omp for schedule(static)
            for (uint32_t i = 0; i < n; ++i)
            {
                x[i] = x[i] + (y[i] - jmtxz_matrix_crs_vector_multiply_row(mtx, x, i)) * aux_vector[i];
            }

#pragma omp master
            {
                err = 0;
            }
#pragma omp barrier

#pragma omp for reduction(+:err) schedule(static)
            for (uint32_t i = 0; i < n; ++i)
            {
                _Complex double val;
                val = jmtxz_matrix_crs_vector_multiply_row(mtx, x, i) - y[i];
                err += val * val;
            }

#pragma omp master
            {
                err = sqrt(err) / mag_y;
                if (p_error_evolution)
                {
                    p_error_evolution[n_iterations] = err;
                }

                n_iterations += 1;
            }
#pragma omp barrier
        } while(err > convergence_dif && n_iterations < n_max_iter);
    }

    args->out_last_iteration = n_iterations;
    args->out_last_error = err;
    return n_iterations == n_max_iter || !isfinite(err) ? JMTX_RESULT_NOT_CONVERGED : JMTX_RESULT_SUCCESS;
}
