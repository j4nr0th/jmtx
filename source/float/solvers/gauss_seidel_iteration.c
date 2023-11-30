//
// Created by jan on 16.6.2022.
//

#include "../../../include/jmtx/float/solvers/gauss_seidel_iteration.h"
#include "../matrices/sparse_row_compressed_internal.h"

#include <math.h>



jmtx_result jmtx_gauss_seidel_crs(const jmtx_matrix_crs* mtx, const float* restrict y, float* restrict x, jmtx_solver_arguments* args)
{
    if (!mtx)
    {
//        REPORT_ERROR_MESSAGE("Matrix pointer was null");
//        LEAVE_FUNCTION();
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTX_TYPE_CRS)
    {
//        REPORT_ERROR_MESSAGE("Matrix was not compressed row sparse");
//        LEAVE_FUNCTION();
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (!y)
    {
//        REPORT_ERROR_MESSAGE("Vector y pointer was null");
//        LEAVE_FUNCTION();
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!x)
    {
//        REPORT_ERROR_MESSAGE("Vector x pointer was null");
//        LEAVE_FUNCTION();
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
        float d = jmtx_matrix_crs_get_entry(mtx, i, i);
        x[i] /= d;
    }


    float err;
    uint32_t n_iterations = 0;
    do
    {

        for (uint32_t i = 0; i < n; ++i)
        {
            float* row_ptr;
            uint32_t* index_ptr;
            uint32_t n_elements = jmtx_matrix_crs_get_row(mtx, i, &index_ptr, &row_ptr);
            float res = (float)0.0;
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
            const float new_x = (y[i] - res) / row_ptr[k];
            x[i] = new_x;
        }

        err = 0;
        for (uint32_t i = 0; i < n; ++i)
        {
            float val;
            val = jmtx_matrix_crs_vector_multiply_row(mtx, x, i);
            val -= y[i];
            err += val * val;
        }
        err = sqrtf(err) / (float)n;
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

jmtx_result jmtx_gauss_seidel_crs_parallel(
        const jmtx_matrix_crs* const mtx, const float* restrict y, float* restrict x, float* restrict aux_vector,
        jmtx_solver_arguments* args)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTX_TYPE_CRS)
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
    float err = 0;
    uint32_t n_iterations = 0;

    float mag_y = 0;
    const float convergence_dif = args->in_convergence_criterion;
    float* const p_error_evolution = args->opt_error_evolution;
    const uint32_t n_max_iter = args->in_max_iterations;
#pragma omp parallel default(none) shared(n, err, n_iterations, x, y, mtx, aux_vector, convergence_dif, p_error_evolution, n_max_iter, mag_y)
    {
#pragma omp for reduction(+:mag_y) schedule(static)
        for (uint32_t i = 0; i < n; ++i)
        {
            float d = jmtx_matrix_crs_get_entry(mtx, i, i);
            aux_vector[i] = 1.0f/d;
            mag_y += y[i] * y[i];
        }

#pragma omp single
        {
            mag_y = sqrtf(mag_y);
            err = 0;
        }

        do
        {
//#pragma omp barrier

#pragma omp for schedule(static)
            for (uint32_t i = 0; i < n; ++i)
            {
                x[i] = x[i] + (y[i] - jmtx_matrix_crs_vector_multiply_row(mtx, x, i)) * aux_vector[i];
            }

#pragma omp master
            {
                err = 0;
            }
#pragma omp barrier

#pragma omp for reduction(+:err) schedule(static)
            for (uint32_t i = 0; i < n; ++i)
            {
                float val;
                val = jmtx_matrix_crs_vector_multiply_row(mtx, x, i) - y[i];
                err += val * val;
            }

#pragma omp master
            {
                err = sqrtf(err) / mag_y;
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
