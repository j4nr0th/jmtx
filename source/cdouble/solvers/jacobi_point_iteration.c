// Automatically generated from source/cfloat/solvers/jacobi_point_iteration.c on Fri Dec  1 18:48:06 2023
// Automatically generated from source/cdouble/solvers/jacobi_point_iteration.c on Fri Dec  1 17:36:03 2023
//
// Created by jan on 15.6.2022.
//

#include "../../../include/jmtx/cdouble/solvers/jacobi_point_iteration.h"
#include "../matrices/sparse_row_compressed_internal.h"
#include <math.h>
#include <stdio.h>
#include <assert.h>

#include <omp.h>
#include <complex.h>

jmtx_result jmtxz_jacobi_crs(
        const jmtxz_matrix_crs* mtx, const _Complex double* restrict y, _Complex double* restrict x, _Complex double* restrict aux_vec1, _Complex double* restrict aux_vec2,
        jmtxd_solver_arguments* args)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.rows != mtx->base.cols)
    {
        //  I am only doing square matrices!!!
        return JMTX_RESULT_BAD_MATRIX;
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
    if (!aux_vec1)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!aux_vec2)
    {
        return JMTX_RESULT_NULL_PARAM;
    }



    //  Length of x and y
    const uint32_t n = mtx->base.cols;
    _Complex double* const div_factor = aux_vec1;
    double y_mag = 0;
    //  Initial guess by assuming that mtx is a diagonal matrix
    for (uint32_t i = 0; i < n; ++i)
    {
        _Complex double d = jmtxz_matrix_crs_get_entry(mtx, i, i);
        if (d == 0.0f)
        {
            //  Diagonal entry is zero!
            //  Can't solve this one with Jacobi
            return JMTX_RESULT_BAD_MATRIX;
        }
        x[i] = y[i] / d;
        div_factor[i] = 1.0f / d;
        y_mag += y[i] * conjf(y[i]);
    }
    y_mag = sqrt(y_mag);

    //  Memory used to store result of the current iteration
    _Complex double* const auxiliary_x = aux_vec2;

    _Complex double* x0 = auxiliary_x;
    _Complex double* x1 = x;
    double err;
    uint32_t n_iterations = 0;
    do
    {
        err = 0.0f;
        {
            _Complex double* tmp = x1;
            x1 = x0;
            x0 = tmp;
        }

        //  For each entry, find the corresponding row in matrix A - D and compute the dot product between x and that row
        for (uint32_t i = 0; i < n; ++i)
        {
            _Complex double* row_ptr;
            uint32_t* index_ptr;
            uint32_t n_elements = jmtxz_matrix_crs_get_row(mtx, i, &index_ptr, &row_ptr);
            _Complex double res = 0;
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
            _Complex double val = jmtxz_matrix_crs_vector_multiply_row(mtx, x1, i);
            val -= y[i];
            err += conjf(val) * val;
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

jmtx_result jmtxz_jacobi_relaxed_crs(
        const jmtxz_matrix_crs* mtx, const _Complex double* restrict y, _Complex double* restrict x, double relaxation_factor, _Complex double* restrict aux_vec1,
        _Complex double* restrict aux_vec2, jmtxd_solver_arguments* args)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.rows != mtx->base.cols)
    {
        //  I am only doing square matrices!!!
        return JMTX_RESULT_BAD_MATRIX;
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
    if (relaxation_factor <= 0.0f)
    {
        //  Relaxation factor must be strictly larger than 0
        return JMTX_RESULT_BAD_PARAM;
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
    const uint32_t n = mtx->base.cols;
    _Complex double* div_factor = aux_vec1;
    double y_mag = 0;
    //  Initial guess by assuming that mtx is a diagonal matrix
    for (uint32_t i = 0; i < n; ++i)
    {
        _Complex double d = jmtxz_matrix_crs_get_entry(mtx, i, i);
        if (d == 0.0f)
        {
            //  Diagonal entry is zero!
            //  Can't solve this one with Jacobi
            return JMTX_RESULT_BAD_MATRIX;
        }
        x[i] = y[i] / d;
        div_factor[i] = relaxation_factor / d;
        y_mag += conjf(y[i]) * y[i];
    }
    y_mag = sqrt(y_mag);

    //  Memory used to store result of the current iteration
    _Complex double* const auxiliary_x = aux_vec2;

    _Complex double* x0 = auxiliary_x;
    _Complex double* x1 = x;
    double err;
    uint32_t n_iterations = 0;
    do
    {
        err = 0.0f;
        {
            _Complex double* tmp = x1;
            x1 = x0;
            x0 = tmp;
        }

        //  For each entry, find the corresponding row in matrix A - D and compute the dot product between x and that row
        for (uint32_t i = 0; i < n; ++i)
        {
            _Complex double* row_ptr;
            uint32_t* index_ptr;
            uint32_t n_elements = jmtxz_matrix_crs_get_row(mtx, i, &index_ptr, &row_ptr);
            _Complex double res = 0;
            for (uint32_t j = 0; j < n_elements; ++j)
            {
                res += row_ptr[j] * x0[index_ptr[j]];
            }
            x1[i] = x0[i] + (y[i] - res) * div_factor[i];
        }

        for (uint32_t i = 0; i < n; ++i)
        {
            _Complex double val = jmtxz_matrix_crs_vector_multiply_row(mtx, x1, i);
            val -= y[i];
            err += val * conjf(val);
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
//    LEAVE_FUNCTION();
    return n_iterations == args->in_max_iterations ? JMTX_RESULT_NOT_CONVERGED : JMTX_RESULT_SUCCESS;
}


jmtx_result jmtxz_jacobi_crs_parallel(
        const jmtxz_matrix_crs* mtx, const _Complex double* restrict y, _Complex double* restrict x, _Complex double* restrict aux_vector1, _Complex double* restrict aux_vector2,
        jmtxd_solver_arguments* args)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.rows != mtx->base.cols)
    {
        //  I am only doing square matrices!!!
        return JMTX_RESULT_BAD_MATRIX;
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
    if (!aux_vector1)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!aux_vector2)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!args)
    {
        return JMTX_RESULT_NULL_PARAM;
    }

    //  Length of x and y
    const uint32_t n = mtx->base.cols;
    double y_mag = 0;
    //  Initial guess by assuming that mtx is a diagonal matrix
    int zero_diag = 0;
#pragma omp parallel for shared(mtx, zero_diag, n, x, y, aux_vector1) reduction(+:y_mag) default(none)
    for (uint32_t i = 0; i < n; ++i)
    {
        _Complex double d;
        if ((d = jmtxz_matrix_crs_get_entry(mtx, i, i)) == 0)
        {
            //  Diagonal entry is zero!
            //  Can't solve this one with Jacobi
            zero_diag += 1;
        }
        aux_vector1[i] = 1.0f / d;
        const _Complex double mag = conjf(y[i]) * y[i];
#pragma omp atomic
        y_mag += mag;
    }
    if (zero_diag != 0)
    {
        //  Either diagonal was zero, or matrix could not find the diagonal element
        return JMTX_RESULT_BAD_MATRIX;
    }
    y_mag = sqrt(y_mag);
    //  Memory used to store result of the current iteration
    _Complex double* x0 = aux_vector2;
    _Complex double* x1 = x;
    double err = 0.0f;
    uint32_t n_iterations = 0;
    double* const p_error_evolution = args->opt_error_evolution;
    const double convergence_dif = args->in_convergence_criterion;
    const uint32_t max_iterations = args->in_max_iterations;
#pragma omp parallel default(none) shared(err, x0, x1, y, mtx, aux_vector1, y_mag, p_error_evolution, n_iterations,\
    convergence_dif, max_iterations, n)
    {
        do
        {
#pragma omp barrier
#pragma omp master
            {
                err = 0.0f;
                _Complex double* tmp = x1;
                x1 = x0;
                x0 = tmp;
            }
#pragma omp barrier

            //  For each entry, find the corresponding row in matrix A - D and compute the dot product between x and that row
#pragma omp for schedule(static)
            for (uint32_t i = 0; i < n; ++i)
            {
                x1[i] = x0[i] + (y[i] - jmtxz_matrix_crs_vector_multiply_row(mtx, x0, i)) * aux_vector1[i];
            }

#pragma omp for reduction(+:err) schedule(static)
            for (uint32_t i = 0; i < n; ++i)
            {
                const _Complex double val = jmtxz_matrix_crs_vector_multiply_row(mtx, x1, i) - y[i];
                err += conjf(val) * val;
            }


#pragma omp master
            {
                err = sqrt(err) / (_Complex double) y_mag;
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
//    LEAVE_FUNCTION();
    return n_iterations == max_iterations || !isfinite(err) ? JMTX_RESULT_NOT_CONVERGED : JMTX_RESULT_SUCCESS;
}
