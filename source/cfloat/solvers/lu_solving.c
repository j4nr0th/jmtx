// Automatically generated from source/float/solvers/lu_solving.c on Fri Dec  1 17:36:03 2023
//
// Created by jan on 6.11.2023.
//

#include <assert.h>
#include <math.h>
#include "../../../include/jmtx/cfloat/solvers/lu_solving.h"
#include "../matrices/sparse_row_compressed_internal.h"
#include "../matrices/band_row_major_internal.h"
#include "../../../include/jmtx/cfloat/solvers/incomplete_lu_decomposition.h"
#include "../../../include/jmtx/cfloat/matrices/sparse_conversion.h"
#include <complex.h>

void jmtxc_lu_solve_crs(const jmtxc_matrix_crs* l, const jmtxc_matrix_crs* u, const _Complex float* restrict y, _Complex float* restrict x)
{
    const uint32_t n = l->base.cols;
    x[0] = y[0];
    //  First is the forward substitution for L v = y
    for (uint32_t i = 1; i < n; ++i)
    {
        uint32_t* indices;
        _Complex float* values;
        uint32_t count = jmtxc_matrix_crs_get_row(l, i, &indices, &values);
        assert(indices[count - 1] == (uint32_t)i);
        assert(values[count - 1] == 1.0f);

        _Complex float v = 0;
        for (uint32_t j = 0; j < count - 1; ++j)
        {
            assert(indices[j] < i);
            v += values[j] * x[indices[j]];
        }
        x[i] = y[i] - v;
    }
    //  Then the backward substitution for U x = v
    for (int32_t i = (int32_t)n - 1; i >= 0; --i)
    {
        uint32_t* indices;
        _Complex float* values;
        uint32_t count = jmtxc_matrix_crs_get_row(u, i, &indices, &values);
        assert(indices[0] == (uint32_t)i);

        _Complex float v = 0;
        for (uint32_t j = 1; j < count; ++j)
        {
            assert(indices[j] > (uint32_t)i);
            v += values[j] * x[indices[j]];
        }
        x[i] = (x[i] - v) / values[0];
    }
}

void jmtxc_lu_solve_inplace_crs(const jmtxc_matrix_crs* l, const jmtxc_matrix_crs* u, _Complex float* restrict x)
{
    const uint32_t n = l->base.cols;
    //  First is the forward substitution for L v = y
    for (uint32_t i = 1; i < n; ++i)
    {
        uint32_t* indices;
        _Complex float* values;
        uint32_t count = jmtxc_matrix_crs_get_row(l, i, &indices, &values);
        assert(indices[count - 1] == (uint32_t)i);

        _Complex float v = 0;
        for (uint32_t j = 0; j < count - 1; ++j)
        {
            assert(indices[j] < i);
            v += values[j] * x[indices[j]];
        }
        x[i] = x[i] - v;
    }
    //  Then the backward substitution for U x = v
    for (int32_t i = (int32_t)n - 1; i >= 0; --i)
    {
        uint32_t* indices;
        _Complex float* values;
        uint32_t count = jmtxc_matrix_crs_get_row(u, i, &indices, &values);
        assert(indices[0] == (uint32_t)i);

        _Complex float v = 0;
        for (uint32_t j = 1; j < count; ++j)
        {
            assert(indices[j] > (uint32_t)i);
            v += values[j] * x[indices[j]];
        }
        x[i] = (x[i] - v) / values[0];
    }
}

static inline void compute_residual(const uint32_t n, const jmtxc_matrix_crs* mtx, const _Complex float* restrict x,
                                    const _Complex float* restrict y, _Complex float* restrict r)
{
    for (uint32_t i = 0; i < n; ++i)
    {
        r[i] = y[i] - jmtxc_matrix_crs_vector_multiply_row(mtx, x, i);
    }
}

jmtx_result jmtxc_incomplete_lu_decomposition_solve_crs(
        const jmtxc_matrix_crs* mtx, const _Complex float* y, _Complex float* x, _Complex float* aux_vec, jmtx_solver_arguments* args,
        const jmtx_allocator_callbacks* allocator_callbacks)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.rows != mtx->base.cols)
    {
        return JMTX_RESULT_BAD_MATRIX;
    }
    if (mtx->base.type != JMTXC_TYPE_CRS)
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

    if (!allocator_callbacks)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }
    else if (!allocator_callbacks->alloc || !allocator_callbacks->free)
    {
        return JMTX_RESULT_BAD_PARAM;
    }



    jmtxc_matrix_crs* lower;
    jmtxc_matrix_ccs* upper_ccs;
    jmtx_result res = jmtxc_incomplete_lu_crs(
            mtx, &lower, &upper_ccs, allocator_callbacks);
    if (res != JMTX_RESULT_SUCCESS)
    {
        return res;
    }
    jmtxc_matrix_crs* upper;
    res = jmtxc_convert_ccs_to_crs(upper_ccs, &upper, NULL);
    if (res != JMTX_RESULT_SUCCESS)
    {
        jmtxc_matrix_crs_destroy(upper);
        jmtxc_matrix_crs_destroy(lower);
        return res;
    }
    jmtxc_matrix_ccs_destroy(upper_ccs);
    upper_ccs = NULL;

    res = jmtxc_incomplete_lu_decomposition_solve_precomputed_crs(mtx, lower, upper, y, x, aux_vec, NULL);

    jmtxc_matrix_crs_destroy(upper);
    jmtxc_matrix_crs_destroy(lower);


    return res;
}

jmtx_result jmtxc_incomplete_lu_decomposition_solve_precomputed_crs(
        const jmtxc_matrix_crs* mtx, const jmtxc_matrix_crs* l, const jmtxc_matrix_crs* u, const _Complex float* y, _Complex float* x,
        _Complex float* aux_vec, jmtx_solver_arguments* args)
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
    if (mtx->base.type != JMTXC_TYPE_CRS)
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


    const uint32_t n = mtx->base.rows;
    _Complex float* const r = aux_vec;
    compute_residual(n, mtx, x, y, r);
    float residual2 = 0;
    float mag_y2 = 0;
    for (uint32_t i = 0; i < n; ++i)
    {
        residual2 += r[i] * conjf(r[i]);
        mag_y2 += y[i] * conjf(y[i]);
    }

    float err = sqrtf(residual2/mag_y2);
    const float convergence_dif = args->in_convergence_criterion;
    const uint32_t n_max_iter = args->in_max_iterations;
    if (err < convergence_dif)
    {
        //  Converged prior to any iteration
        args->out_last_error = err;
        if (args->opt_error_evolution)
        {
            *args->opt_error_evolution = err;
        }
        args->out_last_iteration = 0;
        return JMTX_RESULT_SUCCESS;
    }



    //  Now begin iterations with x^{k + 1} = x^{k} + U'^{-1} L'^{-1} A x^{k}
    _Complex float* x1 = aux_vec;
    _Complex float* x2 = x;
    const float stop_mag2 = mag_y2 * (convergence_dif * convergence_dif);
    uint32_t n_iter;
    for (n_iter = 0; n_iter < n_max_iter && stop_mag2 < residual2; ++n_iter)
    {
        {
            _Complex float* tmp = x2;
            x2 = x1;
            x1 = tmp;
        }


        jmtxc_lu_solve_inplace_crs(l, u, r);

        for (uint32_t i = 0; i < n; ++i)
        {
            x[i] += r[i];
        }

        compute_residual(n, mtx, x, y, r);

        residual2 = 0;
        for (uint32_t i = 0; i < n; ++i)
        {
            residual2 += r[i] * conjf(r[i]);
        }

        if (args->opt_error_evolution)
        {
            args->opt_error_evolution[n_iter] = sqrtf(residual2 / mag_y2);
        }
    }
    err = sqrtf(residual2 / mag_y2);

    args->out_last_error = err;
    args->out_last_iteration = n_iter;

    return err < convergence_dif ? JMTX_RESULT_SUCCESS : JMTX_RESULT_NOT_CONVERGED;
}

jmtx_result jmtxc_incomplete_lu_decomposition_solve_precomputed_crs_parallel(
        const jmtxc_matrix_crs* mtx, const jmtxc_matrix_crs* l, const jmtxc_matrix_crs* u, const _Complex float* y, _Complex float* x,
        _Complex float* aux_vec, jmtx_solver_arguments* args)
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
    if (mtx->base.type != JMTXC_TYPE_CRS)
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


    const uint32_t n = mtx->base.rows;
    _Complex float* const r = aux_vec;
    compute_residual(n, mtx, x, y, r);
    float residual2 = 0;
    float mag_y2 = 0;
    for (uint32_t i = 0; i < n; ++i)
    {
        residual2 += r[i] * conjf(r[i]);
        mag_y2 += y[i] * conjf(y[i]);
    }

    float err = sqrtf(residual2/mag_y2);
    const float convergence_dif = args->in_convergence_criterion;
    const uint32_t n_max_iter = args->in_max_iterations;
    if (err < convergence_dif)
    {
        //  Converged prior to any iteration
        args->out_last_error = err;
        if (args->opt_error_evolution)
        {
            *args->opt_error_evolution = err;
        }

        args->out_last_iteration = 0;
        return JMTX_RESULT_SUCCESS;
    }



    //  Now begin iterations with x^{k + 1} = x^{k} + U'^{-1} L'^{-1} A x^{k}
    _Complex float* x1 = aux_vec;
    _Complex float* x2 = x;
    const float stop_mag2 = mag_y2 * (convergence_dif * convergence_dif);
    uint32_t n_iter;


    for (n_iter = 0; n_iter < n_max_iter && stop_mag2 < residual2; ++n_iter)
    {
        {
            _Complex float* tmp = x2;
            x2 = x1;
            x1 = tmp;
        }


        jmtxc_lu_solve_inplace_crs(l, u, r);

        for (uint32_t i = 0; i < n; ++i)
        {
            x[i] += r[i];
        }

        compute_residual(n, mtx, x, y, r);

        residual2 = 0;
        for (uint32_t i = 0; i < n; ++i)
        {
            residual2 += r[i] * conjf(r[i]);
        }

        if (args->opt_error_evolution)
        {
            args->opt_error_evolution[n_iter] = sqrtf(residual2 / mag_y2);
        }
    }
    err = sqrtf(residual2 / mag_y2);

    args->out_last_error = err;
    args->out_last_iteration = n_iter;

    return err < convergence_dif ? JMTX_RESULT_SUCCESS : JMTX_RESULT_NOT_CONVERGED;
}

jmtx_result jmtxc_incomplete_lu_decomposition_solve_crs_parallel(
        const jmtxc_matrix_crs* mtx, const _Complex float* y, _Complex float* x, _Complex float* aux_vec, jmtx_solver_arguments* args,
        const jmtx_allocator_callbacks* allocator_callbacks)
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
    if (mtx->base.type != JMTXC_TYPE_CRS)
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

    if (!allocator_callbacks)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }
    else if (!allocator_callbacks->alloc || !allocator_callbacks->free)
    {
        return JMTX_RESULT_BAD_PARAM;
    }



    jmtxc_matrix_crs* lower;
    jmtxc_matrix_ccs* upper_ccs;
    jmtx_result res = jmtxc_incomplete_lu_crs(
            mtx, &lower, &upper_ccs, allocator_callbacks);
    if (res != JMTX_RESULT_SUCCESS)
    {
        return res;
    }
    jmtxc_matrix_crs* upper;
    res = jmtxc_convert_ccs_to_crs(upper_ccs, &upper, NULL);
    if (res != JMTX_RESULT_SUCCESS)
    {
        jmtxc_matrix_crs_destroy(upper);
        jmtxc_matrix_crs_destroy(lower);
        return res;
    }
    jmtxc_matrix_ccs_destroy(upper_ccs);
    upper_ccs = NULL;

    res = jmtxc_incomplete_lu_decomposition_solve_precomputed_crs_parallel(mtx, lower, upper, y, x, aux_vec, args);

    jmtxc_matrix_crs_destroy(upper);
    jmtxc_matrix_crs_destroy(lower);


    return res;
}

void jmtxc_lu_solve_brm(const jmtxc_matrix_brm* l, const jmtxc_matrix_brm* u, const _Complex float* y, _Complex float* x)
{
    const uint_fast32_t n = l->base.cols;
    x[0] = y[0];
    //  First is the forward substitution for L v = y
    for (uint_fast32_t i = 1; i < n; ++i)
    {
        const uint_fast32_t off_j = jmtxc_matrix_brm_first_pos_in_row(l, i);
        _Complex float* values = NULL;
        const uint_fast32_t len = jmtxc_matrix_brm_get_row(l, i, &values);
        _Complex float v = 0;
        for (uint_fast32_t j = 0; j < len - 1; ++j)
        {
            v += values[j] * x[off_j + j];
        }
        x[i] = y[i] - v;
    }
    //  Then the backward substitution for U x = v
    for (int32_t i = (int32_t)n - 1; i >= 0; --i)
    {
        const uint_fast32_t off_j = jmtxc_matrix_brm_first_pos_in_row(u, i);
        _Complex float* values;
        const uint_fast32_t len = jmtxc_matrix_brm_get_row(u, i, &values);
        _Complex float v = 0;
        for (uint32_t j = 1; j < len; ++j)
        {
            v += values[j] * x[off_j + j];
        }
        x[i] = (x[i] - v) / values[0];
    }
}

void jmtxc_lu_solve_inplace_brm(const jmtxc_matrix_brm* l, const jmtxc_matrix_brm* u, _Complex float* x)
{
    const uint_fast32_t n = l->base.cols;
//    x[0] = x[0];
    //  First is the forward substitution for L v = y
    for (uint_fast32_t i = 1; i < n; ++i)
    {
        const uint_fast32_t off_j = jmtxc_matrix_brm_first_pos_in_row(l, i);
        _Complex float* values = NULL;
        const uint_fast32_t len = jmtxc_matrix_brm_get_row(l, i, &values);
        _Complex float v = 0;
        for (uint_fast32_t j = 0; j < len - 1; ++j)
        {
            v += values[j] * x[off_j + j];
        }
        x[i] = x[i] - v;
    }
    //  Then the backward substitution for U x = v
    for (int32_t i = (int32_t)n - 1; i >= 0; --i)
    {
        const uint_fast32_t off_j = jmtxc_matrix_brm_first_pos_in_row(u, i);
        _Complex float* values;
        const uint_fast32_t len = jmtxc_matrix_brm_get_row(u, i, &values);
        _Complex float v = 0;
        for (uint32_t j = 1; j < len; ++j)
        {
            v += values[j] * x[off_j + j];
        }
        x[i] = (x[i] - v) / values[0];
    }
}

jmtx_result jmtxc_lu_solve_iterative_bmr(const jmtxc_matrix_brm* a, const jmtxc_matrix_brm* l, const jmtxc_matrix_brm* u,
                                        const _Complex float y[const restrict], _Complex float x[const restrict],
                                        _Complex float aux_vec[const restrict],
                                        jmtx_solver_arguments* args)
{
    const uint_fast32_t n = l->base.cols;
    float y_magnitude2 = 0;
    uint_fast32_t iteration_count = 0;
    float error;
    for (uint_fast32_t i = 0; i < n; ++i)
    {
        y_magnitude2 += y[i] * conjf(y[i]);
    }
    jmtxc_lu_solve_brm(l, u, y, x);

    for (;;)
    {
        for (uint_fast32_t i = 0; i < n; ++i)
        {
            aux_vec[i] = jmtxc_matrix_brm_vector_multiply_row(a, x, i);
        }

        for (uint_fast32_t i = 0; i < n; ++i)
        {
            aux_vec[i] = y[i] - aux_vec[i];
        }

        _Complex float residual_magnitude2 = 0;
        for (uint_fast32_t i = 0; i < n; ++i)
        {
            residual_magnitude2 += conjf(aux_vec[i]) * aux_vec[i];
        }
        error = sqrtf(residual_magnitude2 / y_magnitude2);
        if (args->opt_error_evolution )
        {
            args->opt_error_evolution[iteration_count] = error;
        }
        ++iteration_count;
        if (error < args->in_convergence_criterion || iteration_count > args->in_max_iterations)
        {
            break;
        }
        jmtxc_lu_solve_inplace_brm(l, u, aux_vec);
        for (uint_fast32_t i = 0; i < n; ++i)
        {
            x[i] += aux_vec[i];
        }
    }
    args->out_last_error = error;
    args->out_last_iteration = iteration_count;
    return iteration_count < args->in_max_iterations ? JMTX_RESULT_SUCCESS : JMTX_RESULT_NOT_CONVERGED;
}

jmtx_result jmtxc_lu_solve_iterative_bmr_parallel(const jmtxc_matrix_brm* a, const jmtxc_matrix_brm* l,
                                                 const jmtxc_matrix_brm* u,  const _Complex float y[const restrict],
                                                 _Complex float x[const restrict], _Complex float aux_vec[const restrict],
                                                 jmtx_solver_arguments* args)
{
    const uint_fast32_t n = l->base.cols;
    float y_magnitude2 = 0;
    uint_fast32_t iteration_count = 0;
    float error = 0;
    for (uint_fast32_t i = 0; i < n; ++i)
    {
        y_magnitude2 += y[i] * conjf(y[i]);
    }
    float residual_magnitude2 = 0;
    jmtxc_lu_solve_brm(l, u, y, x);
#pragma omp parallel default(none) shared(y, x, l, u, aux_vec, args, n, error, y_magnitude2, iteration_count,\
                                          residual_magnitude2, a)
    {
        for (;;)
        {
#pragma omp for
            for (uint_fast32_t i = 0; i < n; ++i)
            {
                aux_vec[i] = jmtxc_matrix_brm_vector_multiply_row(a, x, i);
            }

#pragma omp for
            for (uint_fast32_t i = 0; i < n; ++i)
            {
                aux_vec[i] = y[i] - aux_vec[i];
            }

#pragma omp for reduction(+:residual_magnitude2)
            for (uint_fast32_t i = 0; i < n; ++i)
            {
                residual_magnitude2 += conjf(aux_vec[i]) * aux_vec[i];
            }

#pragma omp single
            {
                error = sqrtf(residual_magnitude2 / y_magnitude2);
                if (args->opt_error_evolution)
                {
                    args->opt_error_evolution[iteration_count] = error;
                }
                ++iteration_count;
            }

            if (error < args->in_convergence_criterion || iteration_count > args->in_max_iterations)
            {
                break;
            }
#pragma omp single
            {
                jmtxc_lu_solve_inplace_brm(l, u, aux_vec);
            }
#pragma omp for
            for (uint_fast32_t i = 0; i < n; ++i)
            {
                x[i] += aux_vec[i];
            }
        }
    }
    args->out_last_error = error;
    args->out_last_iteration = iteration_count;
    return iteration_count < args->in_max_iterations ? JMTX_RESULT_SUCCESS : JMTX_RESULT_NOT_CONVERGED;
}
