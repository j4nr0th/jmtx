//
// Created by jan on 6.11.2023.
//

#include <assert.h>
#include <math.h>
#include "../../include/jmtx/solvers/lu_solving.h"
#include "../matrices/sparse_row_compressed_internal.h"
#include "../matrices/band_row_major_internal.h"
#include "../../include/jmtx/solvers/incomplete_lu_decomposition.h"
#include "../../include/jmtx/matrices/sparse_conversion.h"

void jmtx_lu_solve_crs(const jmtx_matrix_crs* l, const jmtx_matrix_crs* u, const float* restrict y, float* restrict x)
{
    const uint32_t n = l->base.cols;
    x[0] = y[0];
    //  First is the forward substitution for L v = y
    for (uint32_t i = 1; i < n; ++i)
    {
        uint32_t* indices;
        float* values;
        uint32_t count = jmtx_matrix_crs_get_row(l, i, &indices, &values);
        assert(indices[count - 1] == (uint32_t)i);
        assert(values[count - 1] == 1.0f);

        float v = 0;
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
        float* values;
        uint32_t count = jmtx_matrix_crs_get_row(u, i, &indices, &values);
        assert(indices[0] == (uint32_t)i);

        float v = 0;
        for (uint32_t j = 1; j < count; ++j)
        {
            assert(indices[j] > (uint32_t)i);
            v += values[j] * x[indices[j]];
        }
        x[i] = (x[i] - v) / values[0];
    }
}

void jmtx_lu_solve_inplace_crs(const jmtx_matrix_crs* l, const jmtx_matrix_crs* u, float* restrict x)
{
    const uint32_t n = l->base.cols;
    //  First is the forward substitution for L v = y
    for (uint32_t i = 1; i < n; ++i)
    {
        uint32_t* indices;
        float* values;
        uint32_t count = jmtx_matrix_crs_get_row(l, i, &indices, &values);
        assert(indices[count - 1] == (uint32_t)i);

        float v = 0;
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
        float* values;
        uint32_t count = jmtx_matrix_crs_get_row(u, i, &indices, &values);
        assert(indices[0] == (uint32_t)i);

        float v = 0;
        for (uint32_t j = 1; j < count; ++j)
        {
            assert(indices[j] > (uint32_t)i);
            v += values[j] * x[indices[j]];
        }
        x[i] = (x[i] - v) / values[0];
    }
}

static inline void compute_residual(const uint32_t n, const jmtx_matrix_crs* mtx, const float* restrict x,
                                    const float* restrict y, float* restrict r)
{
    for (uint32_t i = 0; i < n; ++i)
    {
        r[i] = y[i] - jmtx_matrix_crs_vector_multiply_row(mtx, x, i);
    }
}

jmtx_result jmtx_incomplete_lu_decomposition_solve_crs(
        const jmtx_matrix_crs* mtx, const float* y, float* x, float* aux_vec, jmtx_solver_arguments* args,
        const jmtx_allocator_callbacks* allocator_callbacks)
{
#ifndef JMTX_NO_VERIFY_PARAMS
    if (!mtx)
    {
//        REPORT_ERROR_MESSAGE("Matrix pointer was null");
//        LEAVE_FUNCTION();
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.rows != mtx->base.cols)
    {
        //  I am only doing square matrices!!!
        return JMTX_RESULT_BAD_MATRIX;
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
#endif

    if (!allocator_callbacks)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }
    else if (!allocator_callbacks->alloc || !allocator_callbacks->free)
    {
        return JMTX_RESULT_BAD_PARAM;
    }



    jmtx_matrix_crs* lower;
    jmtx_matrix_ccs* upper_ccs;
    jmtx_result res = jmtx_incomplete_lu_crs(
            mtx, &lower, &upper_ccs, allocator_callbacks);
    if (res != JMTX_RESULT_SUCCESS)
    {
        return res;
    }
    jmtx_matrix_crs* upper;
    res = jmtx_convert_ccs_to_crs(upper_ccs, &upper, NULL);
    if (res != JMTX_RESULT_SUCCESS)
    {
        jmtx_matrix_crs_destroy(upper);
        jmtx_matrix_crs_destroy(lower);
        return res;
    }
    jmtx_matrix_ccs_destroy(upper_ccs);
    upper_ccs = NULL;

    res = jmtx_incomplete_lu_decomposition_solve_precomputed_crs(mtx, lower, upper, y, x, aux_vec, NULL);

    jmtx_matrix_crs_destroy(upper);
    jmtx_matrix_crs_destroy(lower);


    return res;
}

jmtx_result jmtx_incomplete_lu_decomposition_solve_precomputed_crs(
        const jmtx_matrix_crs* mtx, const jmtx_matrix_crs* l, const jmtx_matrix_crs* u, const float* y, float* x,
        float* aux_vec, jmtx_solver_arguments* args)
{
#ifndef JMTX_NO_VERIFY_PARAMS
    if (!mtx)
    {
//        REPORT_ERROR_MESSAGE("Matrix pointer was null");
//        LEAVE_FUNCTION();
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.rows != mtx->base.cols)
    {
        //  I am only doing square matrices!!!
        return JMTX_RESULT_BAD_MATRIX;
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
#endif


    const uint32_t n = mtx->base.rows;
    float* const r = aux_vec;
    compute_residual(n, mtx, x, y, r);
    float residual2 = 0;
    float mag_y2 = 0;
    for (uint32_t i = 0; i < n; ++i)
    {
        residual2 += r[i] * r[i];
        mag_y2 = y[i] * y[i];
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
    float* x1 = aux_vec;
    float* x2 = x;
    const float stop_mag2 = mag_y2 * (convergence_dif * convergence_dif);
    uint32_t n_iter;
    for (n_iter = 0; n_iter < n_max_iter && stop_mag2 < residual2; ++n_iter)
    {
        {
            float* tmp = x2;
            x2 = x1;
            x1 = tmp;
        }


        jmtx_lu_solve_inplace_crs(l, u, r);

        for (uint32_t i = 0; i < n; ++i)
        {
            x[i] += r[i];
        }

        compute_residual(n, mtx, x, y, r);

        residual2 = 0;
        for (uint32_t i = 0; i < n; ++i)
        {
            residual2 += r[i] * r[i];
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

jmtx_result jmtx_incomplete_lu_decomposition_solve_precomputed_crs_parallel(
        const jmtx_matrix_crs* mtx, const jmtx_matrix_crs* l, const jmtx_matrix_crs* u, const float* y, float* x,
        float* aux_vec, jmtx_solver_arguments* args)
{
#ifndef JMTX_NO_VERIFY_PARAMS
    if (!mtx)
    {
//        REPORT_ERROR_MESSAGE("Matrix pointer was null");
//        LEAVE_FUNCTION();
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.rows != mtx->base.cols)
    {
        //  I am only doing square matrices!!!
        return JMTX_RESULT_BAD_MATRIX;
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
#endif


    const uint32_t n = mtx->base.rows;
    float* const r = aux_vec;
    compute_residual(n, mtx, x, y, r);
    float residual2 = 0;
    float mag_y2 = 0;
    for (uint32_t i = 0; i < n; ++i)
    {
        residual2 += r[i] * r[i];
        mag_y2 = y[i] * y[i];
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
    float* x1 = aux_vec;
    float* x2 = x;
    const float stop_mag2 = mag_y2 * (convergence_dif * convergence_dif);
    uint32_t n_iter = 0;


    for (n_iter = 0; n_iter < n_max_iter && stop_mag2 < residual2; ++n_iter)
    {
        {
            float* tmp = x2;
            x2 = x1;
            x1 = tmp;
        }


        jmtx_lu_solve_inplace_crs(l, u, r);

        for (uint32_t i = 0; i < n; ++i)
        {
            x[i] += r[i];
        }

        compute_residual(n, mtx, x, y, r);

        residual2 = 0;
        for (uint32_t i = 0; i < n; ++i)
        {
            residual2 += r[i] * r[i];
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

jmtx_result jmtx_incomplete_lu_decomposition_solve_crs_parallel(
        const jmtx_matrix_crs* mtx, const float* y, float* x, float* aux_vec, jmtx_solver_arguments* args,
        const jmtx_allocator_callbacks* allocator_callbacks)
{
#ifndef JMTX_NO_VERIFY_PARAMS
    if (!mtx)
    {
//        REPORT_ERROR_MESSAGE("Matrix pointer was null");
//        LEAVE_FUNCTION();
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.rows != mtx->base.cols)
    {
        //  I am only doing square matrices!!!
        return JMTX_RESULT_BAD_MATRIX;
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
#endif

    if (!allocator_callbacks)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }
    else if (!allocator_callbacks->alloc || !allocator_callbacks->free)
    {
        return JMTX_RESULT_BAD_PARAM;
    }



    jmtx_matrix_crs* lower;
    jmtx_matrix_ccs* upper_ccs;
    jmtx_result res = jmtx_incomplete_lu_crs(
            mtx, &lower, &upper_ccs, allocator_callbacks);
    if (res != JMTX_RESULT_SUCCESS)
    {
        return res;
    }
    jmtx_matrix_crs* upper;
    res = jmtx_convert_ccs_to_crs(upper_ccs, &upper, NULL);
    if (res != JMTX_RESULT_SUCCESS)
    {
        jmtx_matrix_crs_destroy(upper);
        jmtx_matrix_crs_destroy(lower);
        return res;
    }
    jmtx_matrix_ccs_destroy(upper_ccs);
    upper_ccs = NULL;

    res = jmtx_incomplete_lu_decomposition_solve_precomputed_crs_parallel(mtx, lower, upper, y, x, aux_vec, args);

    jmtx_matrix_crs_destroy(upper);
    jmtx_matrix_crs_destroy(lower);


    return res;
}

void jmtx_lu_solve_brm(const jmtx_matrix_brm* l, const jmtx_matrix_brm* u, const float* y, float* x)
{
    const uint_fast32_t n = l->base.cols;
    x[0] = y[0];
    //  First is the forward substitution for L v = y
    for (uint_fast32_t i = 1; i < n; ++i)
    {
        const uint_fast32_t off_j = jmtx_matrix_brm_first_pos_in_row(l, i);
        float* values = NULL;
        const uint_fast32_t len = jmtx_matrix_brm_get_row(l, i, &values);
        float v = 0;
        for (uint_fast32_t j = 0; j < len - 1; ++j)
        {
            v += values[j] * x[off_j + j];
        }
        x[i] = y[i] - v;
    }
    //  Then the backward substitution for U x = v
    for (int32_t i = (int32_t)n - 1; i >= 0; --i)
    {
        const uint_fast32_t off_j = jmtx_matrix_brm_first_pos_in_row(u, i);
        float* values;
        const uint_fast32_t len = jmtx_matrix_brm_get_row(u, i, &values);
        float v = 0;
        for (uint32_t j = 1; j < len; ++j)
        {
            v += values[j] * x[off_j + j];
        }
        x[i] = (x[i] - v) / values[0];
    }
}

void jmtx_lu_solve_inplace_brm(const jmtx_matrix_brm* l, const jmtx_matrix_brm* u, float* x)
{
    const uint_fast32_t n = l->base.cols;
//    x[0] = x[0];
    //  First is the forward substitution for L v = y
    for (uint_fast32_t i = 1; i < n; ++i)
    {
        const uint_fast32_t off_j = jmtx_matrix_brm_first_pos_in_row(l, i);
        float* values = NULL;
        const uint_fast32_t len = jmtx_matrix_brm_get_row(l, i, &values);
        float v = 0;
        for (uint_fast32_t j = 0; j < len - 1; ++j)
        {
            v += values[j] * x[off_j + j];
        }
        x[i] = x[i] - v;
    }
    //  Then the backward substitution for U x = v
    for (int32_t i = (int32_t)n - 1; i >= 0; --i)
    {
        const uint_fast32_t off_j = jmtx_matrix_brm_first_pos_in_row(u, i);
        float* values;
        const uint_fast32_t len = jmtx_matrix_brm_get_row(u, i, &values);
        float v = 0;
        for (uint32_t j = 1; j < len; ++j)
        {
            v += values[j] * x[off_j + j];
        }
        x[i] = (x[i] - v) / values[0];
    }
}
