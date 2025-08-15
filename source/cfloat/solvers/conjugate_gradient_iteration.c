// Automatically generated from source/float/solvers/conjugate_gradient_iteration.c on Sun Dec 17 20:57:47 2023
//
// Created by jan on 30.10.2023.
//

#ifndef JMTXC_SOLVER_BASE_H
#include "../../../include/jmtx/solver_base.h"
#endif
#include "../../../include/jmtx/cfloat/solvers/cholesky_solving.h"
#include "../../../include/jmtx/cfloat/solvers/conjugate_gradient_iteration.h"
#include "../matrices/sparse_diagonal_compressed_internal.h"
#include "../matrices/sparse_row_compressed_internal.h"
#include <complex.h>
#include <math.h>

jmtx_result jmtxc_solve_iterative_conjugate_gradient_crs(const jmtxc_matrix_crs *mtx, const _Complex float *restrict y,
                                                         _Complex float *restrict x, _Complex float *restrict aux_vec1,
                                                         _Complex float *restrict aux_vec2,
                                                         _Complex float *restrict aux_vec3, jmtx_solver_arguments *args)
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
    if (!aux_vec1)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!aux_vec2)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!aux_vec3)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!args)
    {
        return JMTX_RESULT_NULL_PARAM;
    }

    const uint32_t n = mtx->base.rows;
    _Complex float *const r = aux_vec1;
    _Complex float *const p = aux_vec2;
    _Complex float *const Ap = aux_vec3;
    uint32_t n_iterations = 0;
    _Complex float alpha, beta;

    float rk_dp = 0;
    float err = 0;
    float mag_y = 0;
    _Complex float pAp_dp = 0;
    float new_rk_dp = 0;

    {
        for (unsigned i = 0; i < n; ++i)
        {
            r[i] = y[i] - jmtxc_matrix_crs_vector_multiply_row(mtx, x, i);
            p[i] = r[i];
        }

        rk_dp = 0;
        mag_y = 0;
        for (unsigned i = 0; i < n; ++i)
        {
            rk_dp += conjf(r[i]) * r[i];
            mag_y += conjf(y[i]) * y[i];
        }

        mag_y = sqrtf(mag_y);
        err = sqrtf(rk_dp) / mag_y;
        if (err < args->in_convergence_criterion)
        {
            args->out_last_error = err;
            args->out_last_iteration = 0;
            return JMTX_RESULT_SUCCESS;
        }

        for (;;)
        {
            //  Compute Ap
            pAp_dp = 0;
            for (uint_fast32_t i = 0; i < n; ++i)
            {
                Ap[i] = jmtxc_matrix_crs_vector_multiply_row(mtx, p, i);
                pAp_dp += conjf(p[i]) * Ap[i];
            }

            //  Once alpha goes too low (which happens when p vectors become more and more co-linear), solution won't
            //  progress go forward anymore

            alpha = rk_dp / pAp_dp;
            err = 0;
            new_rk_dp = 0;

            //  Update guess of x
            for (unsigned i = 0; i < n; ++i)
            {
                x[i] += alpha * p[i];
            }

            //  Update guess of r
            new_rk_dp = 0;
            for (unsigned i = 0; i < n; ++i)
            {
                r[i] = r[i] - alpha * Ap[i];
                new_rk_dp += conjf(r[i]) * r[i];
            }

            err = sqrtf(new_rk_dp) / mag_y;
            if (args->opt_error_evolution)
            {
                args->opt_error_evolution[n_iterations] = err;
            }

            if (err < args->in_convergence_criterion)
            {
                //  Error is low enough
                break;
            }

            if (n_iterations == args->in_max_iterations)
            {
                break;
            }

            beta = new_rk_dp / rk_dp;
            rk_dp = new_rk_dp;

            for (unsigned i = 0; i < n; ++i)
            {
                p[i] = r[i] + beta * p[i];
            }

            n_iterations += 1;
        }
    }

    args->out_last_error = err;
    args->out_last_iteration = n_iterations;
    if (!isfinite(err) || err >= args->in_convergence_criterion)
    {
        return JMTX_RESULT_NOT_CONVERGED;
    }
    return err < args->in_convergence_criterion ? JMTX_RESULT_SUCCESS : JMTX_RESULT_NOT_CONVERGED;
}

jmtx_result jmtxc_solve_iterative_conjugate_gradient_crs_parallel(
    const jmtxc_matrix_crs *mtx, const _Complex float *restrict y, _Complex float *restrict x,
    _Complex float *restrict aux_vec1, _Complex float *restrict aux_vec2, _Complex float *restrict aux_vec3,
    jmtx_solver_arguments *args)
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
    if (!aux_vec1)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!aux_vec2)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!aux_vec3)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!args)
    {
        return JMTX_RESULT_NULL_PARAM;
    }

    const uint32_t n = mtx->base.rows;
    _Complex float *const r = aux_vec1;
    _Complex float *const p = aux_vec2;
    _Complex float *const Ap = aux_vec3;
    _Complex float alpha = 0, beta = 0;

    float rk_dp = 0;
    _Complex float pAp_dp = 0;
    float err = 0;
    float new_rk_dp = 0;
    float mag_y = 0;
    float *const err_evolution = args->opt_error_evolution;
    const float tolerance = args->in_convergence_criterion;
    const uint32_t max_iterations = args->in_max_iterations;
    uint32_t n_iterations = 0;
#pragma omp parallel default(none) shared(n, rk_dp, r, mtx, p, Ap, alpha, err, x, y, n_iterations, err_evolution,      \
                                              tolerance, beta, max_iterations, pAp_dp, new_rk_dp, mag_y)
    {
#pragma omp for
        for (unsigned i = 0; i < n; ++i)
        {
            r[i] = y[i] - jmtxc_matrix_crs_vector_multiply_row(mtx, x, i);
            p[i] = r[i];
        }

#pragma omp for reduction(+ : rk_dp, mag_y)
        for (unsigned i = 0; i < n; ++i)
        {
            rk_dp += conjf(r[i]) * r[i];
            mag_y += conjf(y[i]) * y[i];
        }

#pragma omp single
        {
            mag_y = sqrtf(mag_y);
            err = sqrtf(rk_dp) / mag_y;
        }

        for (;;)
        {
#pragma omp barrier
            //  Compute Ap
#pragma omp for reduction(+ : pAp_dp)
            for (unsigned i = 0; i < n; ++i)
            {
                Ap[i] = jmtxc_matrix_crs_vector_multiply_row(mtx, p, i);
                pAp_dp += conjf(p[i]) * Ap[i];
            }

            //  Once alpha goes too low (which happens when p vectors become more and more co-linear), solution won't
            //  progress go forward anymore

#pragma omp single
            {
                alpha = rk_dp / pAp_dp;
                err = 0;
                new_rk_dp = 0;
                pAp_dp = 0;
            }

            //  Update guess of x
#pragma omp for
            for (unsigned i = 0; i < n; ++i)
            {
                x[i] += alpha * p[i];
            }

            //  Update guess of r

#pragma omp for reduction(+ : new_rk_dp)
            for (unsigned i = 0; i < n; ++i)
            {
                r[i] = r[i] - alpha * Ap[i];
                new_rk_dp += conjf(r[i]) * r[i];
            }

#pragma omp single
            {
                err = sqrtf(new_rk_dp) / mag_y;
                if (err_evolution)
                {
                    err_evolution[n_iterations] = err;
                }
            }

            if (n_iterations == max_iterations)
            {
                break;
            }
            if (err < tolerance)
            {
                //  Error is low enough
                break;
            }

#pragma omp single
            {
                beta = new_rk_dp / rk_dp;
                rk_dp = new_rk_dp;
            }

#pragma omp for
            for (unsigned i = 0; i < n; ++i)
            {
                p[i] = r[i] + beta * p[i];
            }

#pragma omp single
            {
                n_iterations += 1;
            }
#pragma omp barrier
        }
    }

    args->out_last_error = err;
    args->out_last_iteration = n_iterations;
    if (!isfinite(err) || err >= args->in_convergence_criterion)
    {
        return JMTX_RESULT_NOT_CONVERGED;
    }
    return JMTX_RESULT_SUCCESS;
}

jmtx_result jmtxc_incomplete_cholesky_preconditioned_solve_iterative_conjugate_gradient_crs(
    const jmtxc_matrix_crs *mtx, const jmtxc_matrix_crs *cho, const jmtxc_matrix_crs *cho_t,
    const _Complex float *restrict y, _Complex float *restrict x, _Complex float *restrict aux_vec1,
    _Complex float *restrict aux_vec2, _Complex float *restrict aux_vec3, _Complex float *restrict aux_vec4,
    jmtx_solver_arguments *args)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!cho)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!cho_t)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.rows != mtx->base.cols)
    {
        //  I am only doing square matrices!!!
        return JMTX_RESULT_BAD_MATRIX;
    }
    if (mtx->base.rows != cho->base.rows)
    {
        //  I am only doing square matrices!!!
        return JMTX_RESULT_BAD_MATRIX;
    }
    if (cho->base.rows != cho->base.cols)
    {
        //  I am only doing square matrices!!!
        return JMTX_RESULT_BAD_MATRIX;
    }
    if (mtx->base.rows != cho_t->base.rows)
    {
        //  I am only doing square matrices!!!
        return JMTX_RESULT_BAD_MATRIX;
    }
    if (cho_t->base.rows != cho_t->base.cols)
    {
        //  I am only doing square matrices!!!
        return JMTX_RESULT_BAD_MATRIX;
    }
    if (mtx->base.type != JMTXC_TYPE_CRS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (cho->base.type != JMTXC_TYPE_CRS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (cho_t->base.type != JMTXC_TYPE_CRS)
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
    if (!aux_vec1)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!aux_vec2)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!aux_vec3)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!aux_vec4)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!args)
    {
        return JMTX_RESULT_NULL_PARAM;
    }

    const uint32_t n = mtx->base.rows;
    _Complex float *const r = aux_vec1;
    _Complex float *const p = aux_vec2;
    _Complex float *const Ap = aux_vec3;
    _Complex float *const z = aux_vec4;
    uint32_t n_iterations = 0;
    _Complex float alpha, beta;

    float rk_dp = 0;
    _Complex float rkzk_dp = 0;
    float mag_y = 0;
    _Complex float pAp_dp = 0;
    _Complex float new_rkzk_dp = 0;
    float err = 0;

    {
        //  Compute initial residual
        for (unsigned i = 0; i < n; ++i)
        {
            r[i] = y[i] - jmtxc_matrix_crs_vector_multiply_row(mtx, x, i);
        }
        //  Compute initial z
        jmtxc_solve_direct_cholesky_crs(cho, cho_t, r, z);

        rk_dp = 0;
        mag_y = 0;
        for (unsigned i = 0; i < n; ++i)
        {
            p[i] = z[i];
            rk_dp += conjf(r[i]) * r[i];
            rkzk_dp += conjf(r[i]) * z[i];
            mag_y += conjf(y[i]) * y[i];
        }
        mag_y = sqrtf(mag_y);

        for (;;)
        {

            //  Compute Ap
            pAp_dp = 0;
            for (unsigned i = 0; i < n; ++i)
            {
                Ap[i] = jmtxc_matrix_crs_vector_multiply_row(mtx, p, i);
                pAp_dp += p[i] * Ap[i];
            }

            //  Once alpha goes too low (which happens when p vectors become more and more co-linear), solution won't
            //  progress go forward anymore

            alpha = rkzk_dp / pAp_dp;
            rk_dp = 0;

            //  Update guess of x
            for (unsigned i = 0; i < n; ++i)
            {
                x[i] += alpha * p[i];
            }

            //  Update guess of r
            //  Update it implicitly
            rk_dp = 0;
            for (unsigned i = 0; i < n; ++i)
            {
                r[i] = r[i] - alpha * Ap[i];
                rk_dp += conjf(r[i]) * r[i];
            }

            err = sqrtf(rk_dp) / mag_y;
            if (args->opt_error_evolution)
            {
                args->opt_error_evolution[n_iterations] = err;
            }

            if (err < args->in_convergence_criterion)
            {
                //  Error is low enough to stop
                break;
            }
            if (n_iterations == args->in_max_iterations)
            {
                break;
            }

            jmtxc_solve_direct_cholesky_crs(cho, cho_t, r, z);

            new_rkzk_dp = 0;
            for (uint32_t i = 0; i < n; ++i)
            {
                new_rkzk_dp += conjf(r[i]) * z[i];
            }

            beta = new_rkzk_dp / rkzk_dp;
            rkzk_dp = new_rkzk_dp;

            for (unsigned i = 0; i < n; ++i)
            {
                p[i] = z[i] + beta * p[i];
            }

            n_iterations += 1;
        }
    }

    args->out_last_error = sqrtf(rk_dp / mag_y);
    args->out_last_iteration = n_iterations;
    if (!isfinite(err) || err >= args->in_convergence_criterion)
    {
        return JMTX_RESULT_NOT_CONVERGED;
    }
    return JMTX_RESULT_SUCCESS;
}

jmtx_result jmtxc_solve_iterative_conjugate_gradient_cds(const jmtxc_matrix_cds *mtx, const _Complex float *restrict y,
                                                         _Complex float *restrict x, _Complex float *restrict aux_vec1,
                                                         _Complex float *restrict aux_vec2,
                                                         _Complex float *restrict aux_vec3, jmtx_solver_arguments *args)
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
    if (mtx->base.type != JMTXC_TYPE_CDS)
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
    if (!aux_vec1)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!aux_vec2)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!aux_vec3)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!args)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (args->in_convergence_criterion < 0.0f || !isfinite(args->in_convergence_criterion))
    {
        return JMTX_RESULT_BAD_PARAM;
    }
    const uint_fast32_t n = mtx->base.rows;
    for (uint_fast32_t i = 0; i < n; ++i)
    {
        if (!isfinite(crealf(x[i])) || !isfinite(cimagf(x[i])) || !isfinite(crealf(y[i])) || !isfinite(cimagf(y[i])))
        {
            return JMTX_RESULT_BAD_PARAM;
        }
    }

    _Complex float *const r = aux_vec1;
    _Complex float *const p = aux_vec2;
    _Complex float *const Ap = aux_vec3;
    float rkrk, new_rkrk;
    float mag_y = 0;
    _Complex float pkApk = 0;
    uint_fast32_t n_iteration = 0;
    float error;
    _Complex float alfa, beta;

    for (uint_fast32_t i = 0; i < n; ++i)
    {
        mag_y += conjf(y[i]) * y[i];
    }
    mag_y = sqrtf(mag_y);
    for (;;)
    {
        //  Explicitly compute the residual
        jmtxc_matrix_cds_vector_multiply(mtx, x, r);
        rkrk = 0;
        for (uint_fast32_t i = 0; i < n; ++i)
        {
            const _Complex float res = y[i] - r[i];
            r[i] = res;
            p[i] = res;
            rkrk += conjf(res) * res;
        }
        error = sqrtf(rkrk) / mag_y;
        if (n_iteration == args->in_max_iterations)
        {
            break;
        }
        if (error < args->in_convergence_criterion)
        {
            break;
        }
        //  Solve using implicitly computed residual
        for (;;)
        {
            jmtxc_matrix_cds_vector_multiply(mtx, p, Ap);
            pkApk = 0;
            for (uint_fast32_t i = 0; i < n; ++i)
            {
                pkApk += conjf(p[i]) * Ap[i];
            }
            alfa = rkrk / pkApk;
            for (uint_fast32_t i = 0; i < n; ++i)
            {
                x[i] = x[i] + alfa * p[i];
            }
            for (uint_fast32_t i = 0; i < n; ++i)
            {
                r[i] = r[i] - alfa * Ap[i];
            }
            new_rkrk = 0;
            for (uint_fast32_t i = 0; i < n; ++i)
            {
                new_rkrk += conjf(r[i]) * r[i];
            }

            beta = new_rkrk / rkrk;
            rkrk = new_rkrk;

            error = sqrtf(rkrk) / mag_y;
            if (args->opt_error_evolution)
            {
                args->opt_error_evolution[n_iteration] = error;
            }
            if (n_iteration == args->in_max_iterations)
            {
                break;
            }
            n_iteration += 1;
            if (error < args->in_convergence_criterion)
            {
                break;
            }

            for (uint_fast32_t i = 0; i < n; ++i)
            {
                p[i] = r[i] + beta * p[i];
            }
        }
    }

    args->out_last_error = error;
    args->out_last_iteration = n_iteration;

    if (!isfinite(error) || error > args->in_convergence_criterion)
    {
        return JMTX_RESULT_NOT_CONVERGED;
    }

    return JMTX_RESULT_SUCCESS;
}
