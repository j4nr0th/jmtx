#include "../solver_base.h"
#include "../matrices/sparse_diagonal_compressed.h"
#include "../matrices/sparse_row_compressed.h"
#include "cholesky_solving.h"
#include "conjugate_gradient_iteration.h"
#include <math.h>

jmtx_result JMTX_NAME_TYPED(solve_iterative_conjugate_gradient_crs)(
    const JMTX_NAME_TYPED(matrix_crs) * mtx, const JMTX_SCALAR_T *restrict y, JMTX_SCALAR_T *restrict x,
    JMTX_SCALAR_T *restrict aux_vec1, JMTX_SCALAR_T *restrict aux_vec2, JMTX_SCALAR_T *restrict aux_vec3,
    JMTX_NAME_TYPED(solver_arguments) * args)
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
    if (mtx->base.type != JMTXD_TYPE_CRS)
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

    const JMTX_INDEX_T n = mtx->base.rows;
    JMTX_SCALAR_T *const r = aux_vec1;
    JMTX_SCALAR_T *const p = aux_vec2;
    JMTX_SCALAR_T *const Ap = aux_vec3;
    JMTX_INDEX_T n_iterations = 0;
    JMTX_SCALAR_T alpha;
    JMTX_REAL_T beta;

    JMTX_REAL_T rk_dp = 0;
    JMTX_REAL_T err = 0;
    JMTX_REAL_T mag_y = 0;
    JMTX_SCALAR_T pAp_dp = 0;
    JMTX_REAL_T new_rk_dp = 0;

    {
        for (unsigned i = 0; i < n; ++i)
        {
            r[i] = y[i] - JMTX_NAME_TYPED(matrix_crs_vector_multiply_row)(mtx, x, i);
            p[i] = r[i];
        }

        rk_dp = 0;
        mag_y = 0;
        for (unsigned i = 0; i < n; ++i)
        {
            rk_dp += JMTX_DOT(r[i], r[i]);
            mag_y += JMTX_DOT(y[i], y[i]);
        }

        mag_y = JMTX_REAL_ROOT(mag_y);
        err = JMTX_REAL_ROOT(rk_dp) / mag_y;
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
            for (JMTX_FAST_INT_T i = 0; i < n; ++i)
            {
                Ap[i] = JMTX_NAME_TYPED(matrix_crs_vector_multiply_row)(mtx, p, i);
                pAp_dp += JMTX_DOT(p[i], Ap[i]);
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
                new_rk_dp += JMTX_DOT(r[i], r[i]);
            }

            err = JMTX_REAL_ROOT(new_rk_dp) / mag_y;
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

jmtx_result JMTX_NAME_TYPED(solve_iterative_conjugate_gradient_crs_parallel)(
    const JMTX_NAME_TYPED(matrix_crs) * mtx, const JMTX_SCALAR_T *restrict y, JMTX_SCALAR_T *restrict x,
    JMTX_SCALAR_T *restrict aux_vec1, JMTX_SCALAR_T *restrict aux_vec2, JMTX_SCALAR_T *restrict aux_vec3,
    JMTX_NAME_TYPED(solver_arguments) * args)
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
    if (mtx->base.type != JMTXD_TYPE_CRS)
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

    const JMTX_INDEX_T n = mtx->base.rows;
    JMTX_SCALAR_T *const r = aux_vec1;
    JMTX_SCALAR_T *const p = aux_vec2;
    JMTX_SCALAR_T *const Ap = aux_vec3;
    JMTX_SCALAR_T alpha = 0, beta = 0;

    JMTX_REAL_T rk_dp = 0;
    JMTX_SCALAR_T pAp_dp = 0;
    JMTX_REAL_T err = 0;
    JMTX_REAL_T new_rk_dp = 0;
    JMTX_REAL_T mag_y = 0;
    JMTX_REAL_T *const err_evolution = args->opt_error_evolution;
    const JMTX_REAL_T tolerance = args->in_convergence_criterion;
    const JMTX_INDEX_T max_iterations = args->in_max_iterations;
    JMTX_INDEX_T n_iterations = 0;
#pragma omp parallel default(none) shared(n, rk_dp, r, mtx, p, Ap, alpha, err, x, y, n_iterations, err_evolution,      \
                                              tolerance, beta, max_iterations, pAp_dp, new_rk_dp, mag_y)
    {
#pragma omp for
        for (unsigned i = 0; i < n; ++i)
        {
            r[i] = y[i] - JMTX_NAME_TYPED(matrix_crs_vector_multiply_row)(mtx, x, i);
            p[i] = r[i];
        }

#pragma omp for reduction(+ : rk_dp, mag_y)
        for (unsigned i = 0; i < n; ++i)
        {
            rk_dp += JMTX_DOT(r[i], r[i]);
            mag_y += JMTX_DOT(y[i], y[i]);
        }

#pragma omp single
        {
            mag_y = JMTX_REAL_ROOT(mag_y);
            err = JMTX_REAL_ROOT(rk_dp) / mag_y;
        }

        for (;;)
        {
#pragma omp barrier
            //  Compute Ap
#pragma omp for reduction(+ : pAp_dp)
            for (unsigned i = 0; i < n; ++i)
            {
                Ap[i] = JMTX_NAME_TYPED(matrix_crs_vector_multiply_row)(mtx, p, i);
                pAp_dp += JMTX_DOT(p[i], Ap[i]);
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
                new_rk_dp += JMTX_DOT(r[i], r[i]);
            }

#pragma omp single
            {
                err = JMTX_REAL_ROOT(new_rk_dp) / mag_y;
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

jmtx_result JMTX_NAME_TYPED(incomplete_cholesky_preconditioned_solve_iterative_conjugate_gradient_crs)(
    const JMTX_NAME_TYPED(matrix_crs) * mtx, const JMTX_NAME_TYPED(matrix_crs) * cho,
    const JMTX_NAME_TYPED(matrix_crs) * cho_t, const JMTX_SCALAR_T *restrict y, JMTX_SCALAR_T *restrict x,
    JMTX_SCALAR_T *restrict aux_vec1, JMTX_SCALAR_T *restrict aux_vec2, JMTX_SCALAR_T *restrict aux_vec3,
    JMTX_SCALAR_T *restrict aux_vec4, JMTX_NAME_TYPED(solver_arguments) * args)
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
    if (mtx->base.type != JMTXD_TYPE_CRS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (cho->base.type != JMTXD_TYPE_CRS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (cho_t->base.type != JMTXD_TYPE_CRS)
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

    const JMTX_INDEX_T n = mtx->base.rows;
    JMTX_SCALAR_T *const r = aux_vec1;
    JMTX_SCALAR_T *const p = aux_vec2;
    JMTX_SCALAR_T *const Ap = aux_vec3;
    JMTX_SCALAR_T *const z = aux_vec4;
    JMTX_INDEX_T n_iterations = 0;
    JMTX_SCALAR_T alpha;
    JMTX_REAL_T beta;

    JMTX_REAL_T rk_dp = 0;
    JMTX_SCALAR_T rkzk_dp = 0;
    JMTX_REAL_T mag_y = 0;
    JMTX_SCALAR_T pAp_dp = 0;
    JMTX_SCALAR_T new_rkzk_dp = 0;
    JMTX_REAL_T err = 0;

    {
        //  Compute initial residual
        for (unsigned i = 0; i < n; ++i)
        {
            r[i] = y[i] - JMTX_NAME_TYPED(matrix_crs_vector_multiply_row)(mtx, x, i);
        }
        //  Compute initial z
        JMTX_NAME_TYPED(solve_direct_cholesky_crs)(cho, cho_t, r, z);

        rk_dp = 0;
        mag_y = 0;
        for (unsigned i = 0; i < n; ++i)
        {
            p[i] = z[i];
            rk_dp += JMTX_DOT(r[i], r[i]);
            rkzk_dp += JMTX_DOT(r[i], z[i]);
            mag_y += JMTX_DOT(y[i], y[i]);
        }
        mag_y = JMTX_REAL_ROOT(mag_y);

        for (;;)
        {
            //  Compute Ap
            pAp_dp = 0;
            for (unsigned i = 0; i < n; ++i)
            {
                Ap[i] = JMTX_NAME_TYPED(matrix_crs_vector_multiply_row)(mtx, p, i);
                pAp_dp += JMTX_DOT(p[i], Ap[i]);
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
                rk_dp += JMTX_DOT(r[i], r[i]);
            }

            err = JMTX_REAL_ROOT(rk_dp) / mag_y;
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

            JMTX_NAME_TYPED(solve_direct_cholesky_crs)(cho, cho_t, r, z);

            new_rkzk_dp = 0;
            for (JMTX_INDEX_T i = 0; i < n; ++i)
            {
                new_rkzk_dp += JMTX_DOT(r[i], z[i]);
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

    args->out_last_error = JMTX_REAL_ROOT(rk_dp / mag_y);
    args->out_last_iteration = n_iterations;
    if (!isfinite(err) || err >= args->in_convergence_criterion)
    {
        return JMTX_RESULT_NOT_CONVERGED;
    }
    return JMTX_RESULT_SUCCESS;
}

jmtx_result JMTX_NAME_TYPED(solve_iterative_conjugate_gradient_cds)(
    const JMTX_NAME_TYPED(matrix_cds) * mtx, const JMTX_SCALAR_T *restrict y, JMTX_SCALAR_T *restrict x,
    JMTX_SCALAR_T *restrict aux_vec1, JMTX_SCALAR_T *restrict aux_vec2, JMTX_SCALAR_T *restrict aux_vec3,
    JMTX_NAME_TYPED(solver_arguments) * args)
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
    if (mtx->base.type != JMTXD_TYPE_CDS)
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
    const JMTX_FAST_INT_T n = mtx->base.rows;

    JMTX_SCALAR_T *const r = aux_vec1;
    JMTX_SCALAR_T *const p = aux_vec2;
    JMTX_SCALAR_T *const Ap = aux_vec3;
    JMTX_SCALAR_T rkrk, new_rkrk;
    JMTX_REAL_T mag_y = 0;
    JMTX_SCALAR_T pkApk = 0;
    JMTX_FAST_INT_T n_iteration = 0;
    JMTX_SCALAR_T alfa;
    JMTX_REAL_T error, beta;

    for (JMTX_FAST_INT_T i = 0; i < n; ++i)
    {
        mag_y += JMTX_DOT(y[i], y[i]);
    }
    mag_y = JMTX_REAL_ROOT(mag_y);
    for (;;)
    {
        //  Explicitly compute the residual
        JMTX_NAME_TYPED(matrix_cds_vector_multiply)(mtx, x, r);
        rkrk = 0;
        for (JMTX_FAST_INT_T i = 0; i < n; ++i)
        {
            const JMTX_SCALAR_T res = y[i] - r[i];
            r[i] = res;
            p[i] = res;
            rkrk += JMTX_DOT(res, res);
        }
        error = JMTX_REAL_ROOT(rkrk) / mag_y;
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
            JMTX_NAME_TYPED(matrix_cds_vector_multiply)(mtx, p, Ap);
            pkApk = 0;
            for (JMTX_FAST_INT_T i = 0; i < n; ++i)
            {
                pkApk += JMTX_DOT(p[i], Ap[i]);
            }
            alfa = rkrk / pkApk;
            for (JMTX_FAST_INT_T i = 0; i < n; ++i)
            {
                x[i] = x[i] + alfa * p[i];
            }
            for (JMTX_FAST_INT_T i = 0; i < n; ++i)
            {
                r[i] = r[i] - alfa * Ap[i];
            }
            new_rkrk = 0;
            for (JMTX_FAST_INT_T i = 0; i < n; ++i)
            {
                new_rkrk += JMTX_DOT(r[i], r[i]);
            }

            beta = new_rkrk / rkrk;
            rkrk = new_rkrk;

            error = JMTX_REAL_ROOT(rkrk) / mag_y;
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

            for (JMTX_FAST_INT_T i = 0; i < n; ++i)
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
