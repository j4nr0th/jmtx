//
// Created by jan on 30.10.2023.
//

#ifndef JMTX_SOLVER_BASE_H
    #include "../../../include/jmtx/solver_base.h"
#endif
#include <math.h>
#include <stdio.h>
#include "../matrices/sparse_row_compressed_internal.h"
#include "../matrices/sparse_diagonal_compressed_internal.h"
#include "../../../include/jmtx/float/solvers/cholesky_solving.h"
#include "../../../include/jmtx/float/solvers/conjugate_gradient_iteration.h"

jmtx_result jmtx_conjugate_gradient_crs(
        const jmtx_matrix_crs* mtx, const float* y, float* x, const float stagnation,
        const uint32_t recalculation_interval, float* restrict aux_vec1, float* restrict aux_vec2,
        float* restrict aux_vec3, jmtx_solver_arguments* args)
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
    if (mtx->base.type != JMTX_TYPE_CRS)
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
    float* const r = aux_vec1;
    float* const p = aux_vec2;
    float* const Ap = aux_vec3;
    uint32_t n_iterations = 0, n_since_calc = 0;
    float alpha, beta;
    int update_residual = 0, stagnated = 0;

    float rk_dp = 0;
    float err = 0;
    float mag_y = 0;
    float prev_err = 0;
    float pAp_dp = 0;
    float new_rk_dp = 0;

    {
        for (unsigned i = 0; i < n; ++i)
        {
            r[i] = y[i] - jmtx_matrix_crs_vector_multiply_row(mtx, x, i);
            p[i] = r[i];
        }

        rk_dp = 0;
        mag_y = 0;
        for (unsigned i = 0; i < n; ++i)
        {
            rk_dp += r[i] * r[i];
            mag_y += y[i] * y[i];
        }

        mag_y = sqrtf(mag_y);
        prev_err = sqrtf(rk_dp);


        do
        {

            //  Compute Ap
            pAp_dp = 0;
            for (unsigned i = 0; i < n; ++i)
            {
                Ap[i] = jmtx_matrix_crs_vector_multiply_row(mtx, p, i);
                pAp_dp += p[i] * Ap[i];
            }

            //  Once alpha goes too low (which happens when p vectors become more and more co-linear), solution won't progress
            //  go forward anymore

            alpha = rk_dp / pAp_dp;
            if ((n_since_calc += 1) == recalculation_interval)
            {
                n_since_calc = 0;
                update_residual = 1;
            }
            else
            {
                update_residual = 0;
            }
            err = 0;
            new_rk_dp = 0;


            //  Update guess of x
            for (unsigned i = 0; i < n; ++i)
            {
                x[i] += alpha * p[i];
            }

            //  Update guess of r
            if (update_residual)
            {
                computing_residual_directly:
                //  Compute it directly
                new_rk_dp = 0;
                for (unsigned i = 0; i < n; ++i)
                {
                    r[i] = y[i] - jmtx_matrix_crs_vector_multiply_row(mtx, x, i);
                    new_rk_dp += r[i] * r[i];
                }
            }
            else
            {
                //  Update it implicitly
                new_rk_dp = 0;
                for (unsigned i = 0; i < n; ++i)
                {
                    r[i] = r[i] - alpha * Ap[i];
                    new_rk_dp += r[i] * r[i];
                }
            }

            err = sqrtf(new_rk_dp) / mag_y;
            if (args->opt_error_evolution)
            {
                args->opt_error_evolution[n_iterations] = err;
            }

            if (err < args->in_convergence_criterion)
            {
                //  Make sure the residual was computed directly
                if (update_residual == 0)
                {
                    err = 0;
                    update_residual = 1;
                    new_rk_dp = 0;
                    goto computing_residual_directly;
                }
                //  Error is low enough
                break;
            }

            beta = new_rk_dp / rk_dp;
            rk_dp = new_rk_dp;

            for (unsigned i = 0; i < n; ++i)
            {
                p[i] = r[i] + beta * p[i];
            }

            n_iterations += 1;
            if (fabsf(prev_err - err) < stagnation)
            {
                stagnated = 1;
            }
            else
            {
                stagnated = 0;
            }
            prev_err = err;
        } while(n_iterations < args->in_max_iterations && stagnated == 0);
    }


    args->out_last_error = err;
    args->out_last_iteration = n_iterations;
    if (stagnated) return JMTX_RESULT_STAGNATED;
    return err < args->in_convergence_criterion ? JMTX_RESULT_SUCCESS : JMTX_RESULT_NOT_CONVERGED;
}

jmtx_result jmtx_conjugate_gradient_crs_parallel(
        const jmtx_matrix_crs* mtx, const float* y, float* x, const float stagnation,
        const uint32_t recalculation_interval, float* restrict aux_vec1, float* restrict aux_vec2,
        float* restrict aux_vec3, jmtx_solver_arguments* args)
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
    if (mtx->base.type != JMTX_TYPE_CRS)
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
    float* const r = aux_vec1;
    float* const p = aux_vec2;
    float* const Ap = aux_vec3;
    uint32_t n_iterations = 0, n_since_calc = 0;
    float alpha = 0, beta = 0;
    int update_residual = 1, stagnated = 0;

    float rk_dp = 0;
    float pAp_dp = 0;
    float err = 0;
    float new_rk_dp = 0;
    float mag_y = 0;
    float prev_err = 0;
    float* const err_evolution = args->opt_error_evolution;
    const float tolerance = args->in_convergence_criterion;
    const uint32_t max_iterations = args->in_max_iterations;
#pragma omp parallel default(none) shared(n, rk_dp, r, mtx, p, Ap, alpha, update_residual, n_iterations, err, x, y,\
                                          err_evolution, tolerance, beta, max_iterations, pAp_dp, new_rk_dp, mag_y,\
                                          stagnation, prev_err, stagnated, recalculation_interval, n_since_calc)
    {
#pragma omp for
        for (unsigned i = 0; i < n; ++i)
        {
            r[i] = y[i] - jmtx_matrix_crs_vector_multiply_row(mtx, x, i);
            p[i] = r[i];
        }

#pragma omp for reduction(+:rk_dp, mag_y)
        for (unsigned i = 0; i < n; ++i)
        {
            rk_dp += r[i] * r[i];
            mag_y += y[i] * y[i];
        }

#pragma omp single
        {
            mag_y = sqrtf(mag_y);
            prev_err = sqrtf(rk_dp);
        }

        do
        {
#pragma omp barrier
            //  Compute Ap
#pragma omp for reduction(+:pAp_dp)
            for (unsigned i = 0; i < n; ++i)
            {
                Ap[i] = jmtx_matrix_crs_vector_multiply_row(mtx, p, i);
                pAp_dp += p[i] * Ap[i];
            }

            //  Once alpha goes too low (which happens when p vectors become more and more co-linear), solution won't progress
            //  go forward anymore

#pragma omp single
            {
                alpha = rk_dp / pAp_dp;
                if ((n_since_calc += 1) == recalculation_interval)
                {
                    n_since_calc = 0;
                    update_residual = 1;
                }
                else
                {
                    update_residual = 0;
                }
                err = 0;
                new_rk_dp = 0;
            };


            //  Update guess of x
#pragma omp for
            for (unsigned i = 0; i < n; ++i)
            {
                x[i] += alpha * p[i];
            }

            //  Update guess of r
            if (update_residual)
            {
computing_residual_directly:
                //  Compute it directly
#pragma omp for reduction(+:new_rk_dp)
                for (unsigned i = 0; i < n; ++i)
                {
                    r[i] = y[i] - jmtx_matrix_crs_vector_multiply_row(mtx, x, i);
                    new_rk_dp += r[i] * r[i];
                }
            }
            else
            {
                //  Update it implicitly
#pragma omp for reduction(+:new_rk_dp)
                for (unsigned i = 0; i < n; ++i)
                {
                    r[i] = r[i] - alpha * Ap[i];
                    new_rk_dp += r[i] * r[i];
                }
            }

#pragma omp single
            {
                err = sqrtf(new_rk_dp) / mag_y;
                if (err_evolution)
                {
                    err_evolution[n_iterations] = err;
                }
            };

            if (err < tolerance)
            {
                //  Make sure the residual was computed directly
                if (update_residual == 0)
                {
#pragma omp single
                    {
                        err = 0;
                        update_residual = 1;
                        new_rk_dp = 0;
                    }
                    goto computing_residual_directly;
                }
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
                if (fabsf(prev_err - err) < stagnation)
                {
                    stagnated = 1;
                }
                else
                {
                    stagnated = 0;
                }
                prev_err = err;
            }
#pragma omp barrier
        } while(n_iterations < max_iterations && stagnated == 0);
    }

    args->out_last_error = err;
    args->out_last_iteration = n_iterations;
    if (stagnated) return JMTX_RESULT_STAGNATED;
    return err < tolerance ? JMTX_RESULT_SUCCESS: JMTX_RESULT_NOT_CONVERGED;
}

jmtx_result jmtx_incomplete_cholesky_preconditioned_conjugate_gradient_crs(
        const jmtx_matrix_crs* mtx, const jmtx_matrix_crs* cho, const jmtx_matrix_crs* cho_t, const float* restrict y,
        float* restrict x, float stagnation, uint32_t recalculation_interval, float* restrict aux_vec1,
        float* restrict aux_vec2, float* restrict aux_vec3, float* restrict aux_vec4, jmtx_solver_arguments* args)
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
    if (mtx->base.type != JMTX_TYPE_CRS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (cho->base.type != JMTX_TYPE_CRS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (cho_t->base.type != JMTX_TYPE_CRS)
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
    float* const r = aux_vec1;
    float* const p = aux_vec2;
    float* const Ap = aux_vec3;
    float* const z = aux_vec4;
    uint32_t n_iterations = 0, n_since_calc = 0;
    float alpha, beta;
    int update_residual = 0, stagnated = 0;

    float rk_dp = 0;
    float rkzk_dp = 0;
    float mag_y = 0;
    float prev_err = 0;
    float pAp_dp = 0;
    float new_rkzk_dp = 0;
    const float stop_criterion = mag_y * (args->in_convergence_criterion * args->in_convergence_criterion);
    float stagnation_criterion;
    float err = 0;

    {
        //  Compute initial residual
        for (unsigned i = 0; i < n; ++i)
        {
            r[i] = y[i] - jmtx_matrix_crs_vector_multiply_row(mtx, x, i);
        }
        //  Compute initial z
        jmtx_cholesky_solve(cho, cho_t, r, z);

        rk_dp = 0;
        mag_y = 0;
        for (unsigned i = 0; i < n; ++i)
        {
            p[i] = z[i];
            rk_dp += r[i] * r[i];
            rkzk_dp += r[i] * z[i];
            mag_y += y[i] * y[i];
        }
        stagnation_criterion = stagnation * sqrtf(mag_y);

        do
        {

            //  Compute Ap
            pAp_dp = 0;
            for (unsigned i = 0; i < n; ++i)
            {
                Ap[i] = jmtx_matrix_crs_vector_multiply_row(mtx, p, i);
                pAp_dp += p[i] * Ap[i];
            }

            //  Once alpha goes too low (which happens when p vectors become more and more co-linear), solution won't progress
            //  go forward anymore

            alpha = rkzk_dp / pAp_dp;
            if ((n_since_calc += 1) == recalculation_interval)
            {
                n_since_calc = 0;
                update_residual = 1;
            }
            else
            {
                update_residual = 0;
            }
            rk_dp = 0;


            //  Update guess of x
            for (unsigned i = 0; i < n; ++i)
            {
                x[i] += alpha * p[i];
            }

            //  Update guess of r
            if (update_residual)
            {
                computing_residual_directly:
                //  Compute it directly
                rk_dp = 0;
                for (unsigned i = 0; i < n; ++i)
                {
                    r[i] = y[i] - jmtx_matrix_crs_vector_multiply_row(mtx, x, i);
                    rk_dp += r[i] * r[i];
                }
            }
            else
            {
                //  Update it implicitly
                rk_dp = 0;
                for (unsigned i = 0; i < n; ++i)
                {
                    r[i] = r[i] - alpha * Ap[i];
                    rk_dp += r[i] * r[i];
                }
            }

            err = sqrtf(rk_dp / mag_y);
            if (args->opt_error_evolution)
            {
                args->opt_error_evolution[n_iterations] = err;
            }

            if (rk_dp < stop_criterion)
            {
                //  Make sure the residual was computed directly
                if (update_residual == 0)
                {
                    update_residual = 1;
                    rk_dp = 0;
                    goto computing_residual_directly;
                }
                //  Error is low enough to stop
                break;
            }

            jmtx_cholesky_solve(cho, cho_t, r, z);

            new_rkzk_dp = 0;
            for (uint32_t i = 0; i < n; ++i)
            {
                new_rkzk_dp += r[i] * z[i];
            }

            beta = new_rkzk_dp / rkzk_dp;
            rkzk_dp = new_rkzk_dp;

            for (unsigned i = 0; i < n; ++i)
            {
                p[i] = z[i] + beta * p[i];
            }

            n_iterations += 1;
            if (fabsf(prev_err - err) < stagnation_criterion)
            {
                stagnated = 1;
            }
            else
            {
                stagnated = 0;
            }
            prev_err = rk_dp;
        } while(n_iterations < args->in_max_iterations && stagnated == 0);
    }


    args->out_last_error = sqrtf(rk_dp / mag_y);
    args->out_last_iteration = n_iterations;
    if (stagnated) return JMTX_RESULT_STAGNATED;
    return args->out_last_error < args->in_convergence_criterion ? JMTX_RESULT_SUCCESS : JMTX_RESULT_NOT_CONVERGED;
}

jmtx_result jmtx_conjugate_gradient_cds(const jmtx_matrix_cds* mtx, const float* restrict y, float* restrict x,
                                        float* restrict aux_vec1, float* restrict aux_vec2, float* restrict aux_vec3,
                                        jmtx_solver_arguments* args)
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
    if (mtx->base.type != JMTX_TYPE_CDS)
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
        if (!isfinite(x[i]) || !isfinite(y[i]))
        {
            return JMTX_RESULT_BAD_PARAM;
        }
    }

    float* const r = aux_vec1;
    float* const p = aux_vec2;
    float* const Ap = aux_vec3;
    float rkrk, new_rkrk;
    float mag_y2 = 0;
    float pkApk = 0;
    uint_fast32_t n_iteration = 0;
    float error, alfa, beta;

    for (uint_fast32_t i = 0; i < n; ++i)
    {
        mag_y2 += y[i] * y[i];
    }
    for (;;)
    {
        //  Explicitly compute the residual
        jmtx_matrix_cds_vector_multiply(mtx, x, r);
        rkrk = 0;
        for (uint_fast32_t i = 0; i < n; ++i)
        {
            const float res = y[i] - r[i];
            r[i] = res;
            p[i] = res;
            rkrk += res * res;
        }
        error = sqrtf(rkrk / mag_y2);
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
            jmtx_matrix_cds_vector_multiply(mtx, p, Ap);
            pkApk = 0;
            for (uint_fast32_t i = 0; i < n; ++i)
            {
                pkApk += p[i] * Ap[i];
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
                new_rkrk += r[i] * r[i];
            }

            beta = new_rkrk / rkrk;
            rkrk = new_rkrk;

            error = sqrtf(rkrk / mag_y2);
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
