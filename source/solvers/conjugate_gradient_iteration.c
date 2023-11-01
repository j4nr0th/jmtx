//
// Created by jan on 30.10.2023.
//

#include "conjugate_gradient_iteration.h"
#include "../matrices/sparse_row_compressed_internal.h"
#include <math.h>
#include <stdio.h>

jmtx_result jmtx_conjugate_gradient_crs(const jmtx_matrix_crs* mtx, const float* restrict y, float* restrict x,
                                        const float tolerance, const float stagnation,
                                        const uint32_t recalculation_interval, const uint32_t max_iterations,
                                        uint32_t* p_final_iteration, float* restrict err_evolution,
                                        float* final_error, float* restrict aux_vec1, float* restrict aux_vec2,
                                        float* restrict aux_vec3)
{
#ifndef JMTX_NO_VERIFY_PARAMS
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
        return JMTX_RESULT_NULL_PARAM;
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
#endif

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
            r[i] = y[i] - jmtx_matrix_crs_vector_multiply_row_raw(mtx, x, i);
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
                Ap[i] = jmtx_matrix_crs_vector_multiply_row_raw(mtx, p, i);
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
                    r[i] = y[i] - jmtx_matrix_crs_vector_multiply_row_raw(mtx, x, i);
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
            if (err_evolution)
            {
                err_evolution[n_iterations] = err;
            }

            if (err < tolerance)
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
        } while(n_iterations < max_iterations && stagnated == 0);
    }


    *final_error = err;
    *p_final_iteration = n_iterations;
    if (stagnated) return JMTX_RESULT_STAGNATED;
    return err < tolerance ? JMTX_RESULT_SUCCESS : JMTX_RESULT_NOT_CONVERGED;
}

jmtx_result jmtx_conjugate_gradient_crs_parallel(const jmtx_matrix_crs* mtx, const float* restrict y, float* restrict x,
                                                 const float tolerance, const float stagnation,
                                                 const uint32_t recalculation_interval, const uint32_t max_iterations,
                                                 uint32_t* p_final_iteration, float* restrict err_evolution,
                                                 float* final_error, float* restrict aux_vec1, float* restrict aux_vec2,
                                                 float* restrict aux_vec3)
{
#ifndef JMTX_NO_VERIFY_PARAMS
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
        return JMTX_RESULT_NULL_PARAM;
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
#endif

    const uint32_t n = mtx->base.rows;
    float* const r = aux_vec1;
    float* const p = aux_vec2;
    float* const Ap = aux_vec3;
    uint32_t n_iterations = 0, n_since_calc = 0;
    float alpha, beta;
    int update_residual, stagnated = 0;

    float rk_dp = 0;
    float pAp_dp = 0;
    float err = 0;
    float new_rk_dp = 0;
    float mag_y = 0;
    float prev_err = 0;

#pragma omp parallel default(none) shared(n, rk_dp, r, mtx, p, Ap, alpha, update_residual, n_iterations, err, x, y,\
                                          err_evolution, tolerance, beta, max_iterations, pAp_dp, new_rk_dp, mag_y,\
                                          stagnation, prev_err, stagnated, recalculation_interval, n_since_calc)
    {
#pragma omp for
        for (unsigned i = 0; i < n; ++i)
        {
            r[i] = y[i] - jmtx_matrix_crs_vector_multiply_row_raw(mtx, x, i);
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
                Ap[i] = jmtx_matrix_crs_vector_multiply_row_raw(mtx, p, i);
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
                    r[i] = y[i] - jmtx_matrix_crs_vector_multiply_row_raw(mtx, x, i);
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

    *final_error = err;
    *p_final_iteration = n_iterations;
    if (stagnated) return JMTX_RESULT_STAGNATED;
    return err < tolerance ? JMTX_RESULT_SUCCESS: JMTX_RESULT_NOT_CONVERGED;
}
