//
// Created by jan on 9.1.2024.
//
#include "../../../include/jmtx/float/solvers/induced_dimension_reduction_iteration.h"
#include "../../../include/jmtx/float/decompositions/dense_lu_decomposition.h"
#include "../../../include/jmtx/float/solvers/lu_solving.h"
#include "../matrices/dense_row_major_internal.h"
#include "../matrices/sparse_row_compressed_internal.h"

#include <math.h>

jmtx_result jmtx_solve_iterative_idrs_crs(const jmtx_matrix_crs *mtx, const float *restrict y, float *restrict x,
                                          float *restrict aux_vec1, float *restrict aux_vec2, uint32_t s,
                                          float *restrict aux_vec3, float *restrict aux_vec4, float *restrict aux_vec5,
                                          float *restrict aux_vec6, float *restrict aux_vec7,
                                          const jmtx_matrix_drm *p_mtx, jmtx_matrix_drm *aux_mtx1,
                                          jmtx_matrix_drm *aux_mtx2, jmtx_matrix_drm *aux_mtx3,
                                          jmtx_matrix_drm *aux_mtx4, jmtx_matrix_drm *aux_mtx5,
                                          jmtx_solver_arguments *args)
{
    (void)aux_vec7, (void)aux_mtx5;
    const uint32_t n = mtx->base.rows;

    //  Size n
    float *const r = aux_vec1;
    //  Size n
    float *const v = aux_vec2;
    //  Size s
    float *const pr = aux_vec3;
    //  Size s
    float *const c = aux_vec4;
    //  Size n
    float *const solve = aux_vec5;
    //  Size n
    float *const t = aux_vec6;
    float mag_y = 0, mag_r = 0;
    jmtx_matrix_crs_vector_multiply(mtx, x, r);
    for (uint32_t i = 0; i < n; ++i)
    {
        r[i] = (y[i] - r[i]);
        mag_y += y[i] * y[i];
        mag_r += r[i] * r[i];
    }
    mag_y = sqrtf(mag_y);
    mag_r = sqrtf(mag_r);
    float err = mag_r / mag_y;
    if (err < args->in_convergence_criterion)
    {
        args->out_last_error = err;
        args->out_last_iteration = 0;
        return JMTX_RESULT_SUCCESS;
    }

    //  Size s x n
    jmtx_matrix_drm *const DR = aux_mtx1;
    //  Size s x n
    jmtx_matrix_drm *const DX = aux_mtx2;
    uint32_t n_iter = 0;
    //  reduce the residual, while producing linearly independent vectors in Krylov subspace
    for (uint32_t i = 0; i < s; ++i)
    {
        //  New search direction
        jmtx_matrix_crs_vector_multiply(mtx, r, v);
        float vr_dp = 0, vv_dp = 0;
        for (uint32_t j = 0; j < n; ++j)
        {
            vr_dp += v[i] * r[i];
            vv_dp += v[i] * v[i];
        }

        const float omega = vr_dp / vv_dp;
        float *const row_dx = DX->values + n * i;
        float *const row_dr = DR->values + n * i;
        for (uint32_t j = 0; j < n; ++j)
        {
            const float dxj = omega * r[j];
            row_dx[j] = dxj;
            x[j] += dxj;
            const float drj = -omega * v[j];
            row_dr[j] = drj;
            r[j] += drj;
        }
        n_iter += 1;
    }

    //  Size s x s
    jmtx_matrix_drm *const pdr_mat = aux_mtx3;
    //  Prepare the P^T DR = pdr matrix
    for (uint32_t i = 0; i < s; ++i)
    {
        const float *p = p_mtx->values + n * i;
        for (uint32_t j = 0; j < s; ++j)
        {
            pdr_mat->values[j + i * s] = 0;
        }
        for (uint32_t k = 0; k < n; ++k)
        {
            for (uint32_t j = 0; j < s; ++j)
            {
                const float *dr = DR->values + j * n;
                pdr_mat->values[j + i * s] += p[k] * dr[k];
            }
        }
    }

    //  Size s x s
    jmtx_matrix_drm *const decomposed = aux_mtx4;

    //  DX and DR are now properly assembled for IDR
    float omega = 1.0f;
    for (;;)
    {
        for (uint32_t k = 0; k < s; ++k)
        {
            jmtx_matrix_drm_vector_multiply(p_mtx, r, pr);
            jmtx_decompose_lu_drm(pdr_mat, decomposed);
            jmtx_solve_direct_lu_drm(decomposed, pr, c, solve);
            float *const dr_new = DR->values + n * k;
            for (uint32_t i = 0; i < n; ++i)
            {
                solve[i] = c[0] * DR->values[i];
                for (uint32_t j = 1; j < s; ++j)
                {
                    const float *dr = DR->values + j * n;
                    solve[i] += c[j] * dr[i];
                }
                dr_new[i] = -r[i] + solve[i];
            }
            //  "dr_new" now contains - DR c

            if (k == 0)
            {
                jmtx_matrix_crs_vector_multiply(mtx, v, t);
                float num = 0, denom = 0;
                for (uint32_t i = 0; i < n; ++i)
                {
                    num += t[i] * v[i];
                    denom += t[i] * t[i];
                }
                omega = num / denom;
                for (uint32_t i = 0; i < n; ++i)
                {
                    dr_new[i] = dr_new[i] - omega * t[i];
                }
                //  Solve contains dr
                //  Update column of pdr_mat corresponding to dr
                for (uint32_t i = 0; i < s; ++i)
                {
                    const float *p = p_mtx->values + n * i;
                    pdr_mat->values[k + i * s] = 0;
                    for (uint32_t l = 0; l < n; ++l)
                    {
                        const float *dr = dr_new;
                        pdr_mat->values[k + i * s] += p[l] * dr[l];
                    }
                }
                //  Update residual
                mag_r = 0;
                for (uint32_t i = 0; i < n; ++i)
                {
                    r[i] = r[i] + dr_new[i];
                    mag_r += r[i] * r[i];
                }
                mag_r = sqrtf(mag_r);
                //  Compute DX
                float *dx_new = DX->values + k * n;
                for (uint32_t i = 0; i < n; ++i)
                {
                    float sum = omega * v[i];
                    for (uint32_t j = 0; j < s; ++j)
                    {
                        const float *dx = DX->values + j * n;
                        sum -= c[j] * dx[i];
                    }
                    dx_new[i] = sum;
                    x[i] += sum;
                }
            }
            else
            {
                //  Compute DX
                float *dx_new = DX->values + k * n;
                for (uint32_t i = 0; i < n; ++i)
                {
                    float sum = omega * v[i];
                    for (uint32_t j = 0; j < s; ++j)
                    {
                        const float *dx = DX->values + j * n;
                        sum -= c[j] * dx[i];
                    }
                    dx_new[i] = sum;
                    x[i] += sum;
                }
                jmtx_matrix_crs_vector_multiply(mtx, dx_new, dr_new);
                mag_r = 0;
                for (uint32_t i = 0; i < n; ++i)
                {
                    r[i] = r[i] + dr_new[i];
                    mag_r += r[i] * r[i];
                }
                mag_r = sqrtf(mag_r);
            }

            err = mag_r / mag_y;
            if (args->opt_error_evolution)
            {
                args->opt_error_evolution[n_iter] = err;
            }
            n_iter += 1;
            if (!isfinite(err) || err < args->in_convergence_criterion || n_iter == args->in_max_iterations)
            {
                break;
            }
        }
        if (!isfinite(err) || err < args->in_convergence_criterion || n_iter == args->in_max_iterations)
        {
            break;
        }
    }

    args->out_last_iteration = n_iter;
    args->out_last_error = err;

    if (!isfinite(err) || err > args->in_convergence_criterion)
    {
        return JMTX_RESULT_NOT_CONVERGED;
    }

    return JMTX_RESULT_SUCCESS;
}
