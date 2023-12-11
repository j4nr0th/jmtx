//
// Created by jan on 17.6.2022.
//

#include <math.h>
#include <stdio.h>
#include "bicgstab_iteration.h"
#include "../matrices/sparse_row_compressed_internal.h"
#include <assert.h>


jmtx_result jmtx_bicgstab_crs(
        const jmtx_matrix_crs* mtx, const float* y, float* x, float convergence_dif,
        uint32_t n_max_iter, uint32_t* p_iter, float* p_final_error, float* p_error,
        const jmtx_allocator_callbacks* allocator_callbacks)
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

    if (!allocator_callbacks)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }
    else if (!allocator_callbacks->alloc || !allocator_callbacks->realloc || !allocator_callbacks->free)
    {
        return JMTX_RESULT_BAD_PARAM;
    }
    
    
    const uint32_t n = mtx->base.rows;
    float* const joint_buffer = allocator_callbacks->alloc(allocator_callbacks->state, (n + n + n + n + n + n) * sizeof(*joint_buffer));
    if (!joint_buffer)
    {
//        CALLOC_FAILED((n + n + n + n + n + n) * sizeof(*joint_buffer));
//        LEAVE_FUNCTION();
        return JMTX_RESULT_BAD_ALLOC;
    }
    float rho_new = 1, rho_old = 1, a = 1, omega_old = 1, omega_new = 1;

    float* const p = joint_buffer;
    float* const v = joint_buffer + n;
    float* const r0 = joint_buffer + 2 * n;
    float* const r_new = joint_buffer + 3 * n;
    float* const s = joint_buffer + 4 * n;
    float* const t = joint_buffer + 5 * n;

    memcpy(x, y, n * sizeof(float));
    for (uint32_t i = 0; i < n; ++i)
    {
        float d = jmtx_matrix_crs_get_entry(mtx, i, i);
        x[i] /= d;
    }

    jmtx_matrix_crs_vector_multiply(mtx, x, r0);
    for (uint32_t i = 0; i < n; ++i)
    {
        r0[i] = y[i] - r0[i];
    }
    memcpy(r_new, r0, sizeof*r0 * n);

    float err_rms;
    uint32_t iter_count = 0;
    do
    {
        rho_old = rho_new;
        omega_old = omega_new;

        rho_new = 0;
        for (uint32_t i = 0; i < n; ++i)
        {
            rho_new += r0[i] * r_new[i];
        }
        const float beta = rho_new / rho_old * a / omega_old;
        for (uint32_t i = 0; i < n; ++i)
        {
            p[i] = r_new[i] + beta * (p[i] - omega_old * v[i]);
        }
        jmtx_matrix_crs_vector_multiply(mtx, p, v);
        float tmp = 0;
        for (uint32_t i = 0; i < n; ++i)
        {
            tmp += r0[i] * v[i];
        }
        a = rho_new / tmp;
        for (uint32_t i = 0; i < n; ++i)
        {
            const float dx = a * p[i];
//            err_rms += dx * dx;
            x[i] = x[i] + dx;
        }
        err_rms = 0;
        for (uint32_t i = 0; i < n; ++i)
        {
            float val = jmtx_matrix_crs_vector_multiply_row(mtx, x, i);
            val -= y[i];
            err_rms += val * val;
        }

        err_rms = sqrtf(err_rms) / (float)n;
        if (p_error)
        {
            p_error[n] = err_rms;
        }
        if (err_rms < convergence_dif) break;
        for (uint32_t i = 0; i < n; ++i)
        {
            s[i] = r_new[i] - a * v[i];
        }
        jmtx_matrix_crs_vector_multiply(mtx, s, t);
        omega_new = 0;
        tmp = 0;
        for (uint32_t i = 0; i < n; ++i)
        {
            omega_new += t[i] * s[i];
            tmp += t[i] * t[i];
        }
        omega_new /= tmp;
        for (uint32_t i = 0; i < n; ++i)
        {
            const float new_x = x[i] + omega_new * s[i];
            x[i] = new_x;
        }
        err_rms = 0;
        for (uint32_t i = 0; i < n; ++i)
        {
            float val = jmtx_matrix_crs_vector_multiply_row(mtx, x, i);
            val -= y[i];
            err_rms += val * val;
        }

        err_rms = sqrtf(err_rms) / (float)n;
//        printf("RMS error on iteration %u is %g\n", iter_count, err_rms);
        if (p_error)
        {
            p_error[n] = err_rms;
        }
        if ((err_rms) < convergence_dif) break;

        for (uint32_t i = 0; i < n; ++i)
        {
            r_new[i] = s[i] - omega_new * t[i];
        }

        ++iter_count;
    } while(iter_count < n_max_iter);
    allocator_callbacks->free(allocator_callbacks->state, joint_buffer);
    if (p_iter) *p_iter = iter_count;
    if (p_final_error) *p_final_error = err_rms;
//    LEAVE_FUNCTION();
    return iter_count == n_max_iter ? JMTX_RESULT_NOT_CONVERGED : JMTX_RESULT_SUCCESS;
}

