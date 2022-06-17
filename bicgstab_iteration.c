//
// Created by jan on 17.6.2022.
//

#include <math.h>
#include <stdio.h>
#include "bicgstab_iteration.h"

mtx_res_t bicgstab_crs(
        const CrsMatrix* mtx, const scalar_t* y, scalar_t* x, scalar_t convergence_dif, uint n_max_iter, uint* n_iter)
{
    const uint n = mtx->rows;
    scalar_t* const joint_buffer = calloc(n + n + n + n + n + n + n, sizeof(*joint_buffer));
    if (!joint_buffer) return -1;
    scalar_t rho_new = 1, rho_old = 1, a = 1, omega_old = 1, omega_new = 1;

    scalar_t* const p = joint_buffer;
    scalar_t* const v = joint_buffer + n;
    scalar_t* const r0 = joint_buffer + 2 * n;
    scalar_t* const r_new = joint_buffer + 3 * n;
    scalar_t* const s = joint_buffer + 4 * n;
    scalar_t* const t = joint_buffer + 5 * n;
    scalar_t* const h = joint_buffer + 6 * n;

    memcpy(x, y, n * sizeof(scalar_t));
    for (uint i = 0; i < n; ++i)
    {
        scalar_t d;
        matrix_crs_get_element(mtx, i, i , &d);
        x[i] /= d;
    }

    matrix_crs_vector_multiply(mtx, x, r0);
    for (uint i = 0; i < n; ++i)
    {
        r0[i] = y[i] - r0[i];
    }
    memcpy(r_new, r0, sizeof*r0 * n);

    scalar_t err_rms;
    uint iter_count = 0;
    do
    {
        rho_old = rho_new;
        omega_old = omega_new;
        err_rms = 0;

        rho_new = 0;
        for (uint i = 0; i < n; ++i)
        {
            rho_new += r0[i] * r_new[i];
        }
        const scalar_t beta = rho_new / rho_old * a / omega_old;
        for (uint i = 0; i < n; ++i)
        {
            p[i] = r_new[i] + beta * (p[i] - omega_old * v[i]);
        }
        matrix_crs_vector_multiply(mtx, p, v);
        scalar_t tmp = 0;
        for (uint i = 0; i < n; ++i)
        {
            tmp += r0[i] * v[i];
        }
        a = rho_new / tmp;
        for (uint i = 0; i < n; ++i)
        {
            const scalar_t dx = a * p[i];
            err_rms += dx * dx;
            h[i] = x[i] + dx;
        }
        err_rms = sqrtf(err_rms);
        if (err_rms < convergence_dif) break;
        for (uint i = 0; i < n; ++i)
        {
            s[i] = r_new[i] - a * v[i];
        }
        matrix_crs_vector_multiply(mtx, s, t);
        omega_new = 0;
        tmp = 0;
        for (uint i = 0; i < n; ++i)
        {
            omega_new += t[i] * s[i];
            tmp += t[i] * t[i];
        }
        omega_new /= tmp;
        err_rms = 0;
        for (uint i = 0; i < n; ++i)
        {
            const scalar_t new_x = h[i] + omega_new * s[i];
            const scalar_t dx = new_x - x[i];
            err_rms += dx * dx;
            x[i] = new_x;
        }
        err_rms = sqrtf(err_rms);
//        printf("RMS error on iteration %u is %g\n", iter_count, err_rms);
        if ((err_rms) < convergence_dif) break;

        for (uint i = 0; i < n; ++i)
        {
            r_new[i] = s[i] - omega_new * t[i];
        }

        ++iter_count;
    } while(iter_count < n_max_iter);
    free(joint_buffer);
    if (n_iter) *n_iter = iter_count;
    return 0;
}
