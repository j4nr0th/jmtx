//
// Created by jan on 17.6.2022.
//

#include <math.h>
#include <stdio.h>
#include "bicgstab_iteration.h"
#include <pthread.h>
#include <stdatomic.h>

mtx_res_t bicgstab_crs(
        const CrsMatrix* mtx, const scalar_t* y, scalar_t* x, scalar_t convergence_dif, uint n_max_iter, uint* p_iter,
        scalar_t* p_error)
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
    if (p_iter) *p_iter = iter_count;
    if (p_error) *p_error = err_rms;
    return 0;
}


struct bicgstab_crs_args
{
    const CrsMatrix* mtx;
    scalar_t* rho_old, *rho_new, *a, *omega_old, *omega_new;
    scalar_t* x;
    scalar_t* p_common;
    scalar_t* p_common2;
    scalar_t* p_err;
    uint* p_iter;
    uint max_iter;
    scalar_t req_error;
    uint n;
    uint id;
    uint thrd_count;
    uint* done;
    scalar_t* joint_buffer;
    pthread_barrier_t* barrier;
    atomic_flag* once_flag;
};

static void* bicgstab_crs_thrd_fn(void* param)
{
    struct bicgstab_crs_args* args = param;
    scalar_t* const x = args->x;
    scalar_t* const p = args->joint_buffer;
    scalar_t* const v = args->joint_buffer + args->n;
    scalar_t* const r0 = args->joint_buffer + 2 * args->n;
    scalar_t* const r_new = args->joint_buffer + 3 * args->n;
    scalar_t* const s = args->joint_buffer + 4 * args->n;
    scalar_t* const t = args->joint_buffer + 5 * args->n;
    scalar_t* const h = args->joint_buffer + 6 * args->n;
    uint iter_count = 0;
    _Bool flag_v = 0;
    scalar_t err;

    const uint complete_chunk = args->n / args->thrd_count;
    const uint remaining_chunk = args->n % args->thrd_count;
    uint work_chunk = complete_chunk;
    if (remaining_chunk)
    {
        work_chunk += args->id < remaining_chunk;
    }
    uint begin_work = complete_chunk * args->id;
    if (remaining_chunk)
    {
        begin_work += (remaining_chunk > args->id ? args->id : remaining_chunk);
    }
    const uint begin = begin_work;
    const uint end = begin_work + work_chunk;


    while (!*args->done)
    {
        flag_v = atomic_flag_test_and_set(args->once_flag);
        pthread_barrier_wait(args->barrier);
        if (!flag_v)
        {
            *args->rho_old = *args->rho_new;
            *args->omega_old = *args->omega_new;
            atomic_flag_clear(args->once_flag);
        }
//        pthread_barrier_wait(args->barrier);


        // compute new rho
        scalar_t res_rho = 0;

        for (uint i = begin; i < end; ++i)
        {
            res_rho += r0[i] * r_new[i];
        }
        args->p_common[args->id] = res_rho;
        flag_v = atomic_flag_test_and_set(args->once_flag);
        pthread_barrier_wait(args->barrier);
        if (!flag_v)
        {
            res_rho = 0;
            for (uint i = 0; i < args->thrd_count; ++i)
            {
                res_rho += args->p_common[i];
            }
            *args->rho_new = res_rho;
            atomic_flag_clear(args->once_flag);
        }
        pthread_barrier_wait(args->barrier);
        const scalar_t rho_new = *args->rho_new;
        const scalar_t rho_old = *args->rho_old;
        const scalar_t omega_old = *args->omega_old;

        const scalar_t beta = rho_new / rho_old * *args->a / omega_old;

        for (uint i = begin; i < end; ++i)
        {
            p[i] = r_new[i] + beta * (p[i] - omega_old * v[i]);
        }

        pthread_barrier_wait(args->barrier);
//        matrix_crs_vector_multiply(mtx, p, v);
        for (uint i = begin; i < end; ++i)
        {
            scalar_t val = 0, * elements;
            uint n, * indices;

            matrix_crs_get_row(args->mtx, i, &n, &indices, &elements);
            for (uint j = 0; j < n; ++j)
            {
                val += elements[j] * p[indices[j]];
            }
            v[i] = val;
        }

        //  No barrier yet, because only values of v between begin and end are needed
        scalar_t mag = 0;
        for (uint i = begin; i < end; ++i)
        {
            mag += r0[i] * v[i];
        }
        args->p_common[args->id] = mag;
        flag_v = atomic_flag_test_and_set(args->once_flag);
        pthread_barrier_wait(args->barrier);
        if (!flag_v)
        {
            mag = 0;
            for (uint i = 0; i < args->thrd_count; ++i)
            {
                mag += args->p_common[i];
            }
            *args->a = rho_new / mag;
            atomic_flag_clear(args->once_flag);
        }
        pthread_barrier_wait(args->barrier);

        const scalar_t a = *args->a;
        err = 0;
        for (uint i = begin; i < end; ++i)
        {
            const scalar_t dx = a * p[i];
            err += dx * dx;
            h[i] = x[i] + dx;
        }
        args->p_common[args->id] = err;

        flag_v = atomic_flag_test_and_set(args->once_flag);
        pthread_barrier_wait(args->barrier);
        if (!flag_v)
        {
            err = 0;
            for (uint i = 0; i < args->thrd_count; ++i)
            {
                err += args->p_common[i];
            }
            err = sqrtf(err);
            if (err < args->req_error)
            {
                *args->p_iter = iter_count;
                *args->p_err = err;
                *args->done = 1;
            }
            atomic_flag_clear(args->once_flag);
        }
        pthread_barrier_wait(args->barrier);
        if (*args->done) break;

        for (uint i = begin; i < end; ++i)
        {
            s[i] = r_new[i] - a * v[i];
        }
        pthread_barrier_wait(args->barrier);

        scalar_t o = 0;
        mag = 0;
        for (uint i = begin; i < end; ++i)
        {
            scalar_t val = 0, * elements;
            uint n, * indices;

            matrix_crs_get_row(args->mtx, i, &n, &indices, &elements);
            for (uint j = 0; j < n; ++j)
            {
                val += elements[j] * s[indices[j]];
            }
            t[i] = val;
            o += val * s[i];
            mag += val * val;
        }
        args->p_common[args->id] = o;
        args->p_common2[args->id] = mag;
        flag_v = atomic_flag_test_and_set(args->once_flag);
        pthread_barrier_wait(args->barrier);
        if (!flag_v)
        {
            o = 0, mag = 0;
            for (uint i = 0; i < args->thrd_count; ++i)
            {
                o += args->p_common[i];
                mag += args->p_common2[i];
            }
            *args->omega_new = o / mag;
            atomic_flag_clear(args->once_flag);
        }
        pthread_barrier_wait(args->barrier);
        const scalar_t omega_new = *args->omega_new;
        err = 0;
        for (uint i = begin; i < end; ++i)
        {
            const scalar_t new_x = h[i] + omega_new * s[i];
            const scalar_t dx = new_x - x[i];
            err += dx * dx;
            x[i] = new_x;
        }
        args->p_common[args->id] = err;
        flag_v = atomic_flag_test_and_set(args->once_flag);
        pthread_barrier_wait(args->barrier);
        if (!flag_v)
        {
            err = 0;
            for (uint i = 0; i < args->thrd_count; ++i)
            {
                err += args->p_common[i];
            }
            err = sqrtf(err);
            if (err < args->req_error)
            {
                *args->p_err = err;
                *args->p_iter = iter_count;
                *args->done = 1;
            }
            atomic_flag_clear(args->once_flag);
        }
        pthread_barrier_wait(args->barrier);
        if (*args->done) break;



        for (uint i = begin; i < end; ++i)
        {
            r_new[i] = s[i] - omega_new * t[i];
        }

        iter_count += 1;
//        pthread_barrier_wait(args->barrier);
    }

    pthread_barrier_wait(args->barrier);

    return 0;
}

mtx_res_t bicgstab_crs_mt(
        const CrsMatrix* mtx, const scalar_t* y, scalar_t* x, scalar_t convergence_dif, uint n_max_iter, uint* p_iter,
        scalar_t* p_error, uint n_thrds)
{
    if (n_thrds == 0)
    {
        return bicgstab_crs(mtx, y, x, convergence_dif, n_max_iter, p_iter, p_error);
    }
    const uint n = mtx->rows;
    pthread_t* const threads = calloc(n_thrds, sizeof*threads);
    if (!threads) return -1;
    atomic_flag once_flag = ATOMIC_FLAG_INIT;
    scalar_t* const common_buffer = calloc(2 * n_thrds, sizeof*common_buffer);
    if (!common_buffer)
    {
        free(threads);
        return -1;
    }
    struct bicgstab_crs_args* const args = calloc(n_thrds, sizeof*args);
    if (!args)
    {
        free(common_buffer);
        free(threads);
        return -1;
    }
    pthread_barrier_t barrier;
    pthread_barrier_init(&barrier, NULL, n_thrds);
    scalar_t* const joint_buffer = calloc(n + n + n + n + n + n + n, sizeof(*joint_buffer));
    if (!joint_buffer)
    {
        pthread_barrier_destroy(&barrier);
        free(args);
        free(common_buffer);
        free(threads);
        return -1;
    }

    scalar_t rho_new = 1, rho_old = 1, a = 1, omega_old = 1, omega_new = 1;

//    scalar_t* const p = joint_buffer;
//    scalar_t* const v = joint_buffer + n;
    scalar_t* const r0 = joint_buffer + 2 * n;
    scalar_t* const r_new = joint_buffer + 3 * n;
//    scalar_t* const s = joint_buffer + 4 * n;
//    scalar_t* const t = joint_buffer + 5 * n;
//    scalar_t* const h = joint_buffer + 6 * n;

    //  Initial guess assumes that mtx is a diagonal matrix
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

    // TODO: create threads and send them on their merry way
    uint done = 0;
    scalar_t err_rms;
    uint iter_count = 0;
    for (uint i = 0; i < n_thrds; ++i)
    {
        args[i] = (struct bicgstab_crs_args)
                {
            .barrier = &barrier,
            .once_flag = &once_flag,
            .p_common = common_buffer,
            .p_common2 = common_buffer + n_thrds,
            .id = i,
            .n = n,
            .joint_buffer = joint_buffer,
            .thrd_count = n_thrds,
            .done = &done,
            .omega_new = &omega_new,
            .omega_old = &omega_old,
            .a = &a,
            .rho_new = &rho_new,
            .rho_old = &rho_old,
            .p_iter = &iter_count,
            .p_err = &err_rms,
            .req_error = convergence_dif,
            .max_iter = n_max_iter,
            .mtx = mtx,
            .x = x
                };
        if (pthread_create(threads + i, NULL, bicgstab_crs_thrd_fn, args + i))
        {
            done = 1;
            for (uint j = 0; j < i; ++i)
            {
                pthread_cancel(threads[i]);
            }
            free(joint_buffer);
            pthread_barrier_destroy(&barrier);
            free(args);
            free(common_buffer);
            free(threads);
            return -1;
        }
    }


    for (uint i = 0; i < n_thrds; ++i)
    {
        pthread_join(threads[i], NULL);
    }

    free(joint_buffer);
    pthread_barrier_destroy(&barrier);
    free(args);
    free(common_buffer);
    free(threads);
    if (p_iter) *p_iter = iter_count;
    if (p_error) *p_error = err_rms;
    return 0;
}

