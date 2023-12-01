// Automatically generated from source/float/solvers/bicgstab_iteration.c on Fri Dec  1 06:43:01 2023
//
// Created by jan on 17.6.2022.
//

#include <math.h>
#include <stdio.h>
#include "bicgstab_iteration.h"
#include "../matrices/sparse_row_compressed_internal.h"
#include <pthread.h>
#include <stdatomic.h>
#include <assert.h>


jmtx_result jmtxd_bicgstab_crs(
        const jmtxd_matrix_crs* mtx, const double* y, double* x, double convergence_dif,
        uint32_t n_max_iter, uint32_t* p_iter, double* p_final_error, double* p_error,
        const jmtx_allocator_callbacks* allocator_callbacks)
{
    if (!mtx)
    {
//        REPORT_ERROR_MESSAGE("Matrix pointer was null");
//        LEAVE_FUNCTION();
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTXD_TYPE_CRS)
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
    double* const joint_buffer = allocator_callbacks->alloc(allocator_callbacks->state, (n + n + n + n + n + n) * sizeof(*joint_buffer));
    if (!joint_buffer)
    {
//        CALLOC_FAILED((n + n + n + n + n + n) * sizeof(*joint_buffer));
//        LEAVE_FUNCTION();
        return JMTX_RESULT_BAD_ALLOC;
    }
    double rho_new = 1, rho_old = 1, a = 1, omega_old = 1, omega_new = 1;

    double* const p = joint_buffer;
    double* const v = joint_buffer + n;
    double* const r0 = joint_buffer + 2 * n;
    double* const r_new = joint_buffer + 3 * n;
    double* const s = joint_buffer + 4 * n;
    double* const t = joint_buffer + 5 * n;

    memcpy(x, y, n * sizeof(double));
    for (uint32_t i = 0; i < n; ++i)
    {
        double d = jmtxd_matrix_crs_get_entry(mtx, i, i);
        x[i] /= d;
    }

    jmtxd_matrix_crs_vector_multiply(mtx, x, r0);
    for (uint32_t i = 0; i < n; ++i)
    {
        r0[i] = y[i] - r0[i];
    }
    memcpy(r_new, r0, sizeof*r0 * n);

    double err_rms;
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
        const double beta = rho_new / rho_old * a / omega_old;
        for (uint32_t i = 0; i < n; ++i)
        {
            p[i] = r_new[i] + beta * (p[i] - omega_old * v[i]);
        }
        jmtxd_matrix_crs_vector_multiply(mtx, p, v);
        double tmp = 0;
        for (uint32_t i = 0; i < n; ++i)
        {
            tmp += r0[i] * v[i];
        }
        a = rho_new / tmp;
        for (uint32_t i = 0; i < n; ++i)
        {
            const double dx = a * p[i];
//            err_rms += dx * dx;
            x[i] = x[i] + dx;
        }
        err_rms = 0;
        for (uint32_t i = 0; i < n; ++i)
        {
            double val = jmtxd_matrix_crs_vector_multiply_row(mtx, x, i);
            val -= y[i];
            err_rms += val * val;
        }

        err_rms = sqrt(err_rms) / (double)n;
        if (p_error)
        {
            p_error[n] = err_rms;
        }
        if (err_rms < convergence_dif) break;
        for (uint32_t i = 0; i < n; ++i)
        {
            s[i] = r_new[i] - a * v[i];
        }
        jmtxd_matrix_crs_vector_multiply(mtx, s, t);
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
            const double new_x = x[i] + omega_new * s[i];
            x[i] = new_x;
        }
        err_rms = 0;
        for (uint32_t i = 0; i < n; ++i)
        {
            double val = jmtxd_matrix_crs_vector_multiply_row(mtx, x, i);
            val -= y[i];
            err_rms += val * val;
        }

        err_rms = sqrt(err_rms) / (double)n;
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


struct bicgstab_crs_args
{
    const jmtxd_matrix_crs* mtx;
    double* x;
    const double* y;
    double* p_common;
    double* p_common2;
    double* p_err;
    double* p_error_evolution;
    uint32_t* p_iter;
    uint32_t max_iter;
    double req_error;
    uint32_t n;
    uint32_t id;
    uint32_t thrd_count;
    uint32_t* done;
    double* joint_buffer;
    pthread_barrier_t* barrier;
};

static void* bicgstab_crs_thrd_fn(void* param)
{
    const struct bicgstab_crs_args* args = param;
//    THREAD_BEGIN("Worker %s (%d/%d)", __func__, args->id, args->thrd_count);
    double* const x = args->x;
    double* const p = args->joint_buffer;
    double* const v = args->joint_buffer + args->n;
    double* const r0 = args->joint_buffer + 2 * args->n;
    double* const r_new = args->joint_buffer + 3 * args->n;
    double* const s = args->joint_buffer + 4 * args->n;
    double* const t = args->joint_buffer + 5 * args->n;
    uint32_t iter_count = 0;
    double err = 0;

    const uint32_t complete_chunk = args->n / args->thrd_count;
    const uint32_t remaining_chunk = args->n % args->thrd_count;
    uint32_t work_chunk = complete_chunk;
    if (remaining_chunk)
    {
        work_chunk += args->id < remaining_chunk;
    }
    uint32_t begin_work = complete_chunk * args->id;
    if (remaining_chunk)
    {
        begin_work += (remaining_chunk > args->id ? args->id : remaining_chunk);
    }
    const uint32_t begin = begin_work;
    const uint32_t end = begin_work + work_chunk;
    double rho_new = 1, rho_old = 1, a = 1, omega_new = 1, omega_old = 1;

    while (!*args->done && iter_count < args->max_iter)
    {

        omega_old = omega_new;
        rho_old = rho_new;


        // compute new rho
        rho_new = 0;

        for (uint32_t i = begin; i < end; ++i)
        {
            rho_new += r0[i] * r_new[i];
        }
        args->p_common[args->id] = rho_new;
        //  Barrier to compute rho_new
        pthread_barrier_wait(args->barrier);
        rho_new = 0;
        for (uint32_t i = 0; i < args->thrd_count; ++i)
        {
            rho_new += args->p_common[i];
        }


        const double beta = rho_new / rho_old * a / omega_old;

        for (uint32_t i = begin; i < end; ++i)
        {
            p[i] = r_new[i] + beta * (p[i] - omega_old * v[i]);
        }

        //  Barrier to perform matrix multiplication
        pthread_barrier_wait(args->barrier);
        for (uint32_t i = begin; i < end; ++i)
        {
            double val = 0, * elements;
            uint32_t* indices;
            uint32_t n =jmtxd_matrix_crs_get_row(args->mtx, i, &indices, &elements);
            for (uint32_t j = 0; j < n; ++j)
            {
                val += elements[j] * p[indices[j]];
            }
            v[i] = val;
        }

        //  No barrier yet, because only values of v between begin and end are needed
        double mag = 0;
        for (uint32_t i = begin; i < end; ++i)
        {
            mag += r0[i] * v[i];
        }
        args->p_common2[args->id] = mag;
        pthread_barrier_wait(args->barrier);
        mag = 0;
        for (uint32_t i = 0; i < args->thrd_count; ++i)
        {
            mag += args->p_common2[i];
        }
        a = rho_new / mag;


        for (uint32_t i = begin; i < end; ++i)
        {
            x[i] = x[i] + a * p[i];
        }

        //  Barrier to compute residual with matrix operations
        pthread_barrier_wait(args->barrier);
        err = 0;
        for (uint32_t i = begin; i < end; ++i)
        {
            double y0 = jmtxd_matrix_crs_vector_multiply_row(args->mtx, x, i);
            y0 -= args->y[i];
            err += y0 * y0;
        }
        args->p_common[args->id] = err;
        //  Barrier to compute total residual
        pthread_barrier_wait(args->barrier);
        err = 0;
        for (uint32_t i = 0; i < args->thrd_count; ++i)
        {
            err += args->p_common[i];
        }
        err = sqrt(err) / (double)args->n;

        if (args->p_error_evolution)
        {
            args->p_error_evolution[iter_count] = err;
        }
        if (err < args->req_error) break;

        for (uint32_t i = begin; i < end; ++i)
        {
            s[i] = r_new[i] - a * v[i];
        }

        //  Barrier to perform matrix operation
        pthread_barrier_wait(args->barrier);
        double o = 0;
        mag = 0;
        for (uint32_t i = begin; i < end; ++i)
        {
            double val = 0, * elements;
            uint32_t* indices;
            uint32_t n = jmtxd_matrix_crs_get_row(args->mtx, i, &indices, &elements);
            for (uint32_t j = 0; j < n; ++j)
            {
                val += elements[j] * s[indices[j]];
            }
            t[i] = val;
            o += val * s[i];
            mag += val * val;
        }
        args->p_common[args->id] = o;
        args->p_common2[args->id] = mag;
        //  Barrier to compute omega new
        pthread_barrier_wait(args->barrier);
        o = 0;
        mag = 0;
        for (uint32_t i = 0; i < args->thrd_count; ++i)
        {
            o += args->p_common[i];
            mag += args->p_common2[i];
        }
        omega_new = o / mag;

        for (uint32_t i = begin; i < end; ++i)
        {
            x[i] = x[i] + omega_new * s[i];
        }

        err = 0;
        //  Barrier to compute residual with matrix operations
        pthread_barrier_wait(args->barrier);
        for (uint32_t i = begin; i < end; ++i)
        {
            double y0 = jmtxd_matrix_crs_vector_multiply_row(args->mtx, x, i);
            y0 -= args->y[i];
            err += y0 * y0;
        }
        args->p_common2[args->id] = err;
        //  Barrier to compute total residual
        pthread_barrier_wait(args->barrier);
        err = 0;
        for (uint32_t i = 0; i < args->thrd_count; ++i)
        {
            err += args->p_common2[i];
        }
        err = sqrt(err) / (double)args->n;
        if (args->p_error_evolution)
        {
            args->p_error_evolution[iter_count] = err;
        }
        if (err < args->req_error) break;

        for (uint32_t i = begin; i < end; ++i)
        {
            r_new[i] = s[i] - omega_new * t[i];
        }

        iter_count += 1;
    }
    *args->p_iter = iter_count;
    *args->p_err = err;
    pthread_barrier_wait(args->barrier);
    args->p_common[args->id] = err;
    pthread_barrier_wait(args->barrier);
    for (uint32_t i = 0; i < args->thrd_count; ++i)
    {
        assert(err != args->p_common[i]);
        {
//            REPORT_ERROR_MESSAGE("Error computed by thread %u does not match error by this thread (%g vs %g)", i, err, args->p_common[i]);
        }
    }
//    THREAD_END;
    return 0;
}

jmtx_result jmtxd_bicgstab_crs_mt(
        const jmtxd_matrix_crs* mtx, const double* y, double* x, double convergence_dif,
        uint32_t n_max_iter, uint32_t* p_iter, double* p_final_error, double* p_error,
        const jmtx_allocator_callbacks* allocator_callbacks, uint32_t n_thrds)
{

    if (!mtx)
    {
//        REPORT_ERROR_MESSAGE("Matrix pointer was null");
//        LEAVE_FUNCTION();
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTXD_TYPE_CRS)
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

    if (n_thrds < 2)
    {
        const jmtx_result result = (
                jmtxd_bicgstab_crs(
                        mtx, y, x, convergence_dif, n_max_iter, p_iter, p_final_error, p_error, allocator_callbacks));
//        LEAVE_FUNCTION();
        return result;
    }
    const uint32_t n = mtx->base.rows;
    pthread_t* const threads = allocator_callbacks->alloc(allocator_callbacks->state, (n_thrds) * sizeof*threads);
    if (!threads)
    {
//        CALLOC_FAILED(n_thrds * sizeof*threads);
//        LEAVE_FUNCTION();
        return JMTX_RESULT_BAD_ALLOC;
    }
    double* const common_buffer = allocator_callbacks->alloc(allocator_callbacks->state, (2 * n_thrds) * sizeof*common_buffer);
    if (!common_buffer)
    {
        allocator_callbacks->free(allocator_callbacks->state, threads);
//        CALLOC_FAILED(n_thrds * sizeof*common_buffer);
//        LEAVE_FUNCTION();
        return JMTX_RESULT_BAD_ALLOC;
    }
    struct bicgstab_crs_args* const args = allocator_callbacks->alloc(allocator_callbacks->state, (n_thrds) * sizeof*args);
    if (!args)
    {
        allocator_callbacks->free(allocator_callbacks->state, common_buffer);
        allocator_callbacks->free(allocator_callbacks->state, threads);
//        CALLOC_FAILED(n_thrds * sizeof*args);
//        LEAVE_FUNCTION();
        return JMTX_RESULT_BAD_ALLOC;
    }
    pthread_barrier_t barrier;
    pthread_barrier_init(&barrier, NULL, n_thrds);
    double* const joint_buffer = allocator_callbacks->alloc(allocator_callbacks->state, (n + n + n + n + n + n) * sizeof(*joint_buffer));
    if (!joint_buffer)
    {
        pthread_barrier_destroy(&barrier);
        allocator_callbacks->free(allocator_callbacks->state, args);
        allocator_callbacks->free(allocator_callbacks->state, common_buffer);
        allocator_callbacks->free(allocator_callbacks->state, threads);
//        CALLOC_FAILED((n + n + n + n + n + n) * sizeof*joint_buffer);
//        LEAVE_FUNCTION();
        return JMTX_RESULT_BAD_ALLOC;
    }


//    scalar_t* const p = joint_buffer;
//    scalar_t* const v = joint_buffer + n;
    double* const r0 = joint_buffer + 2 * n;
    double* const r_new = joint_buffer + 3 * n;
//    scalar_t* const s = joint_buffer + 4 * n;
//    scalar_t* const t = joint_buffer + 5 * n;

    //  Initial guess assumes that mtx is a diagonal matrix
    memcpy(x, y, n * sizeof(double));
    for (uint32_t i = 0; i < n; ++i)
    {
        double d = jmtxd_matrix_crs_get_entry(mtx, i, i);
        x[i] /= d;
    }

    jmtxd_matrix_crs_vector_multiply(mtx, x, r0);
    for (uint32_t i = 0; i < n; ++i)
    {
        r0[i] = y[i] - r0[i];
    }
    memcpy(r_new, r0, sizeof*r0 * n);


    uint32_t done = 0;
    double err_rms;
    uint32_t iter_count = 0;
    for (uint32_t i = 0; i < n_thrds; ++i)
    {
        args[i] = (struct bicgstab_crs_args)
                {
            .barrier = &barrier,
            .p_common = common_buffer,
            .p_common2 = common_buffer + n_thrds,
            .id = i,
            .n = n,
            .joint_buffer = joint_buffer,
            .thrd_count = n_thrds,
            .done = &done,
            .p_iter = &iter_count,
            .p_err = &err_rms,
            .req_error = convergence_dif,
            .max_iter = n_max_iter,
            .mtx = mtx,
            .x = x,
            .y = y,
            .p_error_evolution = p_error,
                };
        if (pthread_create(threads + i, NULL, bicgstab_crs_thrd_fn, args + i))
        {
            done = 1;
//            REPORT_ERROR_MESSAGE("Could not create a new worker thread (%d out of %d), reason: %s", i + 1, n_thrds,
//                                 strerror(errno));
            for (uint32_t j = 0; j < i; ++i)
            {
                pthread_cancel(threads[i]);
            }
            allocator_callbacks->free(allocator_callbacks->state, joint_buffer);
            pthread_barrier_destroy(&barrier);
            allocator_callbacks->free(allocator_callbacks->state, args);
            allocator_callbacks->free(allocator_callbacks->state, common_buffer);
            allocator_callbacks->free(allocator_callbacks->state, threads);
            return JMTX_RESULT_BAD_THREAD;
        }
    }


    for (uint32_t i = 0; i < n_thrds; ++i)
    {
        pthread_join(threads[i], NULL);
    }

    allocator_callbacks->free(allocator_callbacks->state, joint_buffer);
    pthread_barrier_destroy(&barrier);
    allocator_callbacks->free(allocator_callbacks->state, args);
    allocator_callbacks->free(allocator_callbacks->state, common_buffer);
    allocator_callbacks->free(allocator_callbacks->state, threads);
    if (p_iter) *p_iter = iter_count;
    if (p_final_error) *p_final_error = err_rms;
//    LEAVE_FUNCTION();
    return iter_count == n_max_iter ? JMTX_RESULT_NOT_CONVERGED : JMTX_RESULT_SUCCESS;
}

