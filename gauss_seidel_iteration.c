//
// Created by jan on 16.6.2022.
//

#include "gauss_seidel_iteration.h"

#include <pthread.h>
#include <stdatomic.h>
#include <stdio.h>
#include <math.h>
#include "errors.h"
//#include <unistd.h>
//#include <errno.h>

#undef gauss_seidel_crs
#undef gauss_seidel_crs_mt

mtx_res_t gauss_seidel_crs(
        const CrsMatrix* mtx, const scalar_t* y, scalar_t* x, scalar_t convergence_dif, uint n_max_iter, uint* p_iter,
        scalar_t* p_error)
{
#ifdef MTX_MATRIX_CHECKS
    if (!mtx)
    {
        REPORT_ERROR_MESSAGE("Matrix pointer was null");
        LEAVE_FUNCTION();
        return mtx_bad_param;
    }
    if (mtx->type != mtx_type_crs)
    {
        REPORT_ERROR_MESSAGE("Matrix was not compressed row sparse");
        LEAVE_FUNCTION();
        return mtx_bad_param;
    }
    if (!y)
    {
        REPORT_ERROR_MESSAGE("Vector y pointer was null");
        LEAVE_FUNCTION();
        return mtx_bad_param;
    }
    if (!x)
    {
        REPORT_ERROR_MESSAGE("Vector x pointer was null");
        LEAVE_FUNCTION();
        return mtx_bad_param;
    }
#endif
    //  Length of x and y
    const uint n = mtx->columns;

    //  Initial guess of x is just y, which would be the case if matrix is just I
    memcpy(x, y, n * sizeof *x);
    //  Improve the guess by assuming that mtx is a diagonal matrix
    for (uint i = 0; i < n; ++i)
    {
        scalar_t d;
        matrix_crs_get_element(mtx, i, i, &d);
        x[i] /= d;
    }


    scalar_t err;
    uint n_iterations = 0;
    do
    {
        err = 0;

        for (uint i = 0; i < n; ++i)
        {
            scalar_t* row_ptr;
            uint* index_ptr;
            uint n_elements;
            matrix_crs_get_row(mtx, i, &n_elements, &index_ptr, &row_ptr);
            scalar_t res = (scalar_t)0.0;
            uint k = 0;
            for (uint j = 0; j < n_elements; ++j)
            {
                if (i != index_ptr[j])
                {
                    res += row_ptr[j] * x[index_ptr[j]];
                }
                else
                {
                    k = j;
                }
            }
            const scalar_t new_x = (y[i] - res) / row_ptr[k];
            x[i] = new_x;
        }

        err = 0;
        for (uint i = 0; i < n; ++i)
        {
            scalar_t val;
            matrix_crs_vector_multiply_row(mtx, x, i, &val);
            val -= y[i];
            err += val * val;
        }
        err = sqrtf(err) / (scalar_t)n;

        n_iterations += 1;
    } while(err > convergence_dif & n_iterations < n_max_iter);


    if (p_iter) *p_iter = n_iterations;
    if (p_error) *p_error = err;
    LEAVE_FUNCTION();
    return n_iterations == n_max_iter ? mtx_not_converged : mtx_success;
}


struct gauss_seidel_crs_thread_param
{
    const CrsMatrix* matrix;
    const scalar_t* y;
    scalar_t* x;
    uint* done;
    uint n;
    uint n_thrds;
    scalar_t converge_dif;
    uint n_max_iter;
    pthread_barrier_t* work_barrier;
    atomic_flag* error_flag;
    scalar_t* p_errors;
    uint rounds_done;
    uint id;
};

static void* gauss_seidel_thrd_fn(void* param)
{
    struct gauss_seidel_crs_thread_param* args = (struct gauss_seidel_crs_thread_param*)param;
    THREAD_BEGIN("Worker %s (%d/%d)", __func__, args->id, args->n_thrds);
    uint loop_type = 1;
    const uint chunk_size = args->n / args->n_thrds;
    pthread_barrier_wait(args->work_barrier);
    while (*args->done == 0)
    {
        scalar_t err;
        uint begin, end, step;
        if (loop_type)
        {
            loop_type = 0;
            begin = args->id;
            end = args->n;
            step = args->n_thrds;
        }
        else
        {
            begin = args->id * chunk_size;
            end = begin + chunk_size;
            if (args->id == args->n_thrds - 1)
            {
                end += args->n % args->n_thrds;
            }
            step = 1;
            loop_type = 1;
        }

        for (uint i = begin; i < end; i += step)
        {
            scalar_t* row_ptr;
            uint* index_ptr;
            uint n_elements;
            matrix_crs_get_row(args->matrix, i, &n_elements, &index_ptr, &row_ptr);
            scalar_t res = (scalar_t)0.0;
            uint k = 0;
            for (uint j = 0; j < n_elements; ++j)
            {
                if (i != index_ptr[j])
                {
                    res += row_ptr[j] * args->x[index_ptr[j]];
                }
                else
                {
                    k = j;
                }
            }
            const scalar_t new_x = (args->y[i] - res) / row_ptr[k];
            args->x[i] = new_x;
        }
        pthread_barrier_wait(args->work_barrier);
        err = 0;
        for (uint i = begin; i < end; ++i)
        {
            scalar_t val;
            matrix_crs_vector_multiply_row(args->matrix, args->x, i, &val);
            val -= args->y[i];
            err += val * val;
        }
        args->p_errors[args->id] = err;
        args->rounds_done += 1;
        const _Bool flag_c = atomic_flag_test_and_set(args->error_flag);
        pthread_barrier_wait(args->work_barrier);
        if (!flag_c)
        {
            err = 0;
            for (uint i = 0; i < args->n_thrds; ++i)
            {
                err += args->p_errors[i];
            }
            err = sqrtf(err) / (scalar_t)args->n;
            if (err < args->converge_dif || args->rounds_done > args->n_max_iter) *args->done = 1;
            *args->p_errors = err;
            atomic_flag_clear(args->error_flag);
        }
        pthread_barrier_wait(args->work_barrier);
    }
    THREAD_END;
    return 0;
}

mtx_res_t gauss_seidel_crs_mt(
        const CrsMatrix* mtx, const scalar_t* y, scalar_t* x, scalar_t convergence_dif, uint n_max_iter, uint* p_iter,
        scalar_t* p_error, uint n_thrds)
{
#ifdef MTX_MATRIX_CHECKS
    if (!mtx)
    {
        REPORT_ERROR_MESSAGE("Matrix pointer was null");
        LEAVE_FUNCTION();
        return mtx_bad_param;
    }
    if (mtx->type != mtx_type_crs)
    {
        REPORT_ERROR_MESSAGE("Matrix was not compressed row sparse");
        LEAVE_FUNCTION();
        return mtx_bad_param;
    }
    if (!y)
    {
        REPORT_ERROR_MESSAGE("Vector y pointer was null");
        LEAVE_FUNCTION();
        return mtx_bad_param;
    }
    if (!x)
    {
        REPORT_ERROR_MESSAGE("Vector x pointer was null");
        LEAVE_FUNCTION();
        return mtx_bad_param;
    }
#endif
    if (!n_thrds)
    {
        const mtx_res_t result = CALL_FUNCTION(gauss_seidel_crs(mtx, y, x, convergence_dif, n_max_iter, p_iter, p_error));
        LEAVE_FUNCTION();
        return result;
    }

    //  Length of x and y
    const uint n = mtx->columns;

    //  Initial guess of x is just y, which would be the case if matrix is just I
    memcpy(x, y, n * sizeof *x);
    //  Improve the guess by assuming that mtx is a diagonal matrix
    for (uint i = 0; i < n; ++i)
    {
        scalar_t d;
        matrix_crs_get_element(mtx, i, i, &d);
        x[i] /= d;
    }

    //  Memory used to store result of the current iteration
    pthread_t* const thrds = calloc(n_thrds, sizeof*thrds);
    if (!thrds)
    {
        CALLOC_FAILED(n_thrds * sizeof*thrds);
        LEAVE_FUNCTION();
        return mtx_malloc_fail;
    }

    scalar_t* const errors = calloc(n_thrds, sizeof*errors);
    if (!errors)
    {
        free(thrds);
        CALLOC_FAILED(n_thrds * sizeof*errors);
        LEAVE_FUNCTION();
        return mtx_malloc_fail;
    }
    struct gauss_seidel_crs_thread_param* const args = calloc(n_thrds, sizeof* args);
    if (!args)
    {
        free(thrds);
        free(errors);
        CALLOC_FAILED(n_thrds * sizeof*args);
        LEAVE_FUNCTION();
        return mtx_malloc_fail;
    }

    uint happy = 0;
    scalar_t max_dif = 0.0f;
    pthread_barrier_t work_barrier;
    pthread_barrier_init(&work_barrier, NULL, n_thrds);
    atomic_flag work_flag = ATOMIC_FLAG_INIT;
    struct gauss_seidel_crs_thread_param arg =
            {
                    .matrix = mtx,
                    .n = n,
                    .done = &happy,
                    .x = x,
                    .y = y,
                    .n_thrds = n_thrds,
                    .work_barrier = &work_barrier,
                    .error_flag = &work_flag,
                    .p_errors = errors,
                    .n_max_iter = n_max_iter,
                    .converge_dif = convergence_dif,
            };
    for (uint i = 0; i < n_thrds; ++i)
    {
        arg.id = i;
        args[i] = arg;
        if (pthread_create(thrds + i, NULL, gauss_seidel_thrd_fn, args + i) != 0)
        {
            happy = 1;
            REPORT_ERROR_MESSAGE("Could not create a new worker thread (%d out of %d), reason: %s", i + 1, n_thrds,
                                 strerror(errno));
            for (uint j = 0; j < i; ++j)
                pthread_cancel(thrds[j]);
            pthread_barrier_destroy(&work_barrier);
            free(args);
            free(errors);
            free(thrds);
            return mtx_thread_fail;
        }
    }


    uint n_iterations = 0;
    uint done = 0;

    for (uint i = 0; i < n_thrds; ++i)
    {
        pthread_join(thrds[i], NULL);
    }
    max_dif = errors[0];
    n_iterations = args[0].rounds_done;

    pthread_barrier_destroy(&work_barrier);

    free(errors);
    free(thrds);
    free(args);
    if (p_iter) *p_iter = n_iterations;
    if (p_error) *p_error = max_dif;
    LEAVE_FUNCTION();
    return n_iterations == n_max_iter ? mtx_not_converged : mtx_success;
}


