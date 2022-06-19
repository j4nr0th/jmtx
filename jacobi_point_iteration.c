//
// Created by jan on 15.6.2022.
//

#include "jacobi_point_iteration.h"
#include <math.h>
#include <stdio.h>
#include <pthread.h>

#undef jacobi_crs
#undef jacobi_crs_mt

mtx_res_t jacobi_crs(
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

    //  Memory used to store result of the current iteration
    scalar_t* const auxiliary_x = calloc(n, sizeof*x);
    if (!auxiliary_x)
    {
        CALLOC_FAILED(n * sizeof*x);
        LEAVE_FUNCTION();
        return mtx_malloc_fail;
    }

    scalar_t* x0 = auxiliary_x;
    scalar_t* x1 = x;
    scalar_t err = 0.0f;
    uint n_iterations = 0;
    do
    {
        err = 0.0f;
        {
            scalar_t* tmp = x1;
            x1 = x0;
            x0 = tmp;
        }

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
                    res += row_ptr[j] * x0[index_ptr[j]];
                }
                else
                {
                    k = j;
                }
            }
            x1[i] = (y[i] - res) / row_ptr[k];
        }

        for (uint i = 0; i < n; ++i)
        {
            scalar_t val;
            matrix_crs_vector_multiply_row(mtx, x1, i, &val);
            val -= y[i];
            err += val * val;
        }
        err = sqrtf(err) / (scalar_t)n;
        n_iterations += 1;
    } while(err > convergence_dif & n_iterations < n_max_iter);


    if (x1 == auxiliary_x)
    {
        memcpy(x, auxiliary_x, sizeof*x * n);
    }
    free(auxiliary_x);
    if (p_iter) *p_iter = n_iterations;
    if (p_error) *p_error = err;
    LEAVE_FUNCTION();
    return n_iterations == n_max_iter ? mtx_not_converged : mtx_success;
}

struct jacobi_crs_thread_param
{
    const CrsMatrix* matrix;
    const scalar_t* y;
    scalar_t** x0, **x1;
    uint* done;
    uint n;
    uint n_thrds;
    pthread_barrier_t* barrier;
    scalar_t* p_err;
    uint id;
    _Atomic uint* ready;
};

static void* jacobi_crs_thread_fn(void* param)
{
    struct jacobi_crs_thread_param args = *(struct jacobi_crs_thread_param*)param;
    *args.ready = 1;
    THREAD_BEGIN("Worker %s (%d/%d)", __func__, args.id, args.n_thrds);
    while (!*args.done)
    {
        scalar_t residual_sum = 0.0f;
        scalar_t* const x0 = *args.x0;
        scalar_t* const x1 = *args.x1;

        for (uint i = args.id; i < args.n; i += args.n_thrds)
        {
            scalar_t* row_ptr;
            uint* index_ptr;
            uint n_elements;
            matrix_crs_get_row(args.matrix, i, &n_elements, &index_ptr, &row_ptr);
            scalar_t res = (scalar_t)0.0;
            uint k = 0;
            for (uint j = 0; j < n_elements; ++j)
            {
                if (i != index_ptr[j])
                {
                    res += row_ptr[j] * x0[index_ptr[j]];
                }
                else
                {
                    k = j;
                }
            }
            x1[i] = (args.y[i] - res) / row_ptr[k];
            {
                scalar_t iter_dif = (x1[i] - x0[i]);
                uint f_abs_d = *(uint*)(&iter_dif)& 0x7FFFFFFF;
                iter_dif = *(scalar_t*)&f_abs_d;
                residual_sum += iter_dif;
            }
        }
        *args.p_err = residual_sum;
        pthread_barrier_wait(args.barrier);
        pthread_barrier_wait(args.barrier);
    }
    THREAD_END;
    return 0;
}

static inline void print_vector(const scalar_t* x, const uint len)
{
    printf("\n[");
    for (uint i = 0; i < len; ++i)
    {
        printf(" %g", x[i]);
    }
    printf(" ]\n");
}

mtx_res_t jacobi_crs_mt(
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
    if (!n_thrds || n_thrds == 1)
    {
        return CALL_FUNCTION(jacobi_crs(mtx, y, x, convergence_dif, n_max_iter, p_iter, NULL));
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
    scalar_t* const auxiliary_x = calloc(n, sizeof*x);
    if (!auxiliary_x)
    {
        CALLOC_FAILED(n_thrds * sizeof*auxiliary_x);
        LEAVE_FUNCTION();
        return mtx_malloc_fail;
    }
    pthread_t* const thrds = calloc(n_thrds, sizeof*thrds);
    if (!thrds)
    {
        free(auxiliary_x);
        CALLOC_FAILED(n_thrds * sizeof*thrds);
        LEAVE_FUNCTION();
        return mtx_malloc_fail;
    }

    scalar_t* const errors = calloc(n_thrds, sizeof*errors);
    if (!errors)
    {
        free(auxiliary_x);
        free(thrds);
        CALLOC_FAILED(n_thrds * sizeof*errors);
        LEAVE_FUNCTION();
        return mtx_malloc_fail;
    }
    pthread_barrier_t barrier;
    pthread_barrier_init(&barrier, NULL, n_thrds + 1);

    uint happy = 0;
    scalar_t* x0 = x;
    scalar_t* x1 = auxiliary_x;
    scalar_t err;
    _Atomic uint thrd_is_ready = 0;
    struct jacobi_crs_thread_param arg =
            {
            .barrier = &barrier,
            .matrix = mtx,
            .n = n,
            .done = &happy,
            .x0 = &x0,
            .x1 = &x1,
            .y = y,
            .n_thrds = n_thrds,
            .ready = &thrd_is_ready,
            };

    for (uint i = 0; i < n_thrds; ++i)
    {
        arg.id = i;
        arg.p_err = errors + i;
        if (pthread_create(thrds + i, NULL, jacobi_crs_thread_fn, &arg) != 0)
        {
            happy = 1;
            for (uint j = 0; j < i; ++j)
                pthread_cancel(thrds[j]);
            free(errors);
            pthread_barrier_destroy(&barrier);
            free(auxiliary_x);
            free(thrds);
            REPORT_ERROR_MESSAGE("Could not create a new worker thread (%d out of %d), reason: %s", i + 1, n_thrds,
                                 strerror(errno));
            return mtx_thread_fail;
        }
    }


    uint n_iterations = 0;
    do
    {
        pthread_barrier_wait(&barrier);
        {
            scalar_t* tmp = x1;
            x1 = x0;
            x0 = tmp;
        }

        err = 0;

//        print_vector(errors, n_thrds);
//        max_dif = errors[0];
        for (uint i = 0; i < n_thrds; ++i)
        {
            err += errors[i];
        }
        err /= (scalar_t)n;
        happy = !(err > convergence_dif & n_iterations < n_max_iter);
//        print_vector(x0, n);
//        print_vector(x1, n);
        pthread_barrier_wait(&barrier);
        n_iterations += 1;
    } while(!happy);

    for (uint i = 0; i < n_thrds; ++i)
    {
        pthread_join(thrds[i], NULL);
    }
    pthread_barrier_destroy(&barrier);

    if (x0 == auxiliary_x)
    {
        memcpy(x, auxiliary_x, sizeof*x * n);
    }
    free(errors);
    free(thrds);
    free(auxiliary_x);
    if (n_iterations) *p_iter = n_iterations;
    if (p_error) *p_error = err;
    LEAVE_FUNCTION();
    return n_iterations == n_max_iter ? mtx_not_converged : mtx_success;
}
