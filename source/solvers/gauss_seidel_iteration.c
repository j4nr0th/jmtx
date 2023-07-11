//
// Created by jan on 16.6.2022.
//

#include "gauss_seidel_iteration.h"

#include <pthread.h>
#include <stdatomic.h>
#include <stdio.h>
#include <math.h>



jmtx_result gauss_seidel_crs(
        const jmtx_matrix_crs* mtx, const jmtx_scalar_t* y, jmtx_scalar_t* x, jmtx_scalar_t convergence_dif,
        uint32_t n_max_iter, uint32_t* p_iter, jmtx_scalar_t* p_final_error, jmtx_scalar_t* p_error,
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

    //  Length of x and y
    const uint32_t n = mtx->base.cols;

    //  Initial guess of x is just y, which would be the case if matrix is just I
    memcpy(x, y, n * sizeof *x);
    //  Improve the guess by assuming that mtx is a diagonal matrix
    for (uint32_t i = 0; i < n; ++i)
    {
        jmtx_scalar_t d;
        matrix_crs_get_element(mtx, i, i, &d);
        x[i] /= d;
    }


    jmtx_scalar_t err;
    uint32_t n_iterations = 0;
    do
    {

        for (uint32_t i = 0; i < n; ++i)
        {
            jmtx_scalar_t* row_ptr;
            uint32_t* index_ptr;
            uint32_t n_elements;
            matrix_crs_get_row(mtx, i, &n_elements, &index_ptr, &row_ptr);
            jmtx_scalar_t res = (jmtx_scalar_t)0.0;
            uint32_t k = 0;
            for (uint32_t j = 0; j < n_elements; ++j)
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
            const jmtx_scalar_t new_x = (y[i] - res) / row_ptr[k];
            x[i] = new_x;
        }

        err = 0;
        for (uint32_t i = 0; i < n; ++i)
        {
            jmtx_scalar_t val;
            matrix_crs_vector_multiply_row(mtx, x, i, &val);
            val -= y[i];
            err += val * val;
        }
        err = sqrtf(err) / (jmtx_scalar_t)n;
        if (p_error)
        {
            p_error[n_iterations] = err;
        }

        n_iterations += 1;
    } while(err > convergence_dif & n_iterations < n_max_iter);


    if (p_iter) *p_iter = n_iterations;
    if (p_final_error) *p_final_error = err;
//    LEAVE_FUNCTION();
    return n_iterations == n_max_iter ? JMTX_RESULT_NOT_CONVERGED : JMTX_RESULT_SUCCESS;
}


struct gauss_seidel_crs_thread_param
{
    const jmtx_matrix_crs* matrix;
    const jmtx_scalar_t* y;
    jmtx_scalar_t* x;
    uint32_t* done;
    uint32_t n;
    uint32_t n_thrds;
    jmtx_scalar_t converge_dif;
    uint32_t n_max_iter;
    pthread_barrier_t* work_barrier;
    atomic_flag* error_flag;
    jmtx_scalar_t* p_errors;
    jmtx_scalar_t* p_err_evol;
    uint32_t rounds_done;
    uint32_t id;
};

static void* gauss_seidel_thrd_fn(void* param)
{
    struct gauss_seidel_crs_thread_param* args = (struct gauss_seidel_crs_thread_param*)param;
//    THREAD_BEGIN("Worker %s (%d/%d)", __func__, args->id, args->n_thrds);
    uint32_t loop_type = 1;
    const uint32_t chunk_size = args->n / args->n_thrds;
    pthread_barrier_wait(args->work_barrier);
    while (*args->done == 0)
    {
        jmtx_scalar_t err;
        uint32_t begin, end, step;
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

        for (uint32_t i = begin; i < end; i += step)
        {
            jmtx_scalar_t* row_ptr;
            uint32_t* index_ptr;
            uint32_t n_elements;
            matrix_crs_get_row(args->matrix, i, &n_elements, &index_ptr, &row_ptr);
            jmtx_scalar_t res = (jmtx_scalar_t)0.0;
            uint32_t k = 0;
            for (uint32_t j = 0; j < n_elements; ++j)
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
            const jmtx_scalar_t new_x = (args->y[i] - res) / row_ptr[k];
            args->x[i] = new_x;
        }
        pthread_barrier_wait(args->work_barrier);
        err = 0;
        for (uint32_t i = begin; i < end; ++i)
        {
            jmtx_scalar_t val;
            matrix_crs_vector_multiply_row(args->matrix, args->x, i, &val);
            val -= args->y[i];
            err += val * val;
        }
        args->p_errors[args->id] = err;
        args->rounds_done += 1;
        const bool flag_c = atomic_flag_test_and_set(args->error_flag);
        pthread_barrier_wait(args->work_barrier);
        if (!flag_c)
        {
            err = 0;
            for (uint32_t i = 0; i < args->n_thrds; ++i)
            {
                err += args->p_errors[i];
            }
            err = sqrtf(err) / (jmtx_scalar_t)args->n;
            if (args->p_err_evol)
            {
                args->p_err_evol[args->rounds_done] = err;
            }
            if (err < args->converge_dif || args->rounds_done > args->n_max_iter) *args->done = 1;
            *args->p_errors = err;
            atomic_flag_clear(args->error_flag);
        }
        pthread_barrier_wait(args->work_barrier);
    }

    return 0;
}

jmtx_result gauss_seidel_crs_mt(
        const jmtx_matrix_crs* mtx, const jmtx_scalar_t* y, jmtx_scalar_t* x, jmtx_scalar_t convergence_dif,
        uint32_t n_max_iter, uint32_t* p_iter, jmtx_scalar_t* p_final_error, jmtx_scalar_t* p_error,
        const jmtx_allocator_callbacks* allocator_callbacks, uint32_t n_thrds)
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

    if (n_thrds < 2)
    {
        const jmtx_result result = gauss_seidel_crs(
                        mtx, y, x, convergence_dif, n_max_iter, p_iter, p_final_error, p_error, allocator_callbacks);
        return result;
    }

    //  Length of x and y
    const uint32_t n = mtx->base.cols;

    //  Initial guess of x is just y, which would be the case if matrix is just I
    memcpy(x, y, n * sizeof *x);
    //  Improve the guess by assuming that mtx is a diagonal matrix
    for (uint32_t i = 0; i < n; ++i)
    {
        jmtx_scalar_t d;
        matrix_crs_get_element(mtx, i, i, &d);
        x[i] /= d;
    }

    //  Memory used to store result of the current iteration
    pthread_t* const thrds = allocator_callbacks->alloc(allocator_callbacks->state, n_thrds * sizeof(*thrds));
    if (!thrds)
    {
//        CALLOC_FAILED(n_thrds * sizeof*thrds);
//        LEAVE_FUNCTION();
        return JMTX_RESULT_BAD_ALLOC;
    }

    jmtx_scalar_t* const errors = allocator_callbacks->alloc(allocator_callbacks->state, n_thrds * sizeof*errors);
    if (!errors)
    {
        allocator_callbacks->free(allocator_callbacks->state, thrds);
//        CALLOC_FAILED(n_thrds * sizeof*errors);
//        LEAVE_FUNCTION();
        return JMTX_RESULT_BAD_ALLOC;
    }
    struct gauss_seidel_crs_thread_param* const args = allocator_callbacks->alloc(allocator_callbacks->state, n_thrds * sizeof* args);
    if (!args)
    {
        allocator_callbacks->free(allocator_callbacks->state, errors);
        allocator_callbacks->free(allocator_callbacks->state, thrds);
//        CALLOC_FAILED(n_thrds * sizeof*args);
//        LEAVE_FUNCTION();
        return JMTX_RESULT_BAD_ALLOC;
    }

    uint32_t happy = 0;
    jmtx_scalar_t max_dif = 0.0f;
    pthread_barrier_t work_barrier;
    pthread_barrier_init(&work_barrier, NULL, n_thrds);
    atomic_flag work_flag;
    atomic_flag_clear(&work_flag);
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
                    .p_err_evol = p_error
            };
    for (uint32_t i = 0; i < n_thrds; ++i)
    {
        arg.id = i;
        args[i] = arg;
        if (pthread_create(thrds + i, NULL, gauss_seidel_thrd_fn, args + i) != 0)
        {
            happy = 1;
//            REPORT_ERROR_MESSAGE("Could not create a new worker thread (%d out of %d), reason: %s", i + 1, n_thrds,
//                                 strerror(errno));
            for (uint32_t j = 0; j < i; ++j)
                pthread_cancel(thrds[j]);
            pthread_barrier_destroy(&work_barrier);
            allocator_callbacks->free(allocator_callbacks->state, args);
            allocator_callbacks->free(allocator_callbacks->state, errors);
            allocator_callbacks->free(allocator_callbacks->state, thrds);
            return JMTX_RESULT_BAD_THREAD;
        }
    }


    uint32_t n_iterations = 0;
    uint32_t done = 0;

    for (uint32_t i = 0; i < n_thrds; ++i)
    {
        pthread_join(thrds[i], NULL);
    }
    max_dif = errors[0];
    n_iterations = args[0].rounds_done;

    pthread_barrier_destroy(&work_barrier);

    allocator_callbacks->free(allocator_callbacks->state, args);
    allocator_callbacks->free(allocator_callbacks->state, errors);
    allocator_callbacks->free(allocator_callbacks->state, thrds);
    if (p_iter) *p_iter = n_iterations;
    if (p_final_error) *p_final_error = max_dif;
//    LEAVE_FUNCTION();
    return n_iterations == n_max_iter ? JMTX_RESULT_NOT_CONVERGED : JMTX_RESULT_SUCCESS;
}


