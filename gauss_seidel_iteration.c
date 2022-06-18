//
// Created by jan on 16.6.2022.
//

#include "gauss_seidel_iteration.h"

#include <pthread.h>
#include <stdatomic.h>
#include <stdio.h>

//#include <unistd.h>
//#include <errno.h>


mtx_res_t gauss_seidel_crs(
        const CrsMatrix* mtx, const scalar_t* y, scalar_t* x, scalar_t convergence_dif, uint n_max_iter, uint* p_iter,
        scalar_t* p_error)
{
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


    scalar_t max_dif;
    uint n_iterations = 0;
    do
    {
        max_dif = 0.0f;

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
            {
                scalar_t iter_dif = (x[i] - new_x);
                uint f_abs_d = *(uint*)(&iter_dif)& 0x7FFFFFFF;
                iter_dif = *(scalar_t*)&f_abs_d;
                max_dif += iter_dif;
            }

            x[i] = new_x;
        }
        max_dif /= (scalar_t)n;

        n_iterations += 1;
    } while(max_dif > convergence_dif & n_iterations < n_max_iter);


    if (p_iter) *p_iter = n_iterations;
    if (p_error) *p_error = max_dif;
    return 0;
}

//struct gauss_seidel_crs_thread_param
//{
//    const CrsMatrix* matrix;
//    const scalar_t* y;
//    scalar_t* x;
//    uint* done;
//    uint n;
//    uint n_thrds;
//    pthread_barrier_t* barrier;
//    scalar_t* p_err;
//    uint rounds_done;
//    uint id;
//};
//
//static void* gauss_seidel_thrd_fn(void* param)
//{
//    struct gauss_seidel_crs_thread_param* args = (struct gauss_seidel_crs_thread_param*)param;
//    uint loop_type = 1;
//    const uint chunk_size = args->n / args->n_thrds;
//    while (*args->done == 0)
//    {
//        scalar_t max_dif = 0.0f;
//        uint begin, end, step;
//        if (loop_type)
//        {
//            loop_type = 0;
//            begin = args->id;
//            end = args->n;
//            step = args->n_thrds;
//        }
//        else
//        {
//            begin = args->id * chunk_size;
//            end = begin + chunk_size;
//            if (args->id == args->n_thrds - 1)
//            {
//                end += args->n % args->n_thrds;
//            }
//            step = 1;
//            loop_type = 1;
//        }
//
//        for (uint i = begin; i < end; i += step)
//        {
//            scalar_t* row_ptr;
//            uint* index_ptr;
//            uint n_elements;
//            matrix_crs_get_row(args->matrix, i, &n_elements, &index_ptr, &row_ptr);
//            scalar_t res = (scalar_t)0.0;
//            uint k = 0;
//            for (uint j = 0; j < n_elements; ++j)
//            {
//                if (i != index_ptr[j])
//                {
//                    res += row_ptr[j] * args->x[index_ptr[j]];
//                }
//                else
//                {
//                    k = j;
//                }
//            }
//            const scalar_t new_x = (args->y[i] - res) / row_ptr[k];
//            {
//                scalar_t iter_dif = (args->x[i] - new_x);
//                uint f_abs_d = *(uint*)(&iter_dif)& 0x7FFFFFFF;
//                iter_dif = *(scalar_t*)&f_abs_d;
//                max_dif += iter_dif;
//            }
//            args->x[i] = new_x;
//        }
//        *args->p_err = max_dif;
//        args->rounds_done += 1;
//        pthread_barrier_wait(args->barrier);
//        pthread_barrier_wait(args->barrier);
//
//    }
//    return 0;
//}
//
//mtx_res_t gauss_seidel_crs_mt(
//        const CrsMatrix* mtx, const scalar_t* y, scalar_t* x, scalar_t convergence_dif, uint n_max_iter, uint* n_iter,
//        uint n_thrds)
//{
//    if (!n_thrds)
//    {
//        return gauss_seidel_crs(mtx, y, x, convergence_dif, n_max_iter, n_iter);
//    }
//
//    //  Length of x and y
//    const uint n = mtx->columns;
//
//    //  Initial guess of x is just y, which would be the case if matrix is just I
//    memcpy(x, y, n * sizeof *x);
//    //  Improve the guess by assuming that mtx is a diagonal matrix
//    for (uint i = 0; i < n; ++i)
//    {
//        scalar_t d;
//        matrix_crs_get_element(mtx, i, i, &d);
//        x[i] /= d;
//    }
//
//    //  Memory used to store result of the current iteration
//    pthread_t* const thrds = calloc(n_thrds, sizeof*thrds);
//    if (!thrds)
//    {
//        return -1;
//    }
//
//    scalar_t* const errors = calloc(n_thrds, sizeof*errors);
//    if (!errors)
//    {
//        free(thrds);
//        return -1;
//    }
//    struct gauss_seidel_crs_thread_param* const args = calloc(n_thrds, sizeof* args);
//    if (!args)
//    {
//        free(thrds);
//        free(errors);
//        return -1;
//    }
//
////    if (pthread_cond_init(&condition) != thrd_success)
////    {
////        free(args);
////        free(errors);
////        free(mutexes);
////        free(thrds);
////        return -1;
////    }
//
////    for (uint i = 0; i < n_thrds; ++i)
////    {
////        if (mtx_init(mutexes + i, mtx_plain) != thrd_success)
////        {
////            cnd_destroy(&condition);
////            for (uint j = 0; j < i; ++j)
////            {
////                mtx_destroy(mutexes + j);
////            }
////            free(args);
////            free(mutexes);
////            free(thrds);
////            free(errors);
////            return -1;
////        }
////    }
//
//    uint happy = 0;
//    scalar_t max_dif = 0.0f;
//    pthread_barrier_t barrier;
//    pthread_barrier_init(&barrier, NULL, n_thrds + 1);
//    struct gauss_seidel_crs_thread_param arg =
//            {
//                    .matrix = mtx,
//                    .n = n,
//                    .done = &happy,
//                    .x = x,
//                    .y = y,
//                    .n_thrds = n_thrds,
//                    .barrier = &barrier,
//            };
//    for (uint i = 0; i < n_thrds; ++i)
//    {
//        arg.id = i;
//        arg.p_err = errors + i;
//        args[i] = arg;
//        if (pthread_create(thrds + i, NULL, gauss_seidel_thrd_fn, args + i) != 0)
//        {
//            happy = 1;
////            cnd_broadcast(&condition);
////            for (uint j = 0; j < n_thrds; ++j)
////            {
////                if (j < i)
////                {
////                    thrd_join(thrds[j], NULL);
////                }
////                mtx_destroy(mutexes + j);
////            }
//            for (uint j = 0; j < i; ++j)
//                pthread_cancel(thrds[j]);
//            pthread_barrier_destroy(&barrier);
//            free(args);
//            free(errors);
//            free(thrds);
//            return -1;
//        }
//    }
//
//
//    uint n_iterations = 0;
//    uint done = 0;
//    do
//    {
//        pthread_barrier_wait(&barrier);
//
//        // print_vector(errors, n_thrds);
//        max_dif = 0;
//        for (uint i = 0; i < n_thrds; ++i)
//        {
//            max_dif += errors[i];
//        }
//        max_dif /= (scalar_t)n;
//        done = !(max_dif > convergence_dif & n_iterations < n_max_iter);
//        // print_vector(x0, n);
//        // print_vector(x1, n);
//        happy = done;
//        pthread_barrier_wait(&barrier);
//        n_iterations += 1;
//    } while(!done);
//    happy = 1;
////    pthread_barrier_wait(&barrier);
//    for (uint i = 0; i < n_thrds; ++i)
//    {
//        pthread_join(thrds[i], NULL);
//    }
//    pthread_barrier_destroy(&barrier);
//
//    free(errors);
//    free(thrds);
////    if (n_iterations == 1)
////    {
////        printf("Total residual after the single iteration: %f\n", max_dif);
////        printf("Thread work done:\n");
////        for (uint i = 0; i < n_thrds; ++i)
////        {
////            printf("\tthrd %u - %u\n", i, args[i].rounds_done);
////        }
////    }
//    free(args);
//    printf("Max_dif %g\n", max_dif);
//    if (n_iter) *n_iter = n_iterations;
//    return 0;
//}


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
    uint loop_type = 1;
    const uint chunk_size = args->n / args->n_thrds;
    pthread_barrier_wait(args->work_barrier);
    while (*args->done == 0)
    {
        scalar_t max_dif = 0.0f;
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
            {
                scalar_t iter_dif = (args->x[i] - new_x);
                uint f_abs_d = *(uint*)(&iter_dif)& 0x7FFFFFFF;
                iter_dif = *(scalar_t*)&f_abs_d;
                max_dif += iter_dif;
            }
            args->x[i] = new_x;
        }
        args->p_errors[args->id] = max_dif;
        args->rounds_done += 1;
        const _Bool flag_c = atomic_flag_test_and_set(args->error_flag);
        pthread_barrier_wait(args->work_barrier);
        if (!flag_c)
        {
            max_dif = 0;
            for (uint i = 0; i < args->n_thrds; ++i)
            {
                max_dif += args->p_errors[i];
            }
            max_dif /= (scalar_t)args->n;
            if (max_dif < args->converge_dif || args->rounds_done > args->n_max_iter) *args->done = 1;
            *args->p_errors = max_dif;
            atomic_flag_clear(args->error_flag);
        }
        pthread_barrier_wait(args->work_barrier);
    }
    return 0;
}

mtx_res_t gauss_seidel_crs_mt(
        const CrsMatrix* mtx, const scalar_t* y, scalar_t* x, scalar_t convergence_dif, uint n_max_iter, uint* p_iter,
        scalar_t* p_error, uint n_thrds)
{
    if (!n_thrds)
    {
        return gauss_seidel_crs(mtx, y, x, convergence_dif, n_max_iter, p_iter, NULL);
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
        return -1;
    }

    scalar_t* const errors = calloc(n_thrds, sizeof*errors);
    if (!errors)
    {
        free(thrds);
        return -1;
    }
    struct gauss_seidel_crs_thread_param* const args = calloc(n_thrds, sizeof* args);
    if (!args)
    {
        free(thrds);
        free(errors);
        return -1;
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
            for (uint j = 0; j < i; ++j)
                pthread_cancel(thrds[j]);
            pthread_barrier_destroy(&work_barrier);
            free(args);
            free(errors);
            free(thrds);
            return -1;
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
    return 0;
}


