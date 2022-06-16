//
// Created by jan on 15.6.2022.
//

#include "jacobi_point_iteration.h"

#include <stdio.h>
#include <threads.h>

mtx_res_t jacobi_crs(
        const CrsMatrix* mtx, const scalar_t* y, scalar_t* x, scalar_t convergence_dif, uint n_max_iter, uint* n_iter)
{
    //  Length of x and y
    const uint n = mtx->columns;

    //  Initial guess of x is just y, which would be the case if matrix is just I
    memcpy(x, y, n * sizeof *x);

    //  Memory used to store result of the current iteration
    scalar_t* const auxiliary_x = calloc(n, sizeof*x);
    if (!auxiliary_x) return -1;

    scalar_t* x0 = auxiliary_x;
    scalar_t* x1 = x;
    scalar_t max_dif = 0.0f;
    uint n_iterations = 0;
    do
    {
        max_dif = 0.0f;
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
            scalar_t iter_dif = (x1[i] - x0[i]);
            uint f_abs_d = *(uint*)(&iter_dif)& 0x7FFFFFFF;
            iter_dif = *(scalar_t*)&f_abs_d;
            if (iter_dif > max_dif) max_dif = iter_dif;
        }

//        printf("\n[");
//        for (uint i = 0; i < 10; ++i)
//        {
//            printf(" %g", x0[i]);
//        }
//        printf(" ]\n");
//
//        printf("\n[");
//        for (uint i = 0; i < 10; ++i)
//        {
//            printf(" %g", x1[i]);
//        }
//        printf(" ]\n");


        n_iterations += 1;
    } while(max_dif > convergence_dif & n_iterations < n_max_iter);


    if (x1 == auxiliary_x)
    {
        memcpy(x, auxiliary_x, sizeof*x * n);
    }
    free(auxiliary_x);
    if (n_iterations) *n_iter = n_iterations;
    return 0;
}

struct jacobi_crs_thread_param
{
    const CrsMatrix* matrix;
    const scalar_t* y;
    scalar_t** x0, **x1;
    uint* done;
    uint n;
    uint n_thrds;
    mtx_t* mtx;
    cnd_t* cnd;
    scalar_t* p_err;
    uint id;
    _Atomic uint* ready;
};

static int jacobi_crs_thread_fn(void* param)
{
    struct jacobi_crs_thread_param args = *(struct jacobi_crs_thread_param*)param;
    *args.ready = 1;
    mtx_lock(args.mtx);
    cnd_wait(args.cnd, args.mtx);
    while (!*args.done)
    {
        scalar_t max_dif = 0.0f;
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
                if (iter_dif > max_dif) max_dif = iter_dif;
            }
        }
        *args.p_err = max_dif;
        cnd_wait(args.cnd, args.mtx);
    }
    mtx_unlock(args.mtx);
    return 0;
}

static void print_vector(const scalar_t* x, const uint len)
{
    printf("\n[");
    for (uint i = 0; i < len; ++i)
    {
        printf(" %g", x[i]);
    }
    printf(" ]\n");
}

mtx_res_t jacobi_crs_mt(
        const CrsMatrix* mtx, const scalar_t* y, scalar_t* x, scalar_t convergence_dif, uint n_max_iter, uint* n_iter,
        uint n_thrds)
{
    if (!n_thrds)
    {
        return jacobi_crs(mtx, y, x, convergence_dif, n_max_iter, n_iter);
    }

    //  Length of x and y
    const uint n = mtx->columns;

    //  Initial guess of x is just y, which would be the case if matrix is just I
    memcpy(x, y, n * sizeof *x);

    //  Memory used to store result of the current iteration
    scalar_t* const auxiliary_x = calloc(n, sizeof*x);
    if (!auxiliary_x) return -1;
    thrd_t* const thrds = calloc(n_thrds, sizeof*thrds);
    if (!thrds)
    {
        free(auxiliary_x);
        return -1;
    }
    mtx_t* const mutexes = calloc(n_thrds, sizeof*mutexes);
    if (!mutexes)
    {
        free(auxiliary_x);
        free(thrds);
        return -1;
    }
    scalar_t* const errors = calloc(n_thrds, sizeof*errors);
    if (!errors)
    {
        free(auxiliary_x);
        free(thrds);
        free(mutexes);
        return -1;
    }
    cnd_t condition;
    if (cnd_init(&condition) != thrd_success)
    {
        free(errors);
        free(mutexes);
        free(auxiliary_x);
        free(thrds);
        return -1;
    }

    for (uint i = 0; i < n_thrds; ++i)
    {
        if (mtx_init(mutexes + i, mtx_plain) != thrd_success)
        {
            cnd_destroy(&condition);
            for (uint j = 0; j < i; ++j)
            {
                mtx_destroy(mutexes + j);
            }
            free(mutexes);
            free(auxiliary_x);
            free(thrds);
            free(errors);
            return -1;
        }
    }
    uint happy = 0;
    scalar_t* x0 = auxiliary_x;
    scalar_t* x1 = x;
    scalar_t max_dif = 0.0f;
    _Atomic uint thrd_is_ready = 0;
    struct jacobi_crs_thread_param arg =
            {
            .mtx = 0,
            .cnd = &condition,
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
        arg.mtx = mutexes +i;
        arg.id = i;
        arg.p_err = errors + i;
        if (thrd_create(thrds + i, jacobi_crs_thread_fn, &arg) != thrd_success)
        {
            happy = 1;
            cnd_broadcast(&condition);
            for (uint j = 0; j < n_thrds; ++j)
            {
                if (j < i)
                {
                    thrd_join(thrds[j], NULL);
                }
                mtx_destroy(mutexes + j);
            }
            cnd_destroy(&condition);
            free(errors);
            free(mutexes);
            free(auxiliary_x);
            free(thrds);
            return -1;
        }
        while (!(thrd_is_ready))
        {
            struct timespec sleep_time = {.tv_nsec = 100};
            thrd_sleep(&sleep_time, NULL);
        }
        //  When the loop is over, the thread that was created is waiting on the condition
        thrd_is_ready = 0;
    }


    uint n_iterations = 0;
    do
    {
        {
            scalar_t* tmp = x1;
            x1 = x0;
            x0 = tmp;
        }

        cnd_broadcast(&condition);

        for (uint i = 0; i < n_thrds; ++i)
        {
            mtx_lock(mutexes + i);
        }

        for (uint i = 0; i < n_thrds; ++i)
        {
            mtx_unlock(mutexes + i);
        }

//        print_vector(errors, n_thrds);
        max_dif = errors[0];
        for (uint i = 0; i < n_thrds; ++i)
        {
            if (errors[i] > max_dif) max_dif = errors[i];
        }
        happy = !(max_dif > convergence_dif & n_iterations < n_max_iter);
//        print_vector(x0, n);
//        print_vector(x1, n);

        n_iterations += 1;
    } while(!happy);
    cnd_broadcast(&condition);
    for (uint i = 0; i < n_thrds; ++i)
    {
        thrd_join(thrds[i], NULL);
        mtx_destroy(mutexes + i);
    }
    cnd_destroy(&condition);

    if (x1 == auxiliary_x)
    {
        memcpy(x, auxiliary_x, sizeof*x * n);
    }
    free(errors);
    free(mutexes);
    free(thrds);
    free(auxiliary_x);
    if (n_iterations) *n_iter = n_iterations;
    return 0;
}
