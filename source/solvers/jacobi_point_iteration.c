//
// Created by jan on 15.6.2022.
//

#include "jacobi_point_iteration.h"
#include <math.h>
#include <stdio.h>
#include <pthread.h>
#include <assert.h>


jmtx_result jacobi_crs(
        const jmtx_matrix_crs* mtx, const jmtx_scalar_t* y, jmtx_scalar_t* x, jmtx_scalar_t convergence_dif,
        uint32_t n_max_iter, uint32_t* p_iter, jmtx_scalar_t* p_error, jmtx_scalar_t* p_final_error,
        const jmtx_allocator_callbacks* allocator_callbacks)
{
    if (!mtx)
    {
//        REPORT_ERROR_MESSAGE("Matrix pointer was null");
//        LEAVE_FUNCTION();
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.rows != mtx->base.cols)
    {
        //  I am only doing square matrices!!!
        return JMTX_RESULT_BAD_MATRIX;
    }
    if (mtx->base.type != JMTX_TYPE_CRS)
    {
//        REPORT_ERROR_MESSAGE("Matrix was not compressed row sparse");
//        LEAVE_FUNCTION();
        return JMTX_RESULT_NULL_PARAM;
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
    else if (!allocator_callbacks->alloc || !allocator_callbacks->free)
    {
        return JMTX_RESULT_BAD_PARAM;
    }

    //  Length of x and y
    const uint32_t n = mtx->base.cols;
    jmtx_result mtx_res;
    //  Initial guess by assuming that mtx is a diagonal matrix
    for (uint32_t i = 0; i < n; ++i)
    {
        jmtx_scalar_t d;
        mtx_res = matrix_crs_get_element(mtx, i, i, &d);
        assert(mtx_res == JMTX_RESULT_SUCCESS);
        if (d == 0.0f)
        {
            //  Diagonal entry is zero!
            //  Can't solve this one with Jacobi
            return JMTX_RESULT_BAD_MATRIX;
        }
        x[i] = y[i] / d;
    }

    //  Memory used to store result of the current iteration
    jmtx_scalar_t* const auxiliary_x = allocator_callbacks->alloc(allocator_callbacks->state, n * sizeof*x);
    if (!auxiliary_x)
    {
//        CALLOC_FAILED(n * sizeof*x);
//        LEAVE_FUNCTION();
        return JMTX_RESULT_BAD_ALLOC;
    }

    jmtx_scalar_t* x0 = auxiliary_x;
    jmtx_scalar_t* x1 = x;
    jmtx_scalar_t err;
    uint32_t n_iterations = 0;
    do
    {
        err = 0.0f;
        {
            jmtx_scalar_t* tmp = x1;
            x1 = x0;
            x0 = tmp;
        }

        //  For each entry, find the corresponding row in matrix A - D and compute the dot product between x and that row
        for (uint32_t i = 0; i < n; ++i)
        {
            jmtx_scalar_t* row_ptr;
            uint32_t* index_ptr;
            uint32_t n_elements;
            matrix_crs_get_row(mtx, i, &n_elements, &index_ptr, &row_ptr);
            jmtx_scalar_t res = 0;
            uint32_t k = 0;
            for (uint32_t j = 0; j < n_elements; ++j)
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
            //  Multiplication of vector x by D⁻¹
            x1[i] = (y[i] - res) / row_ptr[k];
        }

        for (uint32_t i = 0; i < n; ++i)
        {
            jmtx_scalar_t val;
            matrix_crs_vector_multiply_row(mtx, x1, i, &val);
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


    if (x1 == auxiliary_x)
    {
        memcpy(x, auxiliary_x, sizeof*x * n);
    }
    allocator_callbacks->free(allocator_callbacks->state, auxiliary_x);
    if (p_iter) *p_iter = n_iterations;
    if (p_final_error) *p_final_error = err;
//    LEAVE_FUNCTION();
    return n_iterations == n_max_iter ? JMTX_RESULT_NOT_CONVERGED : JMTX_RESULT_SUCCESS;
}

struct jacobi_crs_thread_param
{
    const jmtx_matrix_crs* matrix;
    const jmtx_scalar_t* y;
    jmtx_scalar_t** x0, **x1;
    uint32_t* done;
    uint32_t n;
    uint32_t n_thrds;
    pthread_barrier_t* barrier;
    jmtx_scalar_t* p_err;
    jmtx_scalar_t* p_final_err;
    uint32_t id;
    _Atomic uint32_t* ready;
};

static void* jacobi_crs_thread_fn(void* param)
{
    struct jacobi_crs_thread_param args = *(struct jacobi_crs_thread_param*)param;
    *args.ready = 1;
//    THREAD_BEGIN("Worker %s (%d/%d)", __func__, args.id, args.n_thrds);
    while (!*args.done)
    {
        jmtx_scalar_t residual_sum = 0.0f;
        jmtx_scalar_t* const x0 = *args.x0;
        jmtx_scalar_t* const x1 = *args.x1;

        for (uint32_t i = args.id; i < args.n; i += args.n_thrds)
        {
            jmtx_scalar_t* row_ptr;
            uint32_t* index_ptr;
            uint32_t n_elements;
            matrix_crs_get_row(args.matrix, i, &n_elements, &index_ptr, &row_ptr);
            jmtx_scalar_t res = (jmtx_scalar_t)0.0;
            uint32_t k = 0;
            for (uint32_t j = 0; j < n_elements; ++j)
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
                jmtx_scalar_t iter_dif = (x1[i] - x0[i]);
                uint32_t f_abs_d = *(uint32_t*)(&iter_dif)& 0x7FFFFFFF;
                iter_dif = *(jmtx_scalar_t*)&f_abs_d;
                residual_sum += iter_dif;
            }
        }
        *args.p_err = residual_sum;
        pthread_barrier_wait(args.barrier);
        pthread_barrier_wait(args.barrier);
    }
//    THREAD_END;
    return 0;
}

static inline void print_vector(const jmtx_scalar_t* x, const uint32_t len)
{
    printf("\n[");
    for (uint32_t i = 0; i < len; ++i)
    {
        printf(" %g", x[i]);
    }
    printf(" ]\n");
}

jmtx_result jacobi_crs_mt(
        const jmtx_matrix_crs* mtx, const jmtx_scalar_t* y, jmtx_scalar_t* x, jmtx_scalar_t convergence_dif,
        uint32_t n_max_iter, uint32_t* p_iter, jmtx_scalar_t* p_error, jmtx_scalar_t* p_final_error,
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

    if (!n_thrds || n_thrds == 1)
    {
        return jacobi_crs(mtx, y, x, convergence_dif, n_max_iter, p_iter, p_error, p_final_error, allocator_callbacks);
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
    jmtx_scalar_t* const auxiliary_x = allocator_callbacks->alloc(allocator_callbacks->state, n * sizeof*x);
    if (!auxiliary_x)
    {
//        CALLOC_FAILED(n_thrds * sizeof*auxiliary_x);
//        LEAVE_FUNCTION();
        return JMTX_RESULT_BAD_ALLOC;
    }
    pthread_t* const thrds = allocator_callbacks->alloc(allocator_callbacks->state, n_thrds * sizeof*thrds);
    if (!thrds)
    {
        allocator_callbacks->free(allocator_callbacks->state, auxiliary_x);
//        CALLOC_FAILED(n_thrds * sizeof*thrds);
//        LEAVE_FUNCTION();
        return JMTX_RESULT_BAD_ALLOC;
    }

    jmtx_scalar_t* const errors = allocator_callbacks->alloc(allocator_callbacks->state, n_thrds * sizeof*errors);
    if (!errors)
    {
        allocator_callbacks->free(allocator_callbacks->state, auxiliary_x);
        allocator_callbacks->free(allocator_callbacks->state, thrds);
//        CALLOC_FAILED(n_thrds * sizeof*errors);
//        LEAVE_FUNCTION();
        return JMTX_RESULT_BAD_ALLOC;
    }
    pthread_barrier_t barrier;
    pthread_barrier_init(&barrier, NULL, n_thrds + 1);

    uint32_t happy = 0;
    jmtx_scalar_t* x0 = x;
    jmtx_scalar_t* x1 = auxiliary_x;
    jmtx_scalar_t err;
    _Atomic uint32_t thrd_is_ready = 0;
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
            .p_err = p_error,
            .p_final_err = p_final_error
            };

    for (uint32_t i = 0; i < n_thrds; ++i)
    {
        arg.id = i;
        arg.p_err = errors + i;
        if (pthread_create(thrds + i, NULL, jacobi_crs_thread_fn, &arg) != 0)
        {
            happy = 1;
            for (uint32_t j = 0; j < i; ++j)
                pthread_cancel(thrds[j]);
            allocator_callbacks->free(allocator_callbacks->state, errors);
            pthread_barrier_destroy(&barrier);
            allocator_callbacks->free(allocator_callbacks->state, auxiliary_x);
            allocator_callbacks->free(allocator_callbacks->state, thrds);
//            REPORT_ERROR_MESSAGE("Could not create a new worker thread (%d out of %d), reason: %s", i + 1, n_thrds,
//                                 strerror(errno));
            return JMTX_RESULT_BAD_THREAD;
        }
    }


    uint32_t n_iterations = 0;
    do
    {
        pthread_barrier_wait(&barrier);
        {
            jmtx_scalar_t* tmp = x1;
            x1 = x0;
            x0 = tmp;
        }

        err = 0;

//        print_vector(errors, n_thrds);
//        max_dif = errors[0];
        for (uint32_t i = 0; i < n_thrds; ++i)
        {
            err += errors[i];
        }
        err /= (jmtx_scalar_t)n;
        if (*p_error)
        {
            p_error[n_iterations] = err;
        }
        happy = !(err > convergence_dif & n_iterations < n_max_iter);
//        print_vector(x0, n);
//        print_vector(x1, n);
        pthread_barrier_wait(&barrier);
        n_iterations += 1;
    } while(!happy);

    for (uint32_t i = 0; i < n_thrds; ++i)
    {
        pthread_join(thrds[i], NULL);
    }
    pthread_barrier_destroy(&barrier);

    if (x0 == auxiliary_x)
    {
        memcpy(x, auxiliary_x, sizeof*x * n);
    }
    allocator_callbacks->free(allocator_callbacks->state, errors);
    allocator_callbacks->free(allocator_callbacks->state, thrds);
    allocator_callbacks->free(allocator_callbacks->state, auxiliary_x);
    if (n_iterations) *p_iter = n_iterations;
    if (p_final_error) *p_final_error = err;
//    LEAVE_FUNCTION();
    return n_iterations == n_max_iter ? JMTX_RESULT_NOT_CONVERGED : JMTX_RESULT_SUCCESS;
}

jmtx_result jacobi_relaxed_crs(
        const jmtx_matrix_crs* mtx, const jmtx_scalar_t* y, jmtx_scalar_t* x, jmtx_scalar_t relaxation_factor,
        jmtx_scalar_t convergence_dif, uint32_t n_max_iter, uint32_t* p_iter, jmtx_scalar_t* p_error,
        jmtx_scalar_t* p_final_error, const jmtx_allocator_callbacks* allocator_callbacks)
{
    if (!mtx)
    {
//        REPORT_ERROR_MESSAGE("Matrix pointer was null");
//        LEAVE_FUNCTION();
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.rows != mtx->base.cols)
    {
        //  I am only doing square matrices!!!
        return JMTX_RESULT_BAD_MATRIX;
    }
    if (mtx->base.type != JMTX_TYPE_CRS)
    {
//        REPORT_ERROR_MESSAGE("Matrix was not compressed row sparse");
//        LEAVE_FUNCTION();
        return JMTX_RESULT_NULL_PARAM;
    }
    if (relaxation_factor <= 0.0f)
    {
        //  Can't solve it if you move nowhere >:(
        return JMTX_RESULT_BAD_PARAM;
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
    else if (!allocator_callbacks->alloc || !allocator_callbacks->free)
    {
        return JMTX_RESULT_BAD_PARAM;
    }

    //  Length of x and y
    const uint32_t n = mtx->base.cols;
    jmtx_result mtx_res;

    //  Memory used to store result of the current iteration
    jmtx_scalar_t* const x_aux = allocator_callbacks->alloc(allocator_callbacks->state, n * sizeof*x);
    if (!x_aux)
    {
//        CALLOC_FAILED(n * sizeof*x);
//        LEAVE_FUNCTION();
        return JMTX_RESULT_BAD_ALLOC;
    }
    jmtx_scalar_t* const division_factors = allocator_callbacks->alloc(allocator_callbacks->state, n * sizeof*x);
    if (!division_factors)
    {
        allocator_callbacks->free(allocator_callbacks->state, x_aux);
        return JMTX_RESULT_BAD_ALLOC;
    }


    //  Initial guess by assuming that mtx is a diagonal matrix
    for (uint32_t i = 0; i < n; ++i)
    {
        jmtx_scalar_t d;
        mtx_res = matrix_crs_get_element(mtx, i, i, &d);
        assert(mtx_res == JMTX_RESULT_SUCCESS);
        if (d == 0.0f)
        {
            //  Diagonal entry is zero!
            //  Can't solve this one with Jacobi
            return JMTX_RESULT_BAD_MATRIX;
        }
        x[i] = y[i] / d;
        division_factors[i] = relaxation_factor / d;
    }



    jmtx_scalar_t* x0 = x_aux;
    jmtx_scalar_t* x1 = x;
    jmtx_scalar_t err;
    uint32_t n_iterations = 0;
    do
    {
        err = 0.0f;
        {
            //  Swap pointers
            jmtx_scalar_t* tmp = x1;
            x1 = x0;
            x0 = tmp;
        }

        //  For each entry, find the corresponding row in matrix A - D and compute the dot product between x and that row
        for (uint32_t i = 0; i < n; ++i)
        {
            jmtx_scalar_t* row_ptr;
            uint32_t* index_ptr;
            uint32_t n_elements;
            matrix_crs_get_row(mtx, i, &n_elements, &index_ptr, &row_ptr);
            jmtx_scalar_t res = 0;
            uint32_t k = 0;
            for (uint32_t j = 0; j < n_elements; ++j)
            {
                if (i != index_ptr[j])
                {
                    res += row_ptr[j] * x0[index_ptr[j]];
                }
//                else
//                {
//                    k = j;
//                }
            }
            //  Multiplication of vector x by ωD⁻¹
            x1[i] = x0[i] + (y[i] - res) * division_factors[i];
        }

        for (uint32_t i = 0; i < n; ++i)
        {
            jmtx_scalar_t val;
            matrix_crs_vector_multiply_row(mtx, x1, i, &val);
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


    if (x1 == x_aux)
    {
        memcpy(x, x_aux, sizeof*x * n);
    }
    allocator_callbacks->free(allocator_callbacks->state, division_factors);
    allocator_callbacks->free(allocator_callbacks->state, x_aux);
    if (p_iter) *p_iter = n_iterations;
    if (p_final_error) *p_final_error = err;
//    LEAVE_FUNCTION();
    return n_iterations == n_max_iter ? JMTX_RESULT_NOT_CONVERGED : JMTX_RESULT_SUCCESS;
}
