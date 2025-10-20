#include "../matrices/band_row_major.h"
#include "../matrices/dense_row_major.h"
#include "../matrices/sparse_row_compressed.h"
#include "lu_solving.h"
#include "../decompositions/incomplete_lu_decomposition.h"
#include "../matrices/sparse_conversion.h"
#include <assert.h>
#include <math.h>

/**
 * Solves a problem L U x = y, where L is a lower triangular matrix with the diagonal equal to 1 and U is an upper
 * triangular matrix.
 * @param l lower triangular matrix with all entries on its main diagonal equal to 1
 * @param u upper triangular matrix
 * @param y memory containing forcing vector
 * @param x memory which receives the solution
 */
void JMTX_NAME_TYPED(solve_direct_lu_brm)(const JMTX_NAME_TYPED(matrix_brm) * l, const JMTX_NAME_TYPED(matrix_brm) * u,
                                          const JMTX_SCALAR_T *restrict y, JMTX_SCALAR_T *restrict x)
{
    const JMTX_FAST_INT_T n = l->base.cols;
    x[0] = y[0];
    //  First is the forward substitution for L v = y
    for (JMTX_FAST_INT_T i = 1; i < n; ++i)
    {
        const JMTX_FAST_INT_T off_j = JMTX_NAME_TYPED(matrix_brm_first_pos_in_row)(l, i);
        JMTX_SCALAR_T *values = NULL;
        const JMTX_FAST_INT_T len = JMTX_NAME_TYPED(matrix_brm_get_row)(l, i, &values);
        JMTX_SCALAR_T v = 0;
        for (JMTX_FAST_INT_T j = 0; j < len - 1; ++j)
        {
            v += values[j] * x[off_j + j];
        }
        x[i] = y[i] - v;
    }
    //  Then the backward substitution for U x = v
    for (int32_t i = (int32_t)n - 1; i >= 0; --i)
    {
        const JMTX_FAST_INT_T off_j = JMTX_NAME_TYPED(matrix_brm_first_pos_in_row)(u, i);
        JMTX_SCALAR_T *values;
        const JMTX_FAST_INT_T len = JMTX_NAME_TYPED(matrix_brm_get_row)(u, i, &values);
        JMTX_SCALAR_T v = 0;
        for (JMTX_INDEX_T j = 1; j < len; ++j)
        {
            v += values[j] * x[off_j + j];
        }
        x[i] = (x[i] - v) / values[0];
    }
}

/**
 * Solves a problem L U x = y, where L is a lower triangular matrix with the diagonal equal to 1 and U is an upper
 * triangular matrix. This version of the function stores the solution vector x back into the same memory where the
 * forcing vector was.
 * @param l lower triangular matrix with all entries on its main diagonal equal to 1
 * @param u upper triangular matrix
 * @param x memory which contains the forcing vector and receives the solution
 */
void JMTX_NAME_TYPED(solve_direct_lu_brm_inplace)(const JMTX_NAME_TYPED(matrix_brm) * l,
                                                  const JMTX_NAME_TYPED(matrix_brm) * u, JMTX_SCALAR_T *restrict x)
{
    const JMTX_FAST_INT_T n = l->base.cols;
    //    x[0] = x[0];
    //  First is the forward substitution for L v = y
    for (JMTX_FAST_INT_T i = 1; i < n; ++i)
    {
        const JMTX_FAST_INT_T off_j = JMTX_NAME_TYPED(matrix_brm_first_pos_in_row)(l, i);
        JMTX_SCALAR_T *values = NULL;
        const JMTX_FAST_INT_T len = JMTX_NAME_TYPED(matrix_brm_get_row)(l, i, &values);
        JMTX_SCALAR_T v = 0;
        for (JMTX_FAST_INT_T j = 0; j < len - 1; ++j)
        {
            v += values[j] * x[off_j + j];
        }
        x[i] = x[i] - v;
    }
    //  Then the backward substitution for U x = v
    for (int32_t i = (int32_t)n - 1; i >= 0; --i)
    {
        const JMTX_FAST_INT_T off_j = JMTX_NAME_TYPED(matrix_brm_first_pos_in_row)(u, i);
        JMTX_SCALAR_T *values;
        const JMTX_FAST_INT_T len = JMTX_NAME_TYPED(matrix_brm_get_row)(u, i, &values);
        JMTX_SCALAR_T v = 0;
        for (JMTX_INDEX_T j = 1; j < len; ++j)
        {
            v += values[j] * x[off_j + j];
        }
        x[i] = (x[i] - v) / values[0];
    }
}

/**
 * Solves the A x = L U x = y problem by computing the residual, then solving for L U e = r for the error e if residual
 * is too large to further refine the solution. This can help eliminate rounding errors, or alternatively can be used
 * with an incomplete LU decomposition (or ILU) to work as an iterative solver. Must be given the matrices LU.
 * @param a system matrix A
 * @param l lower triangular matrix with all entries on its main diagonal equal to 1
 * @param u upper triangular matrix
 * @param y memory containing forcing vector
 * @param x memory which receives the solution
 * @param aux_vec auxiliary memory for a vector of the same size as x and y
 * @param args::in_convergence_criterion tolerance to determine if the solution is close enough
 * @param args::in_max_iterations number of iterations to stop at
 * @param args::out_last_error receives the value of the error criterion at the final iteration
 * @param args::out_last_iteration receives the number of the final iteration
 * @param args::opt_error_evolution (optional) pointer to an array of length max_iterations, that receives the error
 * value of each iteration
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_NOT_CONVERGED if it hasn't reached given stopping criterion,
 * in case of failure it returns the associated error code
 */
jmtx_result JMTX_NAME_TYPED(solve_iterative_lu_brm_refine)(
    const JMTX_NAME_TYPED(matrix_brm) * a, const JMTX_NAME_TYPED(matrix_brm) * l, const JMTX_NAME_TYPED(matrix_brm) * u,
    const JMTX_SCALAR_T y[JMTX_ARRAY_ATTRIB(restrict)], JMTX_SCALAR_T x[JMTX_ARRAY_ATTRIB(restrict)],
    JMTX_SCALAR_T aux_vec[JMTX_ARRAY_ATTRIB(restrict)], JMTX_NAME_TYPED(solver_arguments) * args)
{
    const JMTX_FAST_INT_T n = l->base.cols;
    JMTX_REAL_T y_magnitude2 = 0;
    JMTX_FAST_INT_T iteration_count = 0;
    JMTX_REAL_T error = 0;
    for (JMTX_FAST_INT_T i = 0; i < n; ++i)
    {
        y_magnitude2 += JMTX_DOT(y[i], y[i]);
    }
    JMTX_NAME_TYPED(solve_direct_lu_brm)(l, u, y, x);

    for (;;)
    {
        for (JMTX_FAST_INT_T i = 0; i < n; ++i)
        {
            aux_vec[i] = JMTX_NAME_TYPED(matrix_brm_vector_multiply_row)(a, x, i);
        }

        for (JMTX_FAST_INT_T i = 0; i < n; ++i)
        {
            aux_vec[i] = y[i] - aux_vec[i];
        }

        JMTX_REAL_T residual_magnitude2 = 0;
        for (JMTX_FAST_INT_T i = 0; i < n; ++i)
        {
            residual_magnitude2 += JMTX_DOT(aux_vec[i], aux_vec[i]);
        }
        error = JMTX_REAL_ROOT(residual_magnitude2 / y_magnitude2);
        if (args->opt_error_evolution)
        {
            args->opt_error_evolution[iteration_count] = error;
        }
        ++iteration_count;
        if (error < args->in_convergence_criterion || iteration_count > args->in_max_iterations)
        {
            break;
        }
        JMTX_NAME_TYPED(solve_direct_lu_brm_inplace)(l, u, aux_vec);
        for (JMTX_FAST_INT_T i = 0; i < n; ++i)
        {
            x[i] += aux_vec[i];
        }
    }
    args->out_last_error = error;
    args->out_last_iteration = iteration_count;
    return iteration_count < args->in_max_iterations ? JMTX_RESULT_SUCCESS : JMTX_RESULT_NOT_CONVERGED;
}

/**
 * Solves the A x = L U x = y problem by computing the residual, then solving for L U e = r for the error e if residual
 * is too large to further refine the solution. This can help eliminate rounding errors, or alternatively can be used
 * with an incomplete LU decomposition (or ILU) to work as an iterative solver. Must be given the matrices LU.
 * This version offers a minor degree of parallelization, where calculation of the residual and applying error
 * correction are done in parallel, but the main bottleneck of dealing with inverting L U e = r is done in series.
 * @param a system matrix A
 * @param l lower triangular matrix with all entries on its main diagonal equal to 1
 * @param u upper triangular matrix
 * @param y memory containing forcing vector
 * @param x memory which receives the solution
 * @param aux_vec auxiliary memory for a vector of the same size as x and y
 * @param args::in_convergence_criterion tolerance to determine if the solution is close enough
 * @param args::in_max_iterations number of iterations to stop at
 * @param args::out_last_error receives the value of the error criterion at the final iteration
 * @param args::out_last_iteration receives the number of the final iteration
 * @param args::opt_error_evolution (optional) pointer to an array of length max_iterations, that receives the error
 * value of each iteration
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_NOT_CONVERGED if it hasn't reached given stopping criterion,
 * in case of failure it returns the associated error code
 */
jmtx_result JMTX_NAME_TYPED(solve_iterative_lu_brm_refine_parallel)(
    const JMTX_NAME_TYPED(matrix_brm) * a, const JMTX_NAME_TYPED(matrix_brm) * l, const JMTX_NAME_TYPED(matrix_brm) * u,
    const JMTX_SCALAR_T y[JMTX_ARRAY_ATTRIB(const restrict)], JMTX_SCALAR_T x[JMTX_ARRAY_ATTRIB(const restrict)],
    JMTX_SCALAR_T aux_vec[JMTX_ARRAY_ATTRIB(const restrict)], JMTX_NAME_TYPED(solver_arguments) * args)
{
    const JMTX_FAST_INT_T n = l->base.cols;
    JMTX_REAL_T y_magnitude2 = 0;
    JMTX_FAST_INT_T iteration_count = 0;
    JMTX_REAL_T error = 0;
    for (JMTX_FAST_INT_T i = 0; i < n; ++i)
    {
        y_magnitude2 += JMTX_DOT(y[i], y[i]);
    }
    JMTX_REAL_T residual_magnitude2 = 0;
    JMTX_NAME_TYPED(solve_direct_lu_brm)(l, u, y, x);
#pragma omp parallel default(none)                                                                                     \
    shared(y, x, l, u, aux_vec, args, n, error, y_magnitude2, iteration_count, residual_magnitude2, a)
    {
        for (;;)
        {
#pragma omp for
            for (JMTX_FAST_INT_T i = 0; i < n; ++i)
            {
                aux_vec[i] = JMTX_NAME_TYPED(matrix_brm_vector_multiply_row)(a, x, i);
            }

#pragma omp for
            for (JMTX_FAST_INT_T i = 0; i < n; ++i)
            {
                aux_vec[i] = y[i] - aux_vec[i];
            }

#pragma omp for reduction(+ : residual_magnitude2)
            for (JMTX_FAST_INT_T i = 0; i < n; ++i)
            {
                residual_magnitude2 += JMTX_DOT(aux_vec[i], aux_vec[i]);
            }

#pragma omp single
            {
                error = JMTX_REAL_ROOT(residual_magnitude2 / y_magnitude2);
                if (args->opt_error_evolution)
                {
                    args->opt_error_evolution[iteration_count] = error;
                }
                ++iteration_count;
            }

            if (error < args->in_convergence_criterion || iteration_count > args->in_max_iterations)
            {
                break;
            }
#pragma omp single
            {
                JMTX_NAME_TYPED(solve_direct_lu_brm_inplace)(l, u, aux_vec);
            }
#pragma omp for
            for (JMTX_FAST_INT_T i = 0; i < n; ++i)
            {
                x[i] += aux_vec[i];
            }
        }
    }
    args->out_last_error = error;
    args->out_last_iteration = iteration_count;
    return iteration_count < args->in_max_iterations ? JMTX_RESULT_SUCCESS : JMTX_RESULT_NOT_CONVERGED;
}

/**
 * Solves a problem L U x = y, where L is a lower triangular matrix with the diagonal equal to 1 and U is an upper
 * triangular matrix.
 * @param l lower triangular matrix with all entries on its main diagonal equal to 1
 * @param u upper triangular matrix
 * @param y memory containing forcing vector
 * @param x memory which receives the solution
 */
void JMTX_NAME_TYPED(solve_direct_lu_crs)(const JMTX_NAME_TYPED(matrix_crs) * l, const JMTX_NAME_TYPED(matrix_crs) * u,
                                          const JMTX_SCALAR_T *restrict y, JMTX_SCALAR_T *restrict x)
{
    const JMTX_INDEX_T n = l->base.cols;
    x[0] = y[0];
    //  First is the forward substitution for L v = y
    for (JMTX_INDEX_T i = 1; i < n; ++i)
    {
        JMTX_INDEX_T *indices;
        JMTX_SCALAR_T *values;
        JMTX_INDEX_T count = JMTX_NAME_TYPED(matrix_crs_get_row)(l, i, &indices, &values);
        assert(indices[count - 1] == (JMTX_INDEX_T)i);
        assert(values[count - 1] == 1.0f);

        JMTX_SCALAR_T v = 0;
        for (JMTX_INDEX_T j = 0; j < count - 1; ++j)
        {
            assert(indices[j] < i);
            v += values[j] * x[indices[j]];
        }
        x[i] = y[i] - v;
    }
    //  Then the backward substitution for U x = v
    for (int32_t i = (int32_t)n - 1; i >= 0; --i)
    {
        JMTX_INDEX_T *indices;
        JMTX_SCALAR_T *values;
        JMTX_INDEX_T count = JMTX_NAME_TYPED(matrix_crs_get_row)(u, i, &indices, &values);
        assert(indices[0] == (JMTX_INDEX_T)i);

        JMTX_SCALAR_T v = 0;
        for (JMTX_INDEX_T j = 1; j < count; ++j)
        {
            assert(indices[j] > (JMTX_INDEX_T)i);
            v += values[j] * x[indices[j]];
        }
        x[i] = (x[i] - v) / values[0];
    }
}

/**
 * Solves a problem L U x = y, where L is a lower triangular matrix with the diagonal equal to 1 and U is an upper
 * triangular matrix. This version of the function stores the solution vector x back into the same memory where the
 * forcing vector was.
 * @param l lower triangular matrix with all entries on its main diagonal equal to 1
 * @param u upper triangular matrix
 * @param x memory which contains the forcing vector and receives the solution
 */
void JMTX_NAME_TYPED(solve_direct_lu_crs_inplace)(const JMTX_NAME_TYPED(matrix_crs) * l,
                                                  const JMTX_NAME_TYPED(matrix_crs) * u, JMTX_SCALAR_T *restrict x)
{
    const JMTX_INDEX_T n = l->base.cols;
    //  First is the forward substitution for L v = y
    for (JMTX_INDEX_T i = 1; i < n; ++i)
    {
        JMTX_INDEX_T *indices;
        JMTX_SCALAR_T *values;
        JMTX_INDEX_T count = JMTX_NAME_TYPED(matrix_crs_get_row)(l, i, &indices, &values);
        assert(indices[count - 1] == (JMTX_INDEX_T)i);

        JMTX_SCALAR_T v = 0;
        for (JMTX_INDEX_T j = 0; j < count - 1; ++j)
        {
            assert(indices[j] < i);
            v += values[j] * x[indices[j]];
        }
        x[i] = x[i] - v;
    }
    //  Then the backward substitution for U x = v
    for (int32_t i = (int32_t)n - 1; i >= 0; --i)
    {
        JMTX_INDEX_T *indices;
        JMTX_SCALAR_T *values;
        JMTX_INDEX_T count = JMTX_NAME_TYPED(matrix_crs_get_row)(u, i, &indices, &values);
        assert(indices[0] == (JMTX_INDEX_T)i);

        JMTX_SCALAR_T v = 0;
        for (JMTX_INDEX_T j = 1; j < count; ++j)
        {
            assert(indices[j] > (JMTX_INDEX_T)i);
            v += values[j] * x[indices[j]];
        }
        x[i] = (x[i] - v) / values[0];
    }
}

static inline void compute_residual(const JMTX_INDEX_T n, const JMTX_NAME_TYPED(matrix_crs) * mtx,
                                    const JMTX_SCALAR_T *restrict x, const JMTX_SCALAR_T *restrict y,
                                    JMTX_SCALAR_T *restrict r)
{
    for (JMTX_INDEX_T i = 0; i < n; ++i)
    {
        r[i] = y[i] - JMTX_NAME_TYPED(matrix_crs_vector_multiply_row)(mtx, x, i);
    }
}

/**
 * Solves the A x = L U x = y problem by computing the residual, then solving for L U e = r for the error e if residual
 * is too large to further refine the solution. This can help eliminate rounding errors, or alternatively can be used
 * with an incomplete LU decomposition (or ILU) to work as an iterative solver. Must be given the matrices LU.
 * @param a system matrix A
 * @param y memory containing forcing vector
 * @param x memory which receives the solution
 * @param aux_vec auxiliary memory for a vector of the same size as x and y
 * @param args::in_convergence_criterion tolerance to determine if the solution is close enough
 * @param args::in_max_iterations number of iterations to stop at
 * @param args::out_last_error receives the value of the error criterion at the final iteration
 * @param args::out_last_iteration receives the number of the final iteration
 * @param args::opt_error_evolution (optional) pointer to an array of length max_iterations, that receives the error
 * value of each iteration
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_NOT_CONVERGED if it hasn't reached given stopping criterion,
 * in case of failure it returns the associated error code
 */
jmtx_result JMTX_NAME_TYPED(solve_iterative_ilu_crs)(const JMTX_NAME_TYPED(matrix_crs) * mtx,
                                                     const JMTX_SCALAR_T *restrict y, JMTX_SCALAR_T *restrict x,
                                                     JMTX_SCALAR_T *restrict aux_vec,
                                                     JMTX_NAME_TYPED(solver_arguments) * args,
                                                     const jmtx_allocator_callbacks *allocator_callbacks)
{
    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }

    JMTX_NAME_TYPED(matrix_crs) * lower;
    JMTX_NAME_TYPED(matrix_ccs) * upper_ccs;
    jmtx_result res = JMTX_NAME_TYPED(decompose_ilu_crs)(mtx, &lower, &upper_ccs, allocator_callbacks);
    if (res != JMTX_RESULT_SUCCESS)
    {
        return res;
    }
    JMTX_NAME_TYPED(matrix_crs) * upper;
    res = JMTX_NAME_TYPED(convert_ccs_to_crs)(upper_ccs, &upper, allocator_callbacks);
    if (res != JMTX_RESULT_SUCCESS)
    {
        JMTX_NAME_TYPED(matrix_crs_destroy)(upper);
        JMTX_NAME_TYPED(matrix_crs_destroy)(lower);
        return res;
    }
    JMTX_NAME_TYPED(matrix_ccs_destroy)(upper_ccs);
    upper_ccs = NULL;

    res = JMTX_NAME_TYPED(solve_iterative_ilu_crs_precomputed)(mtx, lower, upper, y, x, aux_vec, args);

    JMTX_NAME_TYPED(matrix_crs_destroy)(upper);
    JMTX_NAME_TYPED(matrix_crs_destroy)(lower);

    return res;
}

/**
 * Solves the A x = L U x = y problem by computing the residual, then solving for L U e = r for the error e if residual
 * is too large to further refine the solution. This can help eliminate rounding errors, or alternatively can be used
 * with an incomplete LU decomposition (or ILU) to work as an iterative solver. Must be given the matrices LU.
 * @param a system matrix A
 * @param l lower triangular matrix with all entries on its main diagonal equal to 1
 * @param u upper triangular matrix
 * @param n dimensions of L, U, y, and x
 * @param y memory containing forcing vector
 * @param x memory which receives the solution
 * @param aux_vec auxiliary memory for a vector of the same size as x and y
 * @param args::in_convergence_criterion tolerance to determine if the solution is close enough
 * @param args::in_max_iterations number of iterations to stop at
 * @param args::out_last_error receives the value of the error criterion at the final iteration
 * @param args::out_last_iteration receives the number of the final iteration
 * @param args::opt_error_evolution (optional) pointer to an array of length max_iterations, that receives the error
 * value of each iteration
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_NOT_CONVERGED if it hasn't reached given stopping criterion,
 * in case of failure it returns the associated error code
 */
jmtx_result JMTX_NAME_TYPED(solve_iterative_ilu_crs_precomputed)(const JMTX_NAME_TYPED(matrix_crs) * mtx,
                                                                 const JMTX_NAME_TYPED(matrix_crs) * l,
                                                                 const JMTX_NAME_TYPED(matrix_crs) * u,
                                                                 const JMTX_SCALAR_T *restrict y,
                                                                 JMTX_SCALAR_T *restrict x, JMTX_SCALAR_T *aux_vec,
                                                                 JMTX_NAME_TYPED(solver_arguments) * args)
{
    const JMTX_INDEX_T n = mtx->base.rows;
    JMTX_SCALAR_T *const r = aux_vec;
    compute_residual(n, mtx, x, y, r);
    JMTX_REAL_T residual2 = 0;
    JMTX_REAL_T mag_y2 = 0;
    for (JMTX_INDEX_T i = 0; i < n; ++i)
    {
        residual2 += JMTX_DOT(r[i], r[i]);
        mag_y2 += JMTX_DOT(y[i], y[i]);
    }

    JMTX_REAL_T err = JMTX_REAL_ROOT(residual2 / mag_y2);
    const JMTX_REAL_T convergence_dif = args->in_convergence_criterion;
    const JMTX_INDEX_T n_max_iter = args->in_max_iterations;
    if (err < convergence_dif)
    {
        //  Converged prior to any iteration
        args->out_last_error = err;
        if (args->opt_error_evolution)
        {
            *args->opt_error_evolution = err;
        }
        args->out_last_iteration = 0;
        return JMTX_RESULT_SUCCESS;
    }

    //  Now begin iterations with x^{k + 1} = x^{k} + U'^{-1} L'^{-1} A x^{k}
    JMTX_SCALAR_T *x1 = aux_vec;
    JMTX_SCALAR_T *x2 = x;
    const JMTX_REAL_T stop_mag2 = mag_y2 * (convergence_dif * convergence_dif);
    JMTX_INDEX_T n_iter;
    for (n_iter = 0; n_iter < n_max_iter && stop_mag2 < residual2; ++n_iter)
    {
        {
            JMTX_SCALAR_T *tmp = x2;
            x2 = x1;
            x1 = tmp;
        }

        JMTX_NAME_TYPED(solve_direct_lu_crs_inplace)(l, u, r);

        for (JMTX_INDEX_T i = 0; i < n; ++i)
        {
            x[i] += r[i];
        }

        compute_residual(n, mtx, x, y, r);

        residual2 = 0;
        for (JMTX_INDEX_T i = 0; i < n; ++i)
        {
            residual2 += JMTX_DOT(r[i], r[i]);
        }

        if (args->opt_error_evolution)
        {
            args->opt_error_evolution[n_iter] = JMTX_REAL_ROOT(residual2 / mag_y2);
        }
    }
    err = JMTX_REAL_ROOT(residual2 / mag_y2);

    args->out_last_error = err;
    args->out_last_iteration = n_iter;

    return err < convergence_dif ? JMTX_RESULT_SUCCESS : JMTX_RESULT_NOT_CONVERGED;
}

/**
 * Solves the A x = L U x = y problem by computing the residual, then solving for L U e = r for the error e if residual
 * is too large to further refine the solution. This can help eliminate rounding errors, or alternatively can be used
 * with an incomplete LU decomposition (or ILU) to work as an iterative solver. Must be given the matrices LU.
 * This version offers a minor degree of parallelization, where calculation of the residual and applying error
 * correction are done in parallel, but the main bottleneck of dealing with inverting L U e = r is done in series.
 * @param a system matrix A
 * @param l lower triangular matrix with all entries on its main diagonal equal to 1
 * @param u upper triangular matrix
 * @param y memory containing forcing vector
 * @param x memory which receives the solution
 * @param aux_vec auxiliary memory for a vector of the same size as x and y
 * @param args::in_convergence_criterion tolerance to determine if the solution is close enough
 * @param args::in_max_iterations number of iterations to stop at
 * @param args::out_last_error receives the value of the error criterion at the final iteration
 * @param args::out_last_iteration receives the number of the final iteration
 * @param args::opt_error_evolution (optional) pointer to an array of length max_iterations, that receives the error
 * value of each iteration
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_NOT_CONVERGED if it hasn't reached given stopping criterion,
 * in case of failure it returns the associated error code
 */
jmtx_result JMTX_NAME_TYPED(solve_iterative_ilu_crs_precomputed_parallel)(
    const JMTX_NAME_TYPED(matrix_crs) * mtx, const JMTX_NAME_TYPED(matrix_crs) * l,
    const JMTX_NAME_TYPED(matrix_crs) * u, const JMTX_SCALAR_T *restrict y, JMTX_SCALAR_T *restrict x,
    JMTX_SCALAR_T *restrict aux_vec, JMTX_NAME_TYPED(solver_arguments) * args)
{
    const JMTX_INDEX_T n = mtx->base.rows;
    JMTX_SCALAR_T *const r = aux_vec;
    compute_residual(n, mtx, x, y, r);
    JMTX_REAL_T residual2 = 0;
    JMTX_REAL_T mag_y2 = 0;
    for (JMTX_INDEX_T i = 0; i < n; ++i)
    {
        residual2 += JMTX_DOT(r[i], r[i]);
        mag_y2 += JMTX_DOT(y[i], y[i]);
    }

    JMTX_REAL_T err = JMTX_REAL_ROOT(residual2 / mag_y2);
    const JMTX_REAL_T convergence_dif = args->in_convergence_criterion;
    const JMTX_INDEX_T n_max_iter = args->in_max_iterations;
    if (err < convergence_dif)
    {
        //  Converged prior to any iteration
        args->out_last_error = err;
        if (args->opt_error_evolution)
        {
            *args->opt_error_evolution = err;
        }

        args->out_last_iteration = 0;
        return JMTX_RESULT_SUCCESS;
    }

    //  Now begin iterations with x^{k + 1} = x^{k} + U'^{-1} L'^{-1} A x^{k}
    JMTX_SCALAR_T *x1 = aux_vec;
    JMTX_SCALAR_T *x2 = x;
    const JMTX_REAL_T stop_mag2 = mag_y2 * (convergence_dif * convergence_dif);
    JMTX_INDEX_T n_iter = 0;

    for (n_iter = 0; n_iter < n_max_iter && stop_mag2 < residual2; ++n_iter)
    {
        {
            JMTX_SCALAR_T *tmp = x2;
            x2 = x1;
            x1 = tmp;
        }

        JMTX_NAME_TYPED(solve_direct_lu_crs_inplace)(l, u, r);

        for (JMTX_INDEX_T i = 0; i < n; ++i)
        {
            x[i] += r[i];
        }

        compute_residual(n, mtx, x, y, r);

        residual2 = 0;
        for (JMTX_INDEX_T i = 0; i < n; ++i)
        {
            residual2 += JMTX_DOT(r[i], r[i]);
        }

        if (args->opt_error_evolution)
        {
            args->opt_error_evolution[n_iter] = JMTX_REAL_ROOT(residual2 / mag_y2);
        }
    }
    err = JMTX_REAL_ROOT(residual2 / mag_y2);

    args->out_last_error = err;
    args->out_last_iteration = n_iter;

    return err < convergence_dif ? JMTX_RESULT_SUCCESS : JMTX_RESULT_NOT_CONVERGED;
}

/**
 * Solves the A x = L U x = y problem by computing the residual, then solving for L U e = r for the error e if residual
 * is too large to further refine the solution. This can help eliminate rounding errors, or alternatively can be used
 * with an incomplete LU decomposition (or ILU) to work as an iterative solver. Must be given the matrices LU.
 * This version offers a minor degree of parallelization, where calculation of the residual and applying error
 * correction are done in parallel, but the main bottleneck of dealing with inverting L U e = r is done in series.
 * @param a system matrix A
 * @param y memory containing forcing vector
 * @param x memory which receives the solution
 * @param aux_vec auxiliary memory for a vector of the same size as x and y
 * @param args::in_convergence_criterion tolerance to determine if the solution is close enough
 * @param args::in_max_iterations number of iterations to stop at
 * @param args::out_last_error receives the value of the error criterion at the final iteration
 * @param args::out_last_iteration receives the number of the final iteration
 * @param args::opt_error_evolution (optional) pointer to an array of length max_iterations, that receives the error
 * value of each iteration
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_NOT_CONVERGED if it hasn't reached given stopping criterion,
 * in case of failure it returns the associated error code
 */
jmtx_result JMTX_NAME_TYPED(solve_iterative_ilu_crs_parallel)(const JMTX_NAME_TYPED(matrix_crs) * mtx,
                                                              const JMTX_SCALAR_T *restrict y,
                                                              JMTX_SCALAR_T *restrict x,
                                                              JMTX_SCALAR_T *restrict aux_vec,
                                                              JMTX_NAME_TYPED(solver_arguments) * args,
                                                              const jmtx_allocator_callbacks *allocator_callbacks)
{
    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }

    JMTX_NAME_TYPED(matrix_crs) * lower;
    JMTX_NAME_TYPED(matrix_ccs) * upper_ccs;
    jmtx_result res = JMTX_NAME_TYPED(decompose_ilu_crs)(mtx, &lower, &upper_ccs, allocator_callbacks);
    if (res != JMTX_RESULT_SUCCESS)
    {
        return res;
    }
    JMTX_NAME_TYPED(matrix_crs) * upper;
    res = JMTX_NAME_TYPED(convert_ccs_to_crs)(upper_ccs, &upper, allocator_callbacks);
    if (res != JMTX_RESULT_SUCCESS)
    {
        JMTX_NAME_TYPED(matrix_crs_destroy)(upper);
        JMTX_NAME_TYPED(matrix_crs_destroy)(lower);
        return res;
    }
    JMTX_NAME_TYPED(matrix_ccs_destroy)(upper_ccs);
    upper_ccs = NULL;

    res = JMTX_NAME_TYPED(solve_iterative_ilu_crs_precomputed_parallel)(mtx, lower, upper, y, x, aux_vec, args);

    JMTX_NAME_TYPED(matrix_crs_destroy)(upper);
    JMTX_NAME_TYPED(matrix_crs_destroy)(lower);

    return res;
}

/**
 * Solves a problem L U x = y, where L is a lower triangular matrix with the diagonal equal to 1 and U is an upper
 * triangular matrix.
 * @param decomposed LU decomposition of the system to solve
 * @param y memory containing forcing vector
 * @param x memory which receives the solution
 */
void JMTX_NAME_TYPED(solve_direct_lu_drm)(const JMTX_NAME_TYPED(matrix_drm) * decomposed,
                                          const JMTX_SCALAR_T *restrict y, JMTX_SCALAR_T *restrict x,
                                          JMTX_SCALAR_T *restrict aux_vec)
{
    const JMTX_FAST_INT_T n = decomposed->base.rows;
    if (!decomposed->permutations)
    {
        //  Decomposition has no pivoting, we can index directly
        //      Forward substitute for solving the L part
        for (JMTX_FAST_INT_T i = 0; i < n; ++i)
        {
            JMTX_SCALAR_T v = 0;
            for (JMTX_INDEX_T j = 0; j < i; ++j)
            {
                v += x[j] * decomposed->values[i * n + j];
            }
            x[i] = y[i] - v;
        }
        //      Backward substitute for solving the U part
        for (JMTX_FAST_INT_T i = 0; i < n; ++i)
        {
            JMTX_SCALAR_T v = 0;
            const JMTX_INDEX_T row_idx = n - 1 - i;
            for (JMTX_INDEX_T j = row_idx + 1; j < n; ++j)
            {
                v += x[j] * decomposed->values[row_idx * n + j];
            }
            x[row_idx] = (x[row_idx] - v) / decomposed->values[row_idx * n + row_idx];
        }
    }
    else
    {
        //  Decomposition has no pivoting, we can index directly
        //      Forward substitute for solving the L part
        for (JMTX_FAST_INT_T i = 0; i < n; ++i)
        {
            JMTX_SCALAR_T v = 0;
            for (JMTX_INDEX_T j = 0; j < i; ++j)
            {
                v += aux_vec[decomposed->permutations[j]] * decomposed->values[decomposed->permutations[i] * n + j];
            }
            aux_vec[decomposed->permutations[i]] = y[decomposed->permutations[i]] - v;
        }
        //      Backward substitute for solving the U part
        for (JMTX_FAST_INT_T i = 0; i < n; ++i)
        {
            JMTX_SCALAR_T v = 0;
            const JMTX_INDEX_T row_idx = n - 1 - i;
            for (JMTX_INDEX_T j = row_idx + 1; j < n; ++j)
            {
                v += x[j] * decomposed->values[decomposed->permutations[row_idx] * n + j];
            }
            x[row_idx] = (aux_vec[decomposed->permutations[row_idx]] - v) /
                         decomposed->values[decomposed->permutations[row_idx] * n + row_idx];
        }
    }
}
