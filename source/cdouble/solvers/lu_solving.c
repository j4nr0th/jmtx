// Automatically generated from source/float/solvers/lu_solving.c on Sun Dec 17 18:01:42 2023
//
// Created by jan on 6.11.2023.
//

#include "../../../include/jmtx/cdouble/solvers/lu_solving.h"
#include "../../../include/jmtx/cdouble/decompositions/incomplete_lu_decomposition.h"
#include "../../../include/jmtx/cdouble/matrices/sparse_conversion.h"
#include "../../../include/jmtx/cdouble/matrices/sparse_conversion_safe.h"
#include "../matrices/band_row_major_internal.h"
#include "../matrices/dense_row_major_internal.h"
#include "../matrices/sparse_row_compressed_internal.h"
#include <assert.h>
#include <complex.h>
#include <math.h>
/**
 * Solves a problem L U x = y, where L is a lower triangular matrix with the diagonal equal to 1 and U is an upper
 * triangular matrix.
 * @param l lower triangular matrix with all entries on its main diagonal equal to 1
 * @param u upper triangular matrix
 * @param y memory containing forcing vector
 * @param x memory which receives the solution
 */
void jmtxz_solve_direct_lu_brm(const jmtxz_matrix_brm *l, const jmtxz_matrix_brm *u, const _Complex double *restrict y,
                               _Complex double *restrict x)
{
    const uint_fast32_t n = l->base.cols;
    x[0] = y[0];
    //  First is the forward substitution for L v = y
    for (uint_fast32_t i = 1; i < n; ++i)
    {
        const uint_fast32_t off_j = jmtxz_matrix_brm_first_pos_in_row(l, i);
        _Complex double *values = NULL;
        const uint_fast32_t len = jmtxz_matrix_brm_get_row(l, i, &values);
        _Complex double v = 0;
        for (uint_fast32_t j = 0; j < len - 1; ++j)
        {
            v += values[j] * x[off_j + j];
        }
        x[i] = y[i] - v;
    }
    //  Then the backward substitution for U x = v
    for (int32_t i = (int32_t)n - 1; i >= 0; --i)
    {
        const uint_fast32_t off_j = jmtxz_matrix_brm_first_pos_in_row(u, i);
        _Complex double *values;
        const uint_fast32_t len = jmtxz_matrix_brm_get_row(u, i, &values);
        _Complex double v = 0;
        for (uint32_t j = 1; j < len; ++j)
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
void jmtxz_solve_direct_lu_brm_inplace(const jmtxz_matrix_brm *l, const jmtxz_matrix_brm *u,
                                       _Complex double *restrict x)
{
    const uint_fast32_t n = l->base.cols;
    //    x[0] = x[0];
    //  First is the forward substitution for L v = y
    for (uint_fast32_t i = 1; i < n; ++i)
    {
        const uint_fast32_t off_j = jmtxz_matrix_brm_first_pos_in_row(l, i);
        _Complex double *values = NULL;
        const uint_fast32_t len = jmtxz_matrix_brm_get_row(l, i, &values);
        _Complex double v = 0;
        for (uint_fast32_t j = 0; j < len - 1; ++j)
        {
            v += values[j] * x[off_j + j];
        }
        x[i] = x[i] - v;
    }
    //  Then the backward substitution for U x = v
    for (int32_t i = (int32_t)n - 1; i >= 0; --i)
    {
        const uint_fast32_t off_j = jmtxz_matrix_brm_first_pos_in_row(u, i);
        _Complex double *values;
        const uint_fast32_t len = jmtxz_matrix_brm_get_row(u, i, &values);
        _Complex double v = 0;
        for (uint32_t j = 1; j < len; ++j)
        {
            v += values[j] * x[off_j + j];
        }
        x[i] = (x[i] - v) / values[0];
    }
}

static inline int check_vector_overlaps(const unsigned n, const size_t size,
                                        const void *ptrs[JMTX_ARRAY_ATTRIB(static const n)])
{
    for (unsigned i = 0; i < n; ++i)
    {
        const uintptr_t p1 = (uintptr_t)ptrs[i];
        for (unsigned j = i + 1; j < n; ++j)
        {
            const uintptr_t p2 = (uintptr_t)ptrs[j];
            if (p1 > p2)
            {
                if (p2 + size > p1)
                {
                    return 1;
                }
            }
            else if (p1 + size > p2)
            {
                return 1;
            }
        }
    }

    return 0;
}

/**
 * Solves a problem L U x = y, where L is a lower triangular matrix with the diagonal equal to 1 and U is an upper
 * triangular matrix.
 * @param l lower triangular matrix with all entries on its main diagonal equal to 1
 * @param u upper triangular matrix
 * @param y memory containing forcing vector
 * @param x memory which receives the solution
 * @returns JMTX_RESULT_SUCCESS if successful, otherwise an error code indicating error in the input parameters
 */
jmtx_result jmtxzs_solve_direct_lu_brm(const jmtxz_matrix_brm *l, const jmtxz_matrix_brm *u, uint32_t n,
                                       const _Complex double y[JMTX_ARRAY_ATTRIB(static restrict n)],
                                       _Complex double x[JMTX_ARRAY_ATTRIB(static restrict n)])
{
    if (!l)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (l->base.type != JMTXZ_TYPE_BRM)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (l->base.rows != n || l->base.cols != n)
    {
        return JMTX_RESULT_BAD_MATRIX;
    }

    if (!u)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (u->base.type != JMTXZ_TYPE_BRM)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (u->base.rows != n || u->base.cols != n)
    {
        return JMTX_RESULT_BAD_MATRIX;
    }

    const void *ptrs[] = {x, y};
    if (check_vector_overlaps(sizeof(ptrs) / sizeof(*ptrs), sizeof(*x) * n, ptrs))
    {
        return JMTX_RESULT_BAD_PARAM;
    }
    jmtxz_solve_direct_lu_brm(l, u, y, x);
    return JMTX_RESULT_SUCCESS;
}

/**
 * Solves a problem L U x = y, where L is a lower triangular matrix with the diagonal equal to 1 and U is an upper
 * triangular matrix. This version of the function stores the solution vector x back into the same memory where the
 * forcing vector was.
 * @param l lower triangular matrix with all entries on its main diagonal equal to 1
 * @param u upper triangular matrix
 * @param y memory containing forcing vector
 * @param x memory which receives the solution
 * @returns JMTX_RESULT_SUCCESS if successful, otherwise an error code indicating error in the input parameters
 */
jmtx_result jmtxzs_solve_direct_lu_brm_inplace(const jmtxz_matrix_brm *l, const jmtxz_matrix_brm *u, uint32_t n,
                                               _Complex double x[JMTX_ARRAY_ATTRIB(static n)])
{
    if (!l)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (l->base.type != JMTXZ_TYPE_BRM)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (l->base.rows != n || l->base.cols != n)
    {
        return JMTX_RESULT_BAD_MATRIX;
    }

    if (!u)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (u->base.type != JMTXZ_TYPE_BRM)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (u->base.rows != n || u->base.cols != n)
    {
        return JMTX_RESULT_BAD_MATRIX;
    }

    jmtxz_solve_direct_lu_brm_inplace(l, u, x);
    return JMTX_RESULT_SUCCESS;
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
jmtx_result jmtxz_solve_iterative_lu_brm_refine(const jmtxz_matrix_brm *a, const jmtxz_matrix_brm *l,
                                                const jmtxz_matrix_brm *u,
                                                const _Complex double y[JMTX_ARRAY_ATTRIB(restrict)],
                                                _Complex double x[JMTX_ARRAY_ATTRIB(restrict)],
                                                _Complex double aux_vec[JMTX_ARRAY_ATTRIB(restrict)],
                                                jmtxd_solver_arguments *args)
{
    const uint_fast32_t n = l->base.cols;
    double y_magnitude2 = 0;
    uint_fast32_t iteration_count = 0;
    double error = 0;
    for (uint_fast32_t i = 0; i < n; ++i)
    {
        y_magnitude2 += conj(y[i]) * y[i];
    }
    jmtxz_solve_direct_lu_brm(l, u, y, x);

    for (;;)
    {
        for (uint_fast32_t i = 0; i < n; ++i)
        {
            aux_vec[i] = jmtxz_matrix_brm_vector_multiply_row(a, x, i);
        }

        for (uint_fast32_t i = 0; i < n; ++i)
        {
            aux_vec[i] = y[i] - aux_vec[i];
        }

        _Complex double residual_magnitude2 = 0;
        for (uint_fast32_t i = 0; i < n; ++i)
        {
            residual_magnitude2 += conj(aux_vec[i]) * aux_vec[i];
        }
        error = sqrt(residual_magnitude2 / y_magnitude2);
        if (args->opt_error_evolution)
        {
            args->opt_error_evolution[iteration_count] = error;
        }
        ++iteration_count;
        if (error < args->in_convergence_criterion || iteration_count > args->in_max_iterations)
        {
            break;
        }
        jmtxz_solve_direct_lu_brm_inplace(l, u, aux_vec);
        for (uint_fast32_t i = 0; i < n; ++i)
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
jmtx_result jmtxzs_solve_iterative_lu_brm_refine(const jmtxz_matrix_brm *a, const jmtxz_matrix_brm *l,
                                                 const jmtxz_matrix_brm *u, uint32_t n,
                                                 const _Complex double y[JMTX_ARRAY_ATTRIB(restrict static n)],
                                                 _Complex double x[JMTX_ARRAY_ATTRIB(restrict n)],
                                                 _Complex double aux_vec[JMTX_ARRAY_ATTRIB(restrict n)],
                                                 jmtxd_solver_arguments *args)
{
    double y_magnitude2 = 0;
    uint_fast32_t iteration_count = 0;
    double error = 0;
    for (uint_fast32_t i = 0; i < n; ++i)
    {
        y_magnitude2 += conj(y[i]) * y[i];
    }
    jmtx_result res = jmtxzs_solve_direct_lu_brm(l, u, n, y, x);
    if (res != JMTX_RESULT_SUCCESS)
    {
        return res;
    }

    for (;;)
    {
        for (uint_fast32_t i = 0; i < n; ++i)
        {
            aux_vec[i] = jmtxz_matrix_brm_vector_multiply_row(a, x, i);
        }

        for (uint_fast32_t i = 0; i < n; ++i)
        {
            aux_vec[i] = y[i] - aux_vec[i];
        }

        _Complex double residual_magnitude2 = 0;
        for (uint_fast32_t i = 0; i < n; ++i)
        {
            residual_magnitude2 += conj(aux_vec[i]) * aux_vec[i];
        }
        error = sqrt(residual_magnitude2 / y_magnitude2);
        if (args->opt_error_evolution)
        {
            args->opt_error_evolution[iteration_count] = error;
        }
        ++iteration_count;
        if (error < args->in_convergence_criterion || iteration_count > args->in_max_iterations)
        {
            break;
        }
        jmtxz_solve_direct_lu_brm_inplace(l, u, aux_vec);
        for (uint_fast32_t i = 0; i < n; ++i)
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
jmtx_result jmtxz_solve_iterative_lu_brm_refine_parallel(const jmtxz_matrix_brm *a, const jmtxz_matrix_brm *l,
                                                         const jmtxz_matrix_brm *u,
                                                         const _Complex double y[JMTX_ARRAY_ATTRIB(const restrict)],
                                                         _Complex double x[JMTX_ARRAY_ATTRIB(const restrict)],
                                                         _Complex double aux_vec[JMTX_ARRAY_ATTRIB(const restrict)],
                                                         jmtxd_solver_arguments *args)
{
    const uint_fast32_t n = l->base.cols;
    double y_magnitude2 = 0;
    uint_fast32_t iteration_count = 0;
    double error = 0;
    for (uint_fast32_t i = 0; i < n; ++i)
    {
        y_magnitude2 += conj(y[i]) * y[i];
    }
    double residual_magnitude2 = 0;
    jmtxz_solve_direct_lu_brm(l, u, y, x);
#pragma omp parallel default(none)                                                                                     \
    shared(y, x, l, u, aux_vec, args, n, error, y_magnitude2, iteration_count, residual_magnitude2, a)
    {
        for (;;)
        {
#pragma omp for
            for (uint_fast32_t i = 0; i < n; ++i)
            {
                aux_vec[i] = jmtxz_matrix_brm_vector_multiply_row(a, x, i);
            }

#pragma omp for
            for (uint_fast32_t i = 0; i < n; ++i)
            {
                aux_vec[i] = y[i] - aux_vec[i];
            }

#pragma omp for reduction(+ : residual_magnitude2)
            for (uint_fast32_t i = 0; i < n; ++i)
            {
                residual_magnitude2 += conj(aux_vec[i]) * aux_vec[i];
            }

#pragma omp single
            {
                error = sqrt(residual_magnitude2 / y_magnitude2);
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
                jmtxz_solve_direct_lu_brm_inplace(l, u, aux_vec);
            }
#pragma omp for
            for (uint_fast32_t i = 0; i < n; ++i)
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
void jmtxz_solve_direct_lu_crs(const jmtxz_matrix_crs *l, const jmtxz_matrix_crs *u, const _Complex double *restrict y,
                               _Complex double *restrict x)
{
    const uint32_t n = l->base.cols;
    x[0] = y[0];
    //  First is the forward substitution for L v = y
    for (uint32_t i = 1; i < n; ++i)
    {
        uint32_t *indices;
        _Complex double *values;
        uint32_t count = jmtxz_matrix_crs_get_row(l, i, &indices, &values);
        assert(indices[count - 1] == (uint32_t)i);
        assert(values[count - 1] == 1.0f);

        _Complex double v = 0;
        for (uint32_t j = 0; j < count - 1; ++j)
        {
            assert(indices[j] < i);
            v += values[j] * x[indices[j]];
        }
        x[i] = y[i] - v;
    }
    //  Then the backward substitution for U x = v
    for (int32_t i = (int32_t)n - 1; i >= 0; --i)
    {
        uint32_t *indices;
        _Complex double *values;
        uint32_t count = jmtxz_matrix_crs_get_row(u, i, &indices, &values);
        assert(indices[0] == (uint32_t)i);

        _Complex double v = 0;
        for (uint32_t j = 1; j < count; ++j)
        {
            assert(indices[j] > (uint32_t)i);
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
void jmtxz_solve_direct_lu_crs_inplace(const jmtxz_matrix_crs *l, const jmtxz_matrix_crs *u,
                                       _Complex double *restrict x)
{
    const uint32_t n = l->base.cols;
    //  First is the forward substitution for L v = y
    for (uint32_t i = 1; i < n; ++i)
    {
        uint32_t *indices;
        _Complex double *values;
        uint32_t count = jmtxz_matrix_crs_get_row(l, i, &indices, &values);
        assert(indices[count - 1] == (uint32_t)i);

        _Complex double v = 0;
        for (uint32_t j = 0; j < count - 1; ++j)
        {
            assert(indices[j] < i);
            v += values[j] * x[indices[j]];
        }
        x[i] = x[i] - v;
    }
    //  Then the backward substitution for U x = v
    for (int32_t i = (int32_t)n - 1; i >= 0; --i)
    {
        uint32_t *indices;
        _Complex double *values;
        uint32_t count = jmtxz_matrix_crs_get_row(u, i, &indices, &values);
        assert(indices[0] == (uint32_t)i);

        _Complex double v = 0;
        for (uint32_t j = 1; j < count; ++j)
        {
            assert(indices[j] > (uint32_t)i);
            v += values[j] * x[indices[j]];
        }
        x[i] = (x[i] - v) / values[0];
    }
}

/**
 * Solves a problem L U x = y, where L is a lower triangular matrix with the diagonal equal to 1 and U is an upper
 * triangular matrix.
 * @param l lower triangular matrix with all entries on its main diagonal equal to 1
 * @param u upper triangular matrix
 * @param n dimensions of L, U, y, and x
 * @param y memory containing forcing vector
 * @param x memory which receives the solution
 * @returns JMTX_RESULT_SUCCESS if successful, otherwise an error code indicating error in the input parameters
 */
jmtx_result jmtxzs_solve_direct_lu_crs(const jmtxz_matrix_crs *l, const jmtxz_matrix_crs *u, uint32_t n,
                                       const _Complex double y[JMTX_ARRAY_ATTRIB(static restrict n)],
                                       _Complex double x[JMTX_ARRAY_ATTRIB(restrict n)])
{
    if (!l)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (l->base.type != JMTXZ_TYPE_CRS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (l->base.rows != n || l->base.cols != n)
    {
        return JMTX_RESULT_BAD_MATRIX;
    }

    if (!u)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (u->base.type != JMTXZ_TYPE_CRS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (u->base.rows != n || u->base.cols != n)
    {
        return JMTX_RESULT_BAD_MATRIX;
    }

    if (!x)
    {
        return JMTX_RESULT_NULL_PARAM;
    }

    const void *ptrs[] = {x, y};
    if (check_vector_overlaps(sizeof(ptrs) / sizeof(*ptrs), sizeof(*x) * n, ptrs))
    {
        return JMTX_RESULT_BAD_PARAM;
    }
    jmtxz_solve_direct_lu_crs(l, u, y, x);
    return JMTX_RESULT_SUCCESS;
}

/**
 * Solves a problem L U x = y, where L is a lower triangular matrix with the diagonal equal to 1 and U is an upper
 * triangular matrix. This version of the function stores the solution vector x back into the same memory where the
 * forcing vector was.
 * @param l lower triangular matrix with all entries on its main diagonal equal to 1
 * @param u upper triangular matrix
 * @param n dimensions of L, U, y, and x
 * @param x memory which receives the solution
 * @returns JMTX_RESULT_SUCCESS if successful, otherwise an error code indicating error in the input parameters
 */
jmtx_result jmtxzs_solve_direct_lu_crs_inplace(const jmtxz_matrix_crs *l, const jmtxz_matrix_crs *u, uint32_t n,
                                               _Complex double x[JMTX_ARRAY_ATTRIB(static n)])
{
    if (!l)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (l->base.type != JMTXZ_TYPE_CRS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (l->base.rows != n || l->base.cols != n)
    {
        return JMTX_RESULT_BAD_MATRIX;
    }

    if (!u)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (u->base.type != JMTXZ_TYPE_CRS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (u->base.rows != n || u->base.cols != n)
    {
        return JMTX_RESULT_BAD_MATRIX;
    }

    jmtxz_solve_direct_lu_crs_inplace(l, u, x);
    return JMTX_RESULT_SUCCESS;
}

static inline void compute_residual(const uint32_t n, const jmtxz_matrix_crs *mtx, const _Complex double *restrict x,
                                    const _Complex double *restrict y, _Complex double *restrict r)
{
    for (uint32_t i = 0; i < n; ++i)
    {
        r[i] = y[i] - jmtxz_matrix_crs_vector_multiply_row(mtx, x, i);
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
jmtx_result jmtxz_solve_iterative_ilu_crs(const jmtxz_matrix_crs *mtx, const _Complex double *restrict y,
                                          _Complex double *restrict x, _Complex double *restrict aux_vec,
                                          jmtxd_solver_arguments *args,
                                          const jmtx_allocator_callbacks *allocator_callbacks)
{
    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }

    jmtxz_matrix_crs *lower;
    jmtxz_matrix_ccs *upper_ccs;
    jmtx_result res = jmtxz_decompose_ilu_crs(mtx, &lower, &upper_ccs, allocator_callbacks);
    if (res != JMTX_RESULT_SUCCESS)
    {
        return res;
    }
    jmtxz_matrix_crs *upper;
    res = jmtxz_convert_ccs_to_crs(upper_ccs, &upper, allocator_callbacks);
    if (res != JMTX_RESULT_SUCCESS)
    {
        jmtxz_matrix_crs_destroy(upper);
        jmtxz_matrix_crs_destroy(lower);
        return res;
    }
    jmtxz_matrix_ccs_destroy(upper_ccs);
    upper_ccs = NULL;

    res = jmtxz_solve_iterative_ilu_crs_precomputed(mtx, lower, upper, y, x, aux_vec, args);

    jmtxz_matrix_crs_destroy(upper);
    jmtxz_matrix_crs_destroy(lower);

    return res;
}

/**
 * Solves the A x = L U x = y problem by computing the residual, then solving for L U e = r for the error e if residual
 * is too large to further refine the solution. This can help eliminate rounding errors, or alternatively can be used
 * with an incomplete LU decomposition (or ILU) to work as an iterative solver. Must be given the matrices LU.
 * @param a system matrix A
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
jmtx_result jmtxzs_solve_iterative_ilu_crs(const jmtxz_matrix_crs *mtx, uint32_t n,
                                           const _Complex double y[JMTX_ARRAY_ATTRIB(restrict static n)],
                                           _Complex double x[JMTX_ARRAY_ATTRIB(restrict n)],
                                           _Complex double aux_vec[JMTX_ARRAY_ATTRIB(restrict n)],
                                           jmtxd_solver_arguments *args,
                                           const jmtx_allocator_callbacks *allocator_callbacks)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.rows != n || mtx->base.cols != n)
    {
        //  I am only doing square matrices!!!
        return JMTX_RESULT_BAD_MATRIX;
    }
    if (mtx->base.type != JMTXZ_TYPE_CRS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (!x)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!aux_vec)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    {
        const void *ptrs[] = {x, y, aux_vec};
        if (check_vector_overlaps(sizeof(ptrs) / sizeof(*ptrs), n * sizeof(*x), ptrs))
        {
            return JMTX_RESULT_BAD_PARAM;
        }
    }
    if (!args)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }
    else if (allocator_callbacks->alloc == NULL || allocator_callbacks->free == NULL)
    {
        return JMTX_RESULT_BAD_PARAM;
    }

    jmtxz_matrix_crs *lower;
    jmtxz_matrix_ccs *upper_ccs;
    jmtx_result res = jmtxzs_decompose_ilu_crs(mtx, &lower, &upper_ccs, allocator_callbacks);
    if (res != JMTX_RESULT_SUCCESS)
    {
        return res;
    }
    jmtxz_matrix_crs *upper;
    res = jmtxzs_convert_ccs_to_crs(upper_ccs, &upper, allocator_callbacks);
    if (res != JMTX_RESULT_SUCCESS)
    {
        jmtxz_matrix_crs_destroy(upper);
        jmtxz_matrix_crs_destroy(lower);
        return res;
    }
    jmtxz_matrix_ccs_destroy(upper_ccs);
    upper_ccs = NULL;

    res = jmtxzs_solve_iterative_ilu_crs_precomputed(mtx, lower, upper, n, y, x, aux_vec, args);

    jmtxz_matrix_crs_destroy(upper);
    jmtxz_matrix_crs_destroy(lower);

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
jmtx_result jmtxz_solve_iterative_ilu_crs_precomputed(const jmtxz_matrix_crs *mtx, const jmtxz_matrix_crs *l,
                                                      const jmtxz_matrix_crs *u, const _Complex double *restrict y,
                                                      _Complex double *restrict x, _Complex double *aux_vec,
                                                      jmtxd_solver_arguments *args)
{
    const uint32_t n = mtx->base.rows;
    _Complex double *const r = aux_vec;
    compute_residual(n, mtx, x, y, r);
    double residual2 = 0;
    double mag_y2 = 0;
    for (uint32_t i = 0; i < n; ++i)
    {
        residual2 += conj(r[i]) * r[i];
        mag_y2 += conj(y[i]) * y[i];
    }

    double err = sqrt(residual2 / mag_y2);
    const double convergence_dif = args->in_convergence_criterion;
    const uint32_t n_max_iter = args->in_max_iterations;
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
    _Complex double *x1 = aux_vec;
    _Complex double *x2 = x;
    const double stop_mag2 = mag_y2 * (convergence_dif * convergence_dif);
    uint32_t n_iter;
    for (n_iter = 0; n_iter < n_max_iter && stop_mag2 < residual2; ++n_iter)
    {
        {
            _Complex double *tmp = x2;
            x2 = x1;
            x1 = tmp;
        }

        jmtxz_solve_direct_lu_crs_inplace(l, u, r);

        for (uint32_t i = 0; i < n; ++i)
        {
            x[i] += r[i];
        }

        compute_residual(n, mtx, x, y, r);

        residual2 = 0;
        for (uint32_t i = 0; i < n; ++i)
        {
            residual2 += conj(r[i]) * r[i];
        }

        if (args->opt_error_evolution)
        {
            args->opt_error_evolution[n_iter] = sqrt(residual2 / mag_y2);
        }
    }
    err = sqrt(residual2 / mag_y2);

    args->out_last_error = err;
    args->out_last_iteration = n_iter;

    return err < convergence_dif ? JMTX_RESULT_SUCCESS : JMTX_RESULT_NOT_CONVERGED;
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
jmtx_result jmtxzs_solve_iterative_ilu_crs_precomputed(const jmtxz_matrix_crs *mtx, const jmtxz_matrix_crs *l,
                                                       const jmtxz_matrix_crs *u, uint32_t n,
                                                       const _Complex double y[JMTX_ARRAY_ATTRIB(restrict static n)],
                                                       _Complex double x[JMTX_ARRAY_ATTRIB(restrict n)],
                                                       _Complex double aux_vec[JMTX_ARRAY_ATTRIB(restrict n)],
                                                       jmtxd_solver_arguments *args)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!l)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!u)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.rows != n || mtx->base.cols != n)
    {
        //  I am only doing square matrices!!!
        return JMTX_RESULT_BAD_MATRIX;
    }
    if (l->base.rows != n || l->base.cols != n)
    {
        //  I am only doing square matrices!!!
        return JMTX_RESULT_BAD_MATRIX;
    }
    if (u->base.rows != n || u->base.cols != n)
    {
        //  I am only doing square matrices!!!
        return JMTX_RESULT_BAD_MATRIX;
    }
    if (mtx->base.type != JMTXZ_TYPE_CRS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (l->base.type != JMTXZ_TYPE_CRS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (u->base.type != JMTXZ_TYPE_CRS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (!x)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!aux_vec)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    {
        const void *ptrs[] = {x, y, aux_vec};
        if (check_vector_overlaps(sizeof(ptrs) / sizeof(*ptrs), n * sizeof(*x), ptrs))
        {
            return JMTX_RESULT_BAD_PARAM;
        }
    }
    if (!args)
    {
        return JMTX_RESULT_NULL_PARAM;
    }

    _Complex double *const r = aux_vec;
    compute_residual(n, mtx, x, y, r);
    double residual2 = 0;
    double mag_y2 = 0;
    for (uint32_t i = 0; i < n; ++i)
    {
        residual2 += conj(r[i]) * r[i];
        mag_y2 += conj(y[i]) * y[i];
    }

    double err = sqrt(residual2 / mag_y2);
    const double convergence_dif = args->in_convergence_criterion;
    const uint32_t n_max_iter = args->in_max_iterations;
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
    _Complex double *x1 = aux_vec;
    _Complex double *x2 = x;
    const double stop_mag2 = mag_y2 * (convergence_dif * convergence_dif);
    uint32_t n_iter;
    for (n_iter = 0; n_iter < n_max_iter && stop_mag2 < residual2; ++n_iter)
    {
        {
            _Complex double *tmp = x2;
            x2 = x1;
            x1 = tmp;
        }

        jmtxz_solve_direct_lu_crs_inplace(l, u, r);

        for (uint32_t i = 0; i < n; ++i)
        {
            x[i] += r[i];
        }

        compute_residual(n, mtx, x, y, r);

        residual2 = 0;
        for (uint32_t i = 0; i < n; ++i)
        {
            residual2 += conj(r[i]) * r[i];
        }

        if (args->opt_error_evolution)
        {
            args->opt_error_evolution[n_iter] = sqrt(residual2 / mag_y2);
        }
    }
    err = sqrt(residual2 / mag_y2);

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
jmtx_result jmtxz_solve_iterative_ilu_crs_precomputed_parallel(const jmtxz_matrix_crs *mtx, const jmtxz_matrix_crs *l,
                                                               const jmtxz_matrix_crs *u,
                                                               const _Complex double *restrict y,
                                                               _Complex double *restrict x,
                                                               _Complex double *restrict aux_vec,
                                                               jmtxd_solver_arguments *args)
{
    const uint32_t n = mtx->base.rows;
    _Complex double *const r = aux_vec;
    compute_residual(n, mtx, x, y, r);
    double residual2 = 0;
    double mag_y2 = 0;
    for (uint32_t i = 0; i < n; ++i)
    {
        residual2 += conj(r[i]) * r[i];
        mag_y2 += conj(y[i]) * y[i];
    }

    double err = sqrt(residual2 / mag_y2);
    const double convergence_dif = args->in_convergence_criterion;
    const uint32_t n_max_iter = args->in_max_iterations;
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
    _Complex double *x1 = aux_vec;
    _Complex double *x2 = x;
    const double stop_mag2 = mag_y2 * (convergence_dif * convergence_dif);
    uint32_t n_iter = 0;

    for (n_iter = 0; n_iter < n_max_iter && stop_mag2 < residual2; ++n_iter)
    {
        {
            _Complex double *tmp = x2;
            x2 = x1;
            x1 = tmp;
        }

        jmtxz_solve_direct_lu_crs_inplace(l, u, r);

        for (uint32_t i = 0; i < n; ++i)
        {
            x[i] += r[i];
        }

        compute_residual(n, mtx, x, y, r);

        residual2 = 0;
        for (uint32_t i = 0; i < n; ++i)
        {
            residual2 += conj(r[i]) * r[i];
        }

        if (args->opt_error_evolution)
        {
            args->opt_error_evolution[n_iter] = sqrt(residual2 / mag_y2);
        }
    }
    err = sqrt(residual2 / mag_y2);

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
jmtx_result jmtxzs_solve_iterative_ilu_crs_precomputed_parallel(
    const jmtxz_matrix_crs *mtx, const jmtxz_matrix_crs *l, const jmtxz_matrix_crs *u, uint32_t n,
    const _Complex double y[JMTX_ARRAY_ATTRIB(restrict static n)], _Complex double x[JMTX_ARRAY_ATTRIB(restrict n)],
    _Complex double aux_vec[JMTX_ARRAY_ATTRIB(restrict n)], jmtxd_solver_arguments *args)
{

    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!l)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!u)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.rows != n || mtx->base.cols != n)
    {
        //  I am only doing square matrices!!!
        return JMTX_RESULT_BAD_MATRIX;
    }
    if (l->base.rows != n || l->base.cols != n)
    {
        //  I am only doing square matrices!!!
        return JMTX_RESULT_BAD_MATRIX;
    }
    if (u->base.rows != n || u->base.cols != n)
    {
        //  I am only doing square matrices!!!
        return JMTX_RESULT_BAD_MATRIX;
    }
    if (mtx->base.type != JMTXZ_TYPE_CRS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (l->base.type != JMTXZ_TYPE_CRS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (u->base.type != JMTXZ_TYPE_CRS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (!x)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!aux_vec)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    {
        const void *ptrs[] = {x, y, aux_vec};
        if (check_vector_overlaps(sizeof(ptrs) / sizeof(*ptrs), n * sizeof(*x), ptrs))
        {
            return JMTX_RESULT_BAD_PARAM;
        }
    }
    if (!args)
    {
        return JMTX_RESULT_NULL_PARAM;
    }

    _Complex double *const r = aux_vec;
    compute_residual(n, mtx, x, y, r);
    double residual2 = 0;
    double mag_y2 = 0;
    for (uint32_t i = 0; i < n; ++i)
    {
        residual2 += conj(r[i]) * r[i];
        mag_y2 += conj(y[i]) * y[i];
    }

    double err = sqrt(residual2 / mag_y2);
    const double convergence_dif = args->in_convergence_criterion;
    const uint32_t n_max_iter = args->in_max_iterations;
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
    _Complex double *x1 = aux_vec;
    _Complex double *x2 = x;
    const double stop_mag2 = mag_y2 * (convergence_dif * convergence_dif);
    uint32_t n_iter = 0;

    for (n_iter = 0; n_iter < n_max_iter && stop_mag2 < residual2; ++n_iter)
    {
        {
            _Complex double *tmp = x2;
            x2 = x1;
            x1 = tmp;
        }

        jmtxz_solve_direct_lu_crs_inplace(l, u, r);

        for (uint32_t i = 0; i < n; ++i)
        {
            x[i] += r[i];
        }

        compute_residual(n, mtx, x, y, r);

        residual2 = 0;
        for (uint32_t i = 0; i < n; ++i)
        {
            residual2 += conj(r[i]) * r[i];
        }

        if (args->opt_error_evolution)
        {
            args->opt_error_evolution[n_iter] = sqrt(residual2 / mag_y2);
        }
    }
    err = sqrt(residual2 / mag_y2);

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
jmtx_result jmtxz_solve_iterative_ilu_crs_parallel(const jmtxz_matrix_crs *mtx, const _Complex double *restrict y,
                                                   _Complex double *restrict x, _Complex double *restrict aux_vec,
                                                   jmtxd_solver_arguments *args,
                                                   const jmtx_allocator_callbacks *allocator_callbacks)
{
    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }

    jmtxz_matrix_crs *lower;
    jmtxz_matrix_ccs *upper_ccs;
    jmtx_result res = jmtxz_decompose_ilu_crs(mtx, &lower, &upper_ccs, allocator_callbacks);
    if (res != JMTX_RESULT_SUCCESS)
    {
        return res;
    }
    jmtxz_matrix_crs *upper;
    res = jmtxz_convert_ccs_to_crs(upper_ccs, &upper, allocator_callbacks);
    if (res != JMTX_RESULT_SUCCESS)
    {
        jmtxz_matrix_crs_destroy(upper);
        jmtxz_matrix_crs_destroy(lower);
        return res;
    }
    jmtxz_matrix_ccs_destroy(upper_ccs);
    upper_ccs = NULL;

    res = jmtxz_solve_iterative_ilu_crs_precomputed_parallel(mtx, lower, upper, y, x, aux_vec, args);

    jmtxz_matrix_crs_destroy(upper);
    jmtxz_matrix_crs_destroy(lower);

    return res;
}

/**
 * Solves the A x = L U x = y problem by computing the residual, then solving for L U e = r for the error e if residual
 * is too large to further refine the solution. This can help eliminate rounding errors, or alternatively can be used
 * with an incomplete LU decomposition (or ILU) to work as an iterative solver. Must be given the matrices LU.
 * This version offers a minor degree of parallelization, where calculation of the residual and applying error
 * correction are done in parallel, but the main bottleneck of dealing with inverting L U e = r is done in series.
 * @param a system matrix A
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
jmtx_result jmtxzs_solve_iterative_ilu_crs_parallel(const jmtxz_matrix_crs *mtx, uint32_t n,
                                                    const _Complex double y[JMTX_ARRAY_ATTRIB(restrict static n)],
                                                    _Complex double x[JMTX_ARRAY_ATTRIB(restrict n)],
                                                    _Complex double aux_vec[JMTX_ARRAY_ATTRIB(restrict n)],
                                                    jmtxd_solver_arguments *args,
                                                    const jmtx_allocator_callbacks *allocator_callbacks)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.rows != n || mtx->base.cols != n)
    {
        //  I am only doing square matrices!!!
        return JMTX_RESULT_BAD_MATRIX;
    }
    if (mtx->base.type != JMTXZ_TYPE_CRS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (!x)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (!aux_vec)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    {
        const void *ptrs[] = {x, y, aux_vec};
        if (check_vector_overlaps(sizeof(ptrs) / sizeof(*ptrs), n * sizeof(*x), ptrs))
        {
            return JMTX_RESULT_BAD_PARAM;
        }
    }
    if (!args)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (allocator_callbacks == NULL)
    {
        allocator_callbacks = &JMTX_DEFAULT_ALLOCATOR_CALLBACKS;
    }
    else if (allocator_callbacks->alloc == NULL || allocator_callbacks->free == NULL)
    {
        return JMTX_RESULT_BAD_PARAM;
    }

    jmtxz_matrix_crs *lower;
    jmtxz_matrix_ccs *upper_ccs;
    jmtx_result res = jmtxzs_decompose_ilu_crs(mtx, &lower, &upper_ccs, allocator_callbacks);
    if (res != JMTX_RESULT_SUCCESS)
    {
        return res;
    }
    jmtxz_matrix_crs *upper;
    res = jmtxz_convert_ccs_to_crs(upper_ccs, &upper, allocator_callbacks);
    if (res != JMTX_RESULT_SUCCESS)
    {
        jmtxz_matrix_crs_destroy(upper);
        jmtxz_matrix_crs_destroy(lower);
        return res;
    }
    jmtxz_matrix_ccs_destroy(upper_ccs);
    upper_ccs = NULL;

    res = jmtxzs_solve_iterative_ilu_crs_precomputed_parallel(mtx, lower, upper, n, y, x, aux_vec, args);

    jmtxz_matrix_crs_destroy(upper);
    jmtxz_matrix_crs_destroy(lower);

    return res;
}

/**
 * Solves a problem L U x = y, where L is a lower triangular matrix with the diagonal equal to 1 and U is an upper
 * triangular matrix.
 * @param decomposed LU decomposition of the system to solve
 * @param y memory containing forcing vector
 * @param x memory which receives the solution
 */
void jmtxz_solve_direct_lu_drm(const jmtxz_matrix_drm *decomposed, const _Complex double *restrict y,
                               _Complex double *restrict x, _Complex double *restrict aux_vec)
{
    const uint_fast32_t n = decomposed->base.rows;
    if (!decomposed->permutations)
    {
        //  Decomposition has no pivoting, we can index directly
        //      Forward substitute for solving the L part
        for (uint_fast32_t i = 0; i < n; ++i)
        {
            _Complex double v = 0;
            for (uint32_t j = 0; j < i; ++j)
            {
                v += x[j] * decomposed->values[i * n + j];
            }
            x[i] = y[i] - v;
        }
        //      Backward substitute for solving the U part
        for (uint_fast32_t i = 0; i < n; ++i)
        {
            _Complex double v = 0;
            const uint32_t row_idx = n - 1 - i;
            for (uint32_t j = row_idx + 1; j < n; ++j)
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
        for (uint_fast32_t i = 0; i < n; ++i)
        {
            _Complex double v = 0;
            for (uint32_t j = 0; j < i; ++j)
            {
                v += aux_vec[decomposed->permutations[j]] * decomposed->values[decomposed->permutations[i] * n + j];
            }
            aux_vec[decomposed->permutations[i]] = y[decomposed->permutations[i]] - v;
        }
        //      Backward substitute for solving the U part
        for (uint_fast32_t i = 0; i < n; ++i)
        {
            _Complex double v = 0;
            const uint32_t row_idx = n - 1 - i;
            for (uint32_t j = row_idx + 1; j < n; ++j)
            {
                v += x[j] * decomposed->values[decomposed->permutations[row_idx] * n + j];
            }
            x[row_idx] = (aux_vec[decomposed->permutations[row_idx]] - v) /
                         decomposed->values[decomposed->permutations[row_idx] * n + row_idx];
        }
    }
}
