//
// Created by jan on 1.1.2024.
//
#include <math.h>

#include "recursive_generalized_minimum_residual_iteration.h"
#include "../matrices/band_row_major.h"
#include "../matrices/sparse_diagonal_compressed.h"
#include "../matrices/sparse_row_compressed.h"

#include "gmres_internal.h"

/**
 * Applies Generalized Minimum Residual Recursive method (known as GMRESR) to solve a linear system A x = y.
 * Builds up a set of m orthonormal basis for the Krylov subspace in order to find an optimal search direction, then
 * makes its contribution to residual orthogonal to previous search directions. This makes it in general much better
 * behaved than GMRES with a short restart interval. Previously used search directions are saved up to the specified
 * number and removed from the previously applied directions.
 *
 * @param mtx system matrix A
 * @param y the solution of the system A x = y
 * @param x the solution vector which contains the initial guess of the solution
 * @param m the GMRES restart interval
 * @param l the CGR truncation interval
 * @param r an m by m upper triangular matrix (lbw = 0, ubw = m - 1) that is to be used in solving the least squares
 * problem
 * @param aux_vec1 auxiliary memory for a vector of m elements
 * @param aux_vec2 auxiliary memory for a vector of m elements
 * @param aux_vec3 auxiliary memory for a vector of m elements
 * @param aux_vec4 auxiliary memory for a vector of m elements
 * @param aux_vec5 auxiliary memory for a vector of m elements
 * @param aux_vec6 auxiliary memory for a vector of the same size as x and y
 * @param aux_vecs1 auxiliary memory for m vectors of the same size as x and y (n by m)
 * @param aux_vecs2 auxiliary memory for l vectors of the same size as x and y (n by m)
 * @param aux_vecs3 auxiliary memory for l vectors of the same size as x and y (n by m)
 * @param args::in_convergence_criterion tolerance to determine if the solution is close enough
 * @param args::in_max_iterations number of iterations to stop at
 * @param args::out_last_error receives the value of the error criterion at the final iteration
 * @param args::out_last_iteration receives the number of the final iteration
 * @param args::opt_error_evolution (optional) pointer to an array of length max_iterations, that receives the error
 * value of each iteration
 * @return JMTX_RESULT_SUCCESS if solution converged, JMTX_RESULT_NOT_CONVERGED if solution did not converge in the
 * given number of iterations, other error codes for other errors
 */
jmtx_result JMTX_NAME_TYPED(solve_iterative_gmresr_cds)(
    const JMTX_NAME_TYPED(matrix_cds) * mtx, const JMTX_SCALAR_T *restrict y, JMTX_SCALAR_T *restrict x, JMTX_INDEX_T m,
    JMTX_INDEX_T l, JMTX_NAME_TYPED(matrix_brm) * r_mtx, JMTX_SCALAR_T aux_vec1[JMTX_ARRAY_ATTRIB(restrict m)],
    JMTX_SCALAR_T aux_vec2[JMTX_ARRAY_ATTRIB(restrict m)], JMTX_SCALAR_T aux_vec3[JMTX_ARRAY_ATTRIB(restrict m)],
    JMTX_SCALAR_T aux_vec4[JMTX_ARRAY_ATTRIB(restrict m)], JMTX_SCALAR_T aux_vec5[JMTX_ARRAY_ATTRIB(restrict m)],
    JMTX_SCALAR_T *restrict aux_vec6, JMTX_SCALAR_T *restrict aux_vecs1, JMTX_SCALAR_T *restrict aux_vecs2,
    JMTX_SCALAR_T *restrict aux_vecs3, JMTX_NAME_TYPED(solver_arguments) * args)
{
    JMTX_REAL_T err, r_mag, y_mag;
    const JMTX_INDEX_T n = mtx->base.rows;
    JMTX_SCALAR_T *const r = aux_vec6;
    JMTX_NAME_TYPED(matrix_cds_vector_multiply)(mtx, x, r);
    y_mag = 0;
    r_mag = 0;
    for (JMTX_INDEX_T i = 0; i < n; ++i)
    {
        const JMTX_SCALAR_T res = y[i] - r[i];
        y_mag += JMTX_DOT(y[i], y[i]);
        r_mag += JMTX_DOT(res, res);
        r[i] = res;
    }
    y_mag = JMTX_REAL_ROOT(y_mag);
    r_mag = JMTX_REAL_ROOT(r_mag);
    err = r_mag / y_mag;
    if (err < args->in_convergence_criterion || args->in_max_iterations == 0)
    {
        args->out_last_iteration = 0;
        args->out_last_error = err;
        return JMTX_RESULT_SUCCESS;
    }

    JMTX_INDEX_T n_iter = 0, vec_idx = 0, vec_cnt = 0;
    JMTX_SCALAR_T *const s_vectors = aux_vecs2;
    JMTX_SCALAR_T *const p_vectors = aux_vecs3;
    //  Begin iterations
    for (;;)
    {
        JMTX_SCALAR_T *const s = s_vectors + n * vec_idx;
        JMTX_SCALAR_T *const p = p_vectors + n * vec_idx;
        for (JMTX_INDEX_T i = 0; i < n; ++i)
        {
            p[i] = 0;
        }
        //  Use GMRES(m) to generate a correction to the current error
        (void)JMTX_NAME_TYPED(gmresm_round_cds)(mtx, n, m, y_mag, args->in_convergence_criterion, r, p, r_mtx, aux_vec1,
                                                aux_vec2, aux_vec3, aux_vec4, aux_vec5, aux_vecs1);
        //  Generate conjugate vector
        JMTX_NAME_TYPED(matrix_cds_vector_multiply)(mtx, p, s);

        //  Make conjugate vector orthonormal to previous l search vectors
        for (JMTX_INDEX_T i = 0; i < vec_cnt; ++i)
        {
            if (i == vec_idx)
            {
                continue;
            }
            JMTX_SCALAR_T dp = 0;
            const JMTX_SCALAR_T *other_s = s_vectors + i * n;
            const JMTX_SCALAR_T *other_p = p_vectors + i * n;
            for (JMTX_INDEX_T j = 0; j < n; ++j)
            {
                dp += JMTX_DOT(other_s[j], s[j]);
            }
            for (JMTX_INDEX_T j = 0; j < n; ++j)
            {
                s[j] -= dp * other_s[j];
                p[j] -= dp * other_p[j];
            }
        }
        JMTX_REAL_T s_mag = 0;
        for (JMTX_INDEX_T j = 0; j < n; ++j)
        {
            s_mag += JMTX_DOT(s[j], s[j]);
        }
        s_mag = JMTX_REAL_ROOT(s_mag);
        for (JMTX_INDEX_T j = 0; j < n; ++j)
        {
            s[j] /= s_mag;
            p[j] /= s_mag;
        }

        JMTX_SCALAR_T alpha = 0;
        for (JMTX_INDEX_T j = 0; j < n; ++j)
        {
            alpha += JMTX_DOT(r[j], s[j]);
        }

        r_mag = 0;
        for (JMTX_INDEX_T j = 0; j < n; ++j)
        {
            r[j] -= alpha * s[j];
            r_mag += JMTX_DOT(r[j], r[j]);
            x[j] += alpha * p[j];
        }
        r_mag = JMTX_REAL_ROOT(r_mag);
        err = r_mag / y_mag;

        if (vec_cnt != l)
        {
            vec_cnt += 1;
        }
        vec_idx += 1;
        if (vec_idx == vec_cnt)
        {
            vec_idx = 0;
        }

        if (args->opt_error_evolution)
        {
            args->opt_error_evolution[n_iter] = err;
        }
        n_iter += 1;
        if (n_iter == args->in_max_iterations || err < args->in_convergence_criterion)
        {
            break;
        }
    }

    args->out_last_error = err;
    args->out_last_iteration = n_iter;
    if (isfinite(err) && err < args->in_convergence_criterion)
    {
        return JMTX_RESULT_SUCCESS;
    }

    return JMTX_RESULT_NOT_CONVERGED;
}

/**
 * Applies Generalized Minimum Residual Recursive method (known as GMRESR) to solve a linear system A x = y.
 * Builds up a set of m orthonormal basis for the Krylov subspace in order to find an optimal search direction, then
 * makes its contribution to residual orthogonal to previous search directions. This makes it in general much better
 * behaved than GMRES with a short restart interval. Previously used search directions are saved up to the specified
 * number and removed from the previously applied directions.
 *
 * @param mtx system matrix A
 * @param y the solution of the system A x = y
 * @param x the solution vector which contains the initial guess of the solution
 * @param m the GMRES restart interval
 * @param l the CGR truncation interval
 * @param r an m by m upper triangular matrix (lbw = 0, ubw = m - 1) that is to be used in solving the least squares
 * problem
 * @param aux_vec1 auxiliary memory for a vector of m elements
 * @param aux_vec2 auxiliary memory for a vector of m elements
 * @param aux_vec3 auxiliary memory for a vector of m elements
 * @param aux_vec4 auxiliary memory for a vector of m elements
 * @param aux_vec5 auxiliary memory for a vector of m elements
 * @param aux_vec6 auxiliary memory for a vector of the same size as x and y
 * @param aux_vecs1 auxiliary memory for m vectors of the same size as x and y (n by m)
 * @param aux_vecs2 auxiliary memory for l vectors of the same size as x and y (n by m)
 * @param aux_vecs3 auxiliary memory for l vectors of the same size as x and y (n by m)
 * @param args::in_convergence_criterion tolerance to determine if the solution is close enough
 * @param args::in_max_iterations number of iterations to stop at
 * @param args::out_last_error receives the value of the error criterion at the final iteration
 * @param args::out_last_iteration receives the number of the final iteration
 * @param args::opt_error_evolution (optional) pointer to an array of length max_iterations, that receives the error
 * value of each iteration
 * @return JMTX_RESULT_SUCCESS if solution converged, JMTX_RESULT_NOT_CONVERGED if solution did not converge in the
 * given number of iterations, other error codes for other errors
 */
jmtx_result JMTX_NAME_TYPED(solve_iterative_gmresr_crs)(
    const JMTX_NAME_TYPED(matrix_crs) * mtx, const JMTX_SCALAR_T *restrict y, JMTX_SCALAR_T *restrict x, JMTX_INDEX_T m,
    JMTX_INDEX_T l, JMTX_NAME_TYPED(matrix_brm) * r_mtx, JMTX_SCALAR_T aux_vec1[JMTX_ARRAY_ATTRIB(restrict m)],
    JMTX_SCALAR_T aux_vec2[JMTX_ARRAY_ATTRIB(restrict m)], JMTX_SCALAR_T aux_vec3[JMTX_ARRAY_ATTRIB(restrict m)],
    JMTX_SCALAR_T aux_vec4[JMTX_ARRAY_ATTRIB(restrict m)], JMTX_SCALAR_T aux_vec5[JMTX_ARRAY_ATTRIB(restrict m)],
    JMTX_SCALAR_T *restrict aux_vec6, JMTX_SCALAR_T *restrict aux_vecs1, JMTX_SCALAR_T *restrict aux_vecs2,
    JMTX_SCALAR_T *restrict aux_vecs3, JMTX_NAME_TYPED(solver_arguments) * args)
{
    JMTX_REAL_T err, r_mag, y_mag;
    const JMTX_INDEX_T n = mtx->base.rows;
    JMTX_SCALAR_T *const r = aux_vec6;
    JMTX_NAME_TYPED(matrix_crs_vector_multiply)(mtx, x, r);
    y_mag = 0;
    r_mag = 0;
    for (JMTX_INDEX_T i = 0; i < n; ++i)
    {
        const JMTX_SCALAR_T res = y[i] - r[i];
        y_mag += JMTX_DOT(y[i], y[i]);
        r_mag += JMTX_DOT(res, res);
        r[i] = res;
    }
    y_mag = JMTX_REAL_ROOT(y_mag);
    r_mag = JMTX_REAL_ROOT(r_mag);
    err = r_mag / y_mag;
    if (err < args->in_convergence_criterion || args->in_max_iterations == 0)
    {
        args->out_last_iteration = 0;
        args->out_last_error = err;
        return JMTX_RESULT_SUCCESS;
    }

    JMTX_INDEX_T n_iter = 0, vec_idx = 0, vec_cnt = 0;
    JMTX_SCALAR_T *const s_vectors = aux_vecs2;
    JMTX_SCALAR_T *const p_vectors = aux_vecs3;
    //  Begin iterations
    for (;;)
    {
        JMTX_SCALAR_T *const s = s_vectors + n * vec_idx;
        JMTX_SCALAR_T *const p = p_vectors + n * vec_idx;
        for (JMTX_INDEX_T i = 0; i < n; ++i)
        {
            p[i] = 0;
        }
        //  Use GMRES(m) to generate a correction to the current error
        (void)JMTX_NAME_TYPED(gmresm_round_crs)(mtx, n, m, y_mag, args->in_convergence_criterion, r, p, r_mtx, aux_vec1,
                                                aux_vec2, aux_vec3, aux_vec4, aux_vec5, aux_vecs1);
        //  Generate conjugate vector
        JMTX_NAME_TYPED(matrix_crs_vector_multiply)(mtx, p, s);

        //  Make conjugate vector orthonormal to previous l search vectors
        for (JMTX_INDEX_T i = 0; i < vec_cnt; ++i)
        {
            if (i == vec_idx)
            {
                continue;
            }
            JMTX_SCALAR_T dp = 0;
            const JMTX_SCALAR_T *other_s = s_vectors + i * n;
            const JMTX_SCALAR_T *other_p = p_vectors + i * n;
            for (JMTX_INDEX_T j = 0; j < n; ++j)
            {
                dp += JMTX_DOT(other_s[j], s[j]);
            }
            for (JMTX_INDEX_T j = 0; j < n; ++j)
            {
                s[j] -= dp * other_s[j];
                p[j] -= dp * other_p[j];
            }
        }
        JMTX_REAL_T s_mag = 0;
        for (JMTX_INDEX_T j = 0; j < n; ++j)
        {
            s_mag += JMTX_DOT(s[j], s[j]);
        }
        s_mag = JMTX_REAL_ROOT(s_mag);
        for (JMTX_INDEX_T j = 0; j < n; ++j)
        {
            s[j] /= s_mag;
            p[j] /= s_mag;
        }

        JMTX_SCALAR_T alpha = 0;
        for (JMTX_INDEX_T j = 0; j < n; ++j)
        {
            alpha += JMTX_DOT(r[j], s[j]);
        }

        r_mag = 0;
        for (JMTX_INDEX_T j = 0; j < n; ++j)
        {
            r[j] -= alpha * s[j];
            r_mag += JMTX_DOT(r[j], r[j]);
            x[j] += alpha * p[j];
        }
        r_mag = JMTX_REAL_ROOT(r_mag);
        err = r_mag / y_mag;

        if (vec_cnt != l)
        {
            vec_cnt += 1;
        }
        vec_idx += 1;
        if (vec_idx == vec_cnt)
        {
            vec_idx = 0;
        }

        if (args->opt_error_evolution)
        {
            args->opt_error_evolution[n_iter] = err;
        }
        n_iter += 1;
        if (n_iter == args->in_max_iterations || err < args->in_convergence_criterion)
        {
            break;
        }
    }

    args->out_last_error = err;
    args->out_last_iteration = n_iter;
    if (isfinite(err) && err < args->in_convergence_criterion)
    {
        return JMTX_RESULT_SUCCESS;
    }

    return JMTX_RESULT_NOT_CONVERGED;
}
