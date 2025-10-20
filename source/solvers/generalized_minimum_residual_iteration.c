//
// Created by jan on 1.1.2024.
//

#include "generalized_minimum_residual_iteration.h"
#include "../matrices/band_row_major.h"
#include "../matrices/sparse_diagonal_compressed.h"
#include "../matrices/sparse_row_compressed.h"
#include <assert.h>
#include <math.h>

jmtx_result JMTX_NAME_TYPED(solve_iterative_gmresm_crs)(
    const JMTX_NAME_TYPED(matrix_crs) * mtx, const JMTX_SCALAR_T *restrict y, JMTX_SCALAR_T *restrict x, JMTX_INDEX_T m,
    JMTX_NAME_TYPED(matrix_brm) * r, JMTX_SCALAR_T aux_vec1[JMTX_ARRAY_ATTRIB(restrict m)],
    JMTX_SCALAR_T aux_vec2[JMTX_ARRAY_ATTRIB(restrict m)], JMTX_SCALAR_T aux_vec3[JMTX_ARRAY_ATTRIB(restrict m)],
    JMTX_SCALAR_T aux_vec4[JMTX_ARRAY_ATTRIB(restrict m)], JMTX_SCALAR_T aux_vec5[JMTX_ARRAY_ATTRIB(restrict m)],
    JMTX_SCALAR_T *restrict aux_vecs, JMTX_NAME_TYPED(solver_arguments) * args)
{
    JMTX_REAL_T err = 0, y_mag = 0, r_mag = 0;
    JMTX_INDEX_T n_iteration = 0;
    const JMTX_INDEX_T n = mtx->base.rows;
    for (JMTX_INDEX_T i = 0; i < n; ++i)
    {
        y_mag += JMTX_DOT(y[i], y[i]);
    }
    y_mag = JMTX_REAL_ROOT(y_mag);

    const JMTX_INDEX_T round_count = args->in_max_iterations / m;
    for (JMTX_INDEX_T round = 0; round < round_count; ++round)
    {
        JMTX_SCALAR_T *const q = aux_vecs;
        JMTX_SCALAR_T *const ck = aux_vec1;
        JMTX_SCALAR_T *const sk = aux_vec2;
        JMTX_SCALAR_T *const g = aux_vec3;
        JMTX_SCALAR_T *const alpha = aux_vec4;
        JMTX_SCALAR_T *const h = aux_vec5;
        JMTX_SCALAR_T *p = q;
        for (JMTX_INDEX_T i = 0; i < n; ++i)
        {
            const JMTX_SCALAR_T res = y[i] - JMTX_NAME_TYPED(matrix_crs_vector_multiply_row)(mtx, x, i);
            r_mag += JMTX_DOT(res, res);
            p[i] = res;
        }
        r_mag = JMTX_REAL_ROOT(r_mag);
        err = r_mag / y_mag;
        if (err < args->in_convergence_criterion)
        {
            args->out_last_iteration = n_iteration;
            args->out_last_error = err;
            return JMTX_RESULT_SUCCESS;
        }
        for (JMTX_INDEX_T i = 0; i < n; ++i)
        {
            p[i] /= r_mag;
        }
        g[0] = r_mag;
        JMTX_INDEX_T k;
        for (k = 1; k < m; ++k)
        {
            //  Generate new basis vector
            JMTX_NAME_TYPED(matrix_crs_vector_multiply)(mtx, p, p + n);
            p += n;
            //  Make the basis orthogonal to other basis
            for (JMTX_INDEX_T l = 0; l < k; ++l)
            {
                h[l] = 0;
                const JMTX_SCALAR_T *old_p = q + n * l;
                for (JMTX_INDEX_T i = 0; i < n; ++i)
                {
                    h[l] += JMTX_DOT(old_p[i], p[i]);
                }
                for (JMTX_INDEX_T i = 0; i < n; ++i)
                {
                    p[i] -= JMTX_DOT(h[l], old_p[i]);
                }
            }

            //  Find magnitude
            JMTX_REAL_T mag_p = 0;
            for (JMTX_INDEX_T i = 0; i < n; ++i)
            {
                mag_p += JMTX_DOT(p[i], p[i]);
            }
            const JMTX_SCALAR_T mag_p2 = mag_p;
            mag_p = JMTX_REAL_ROOT(mag_p);

            //  Normalize basis vector
            for (JMTX_INDEX_T i = 0; i < n; ++i)
            {
                p[i] /= mag_p;
            }

            //  Apply previous Givens rotations to the new column of R
            for (JMTX_INDEX_T l = 0; l < k - 1; ++l)
            {
                const JMTX_SCALAR_T tmp = ck[l] * h[l] + sk[l] * h[l + 1];
                h[l + 1] = -sk[l] * h[l] + ck[l] * h[l + 1];
                h[l] = tmp;
            }

            //  Compute the new givens rotation
            const JMTX_REAL_T rho = JMTX_REAL_ROOT(mag_p2 + JMTX_DOT(h[k - 1], h[k - 1]));
            const JMTX_SCALAR_T c_new = h[k - 1] / rho;
            const JMTX_REAL_T s_new = mag_p / rho;
            ck[k - 1] = c_new;
            sk[k - 1] = s_new;
            h[k - 1] = c_new * h[k - 1] + s_new * mag_p;

            JMTX_NAME_TYPED(matrix_brm_set_col)(r, k - 1, h);

            g[k] = -s_new * g[k - 1];
            g[k - 1] = c_new * g[k - 1];

            r_mag = JMTX_ABS(g[k]);
            err = r_mag / y_mag;
            if (args->opt_error_evolution)
            {
                args->opt_error_evolution[n_iteration] = err;
            }
            n_iteration += 1;
            if (err < args->in_convergence_criterion || k + 1 == m)
            {
                break;
            }
        }

        //  Solve the least squares problem using back substitution
        for (JMTX_FAST_INT_T row = 0; row < k; ++row)
        {
            const JMTX_FAST_INT_T i = k - 1 - row;
            JMTX_SCALAR_T *elements;
            JMTX_NAME_TYPED(matrix_brm_get_row)(r, i, &elements);
            JMTX_SCALAR_T sum = 0;
            for (JMTX_FAST_INT_T j = 1; j < row + 1; ++j)
            {
                sum += JMTX_DOT(alpha[i + j], elements[j]);
            }
            alpha[i] = (g[i] - sum) / elements[0];
        }

        //  Compute improvement to x
        for (JMTX_INDEX_T i = 0; i < n; ++i)
        {
            for (JMTX_INDEX_T j = 0; j < k; ++j)
            {
                x[i] += alpha[j] * q[j * n + i];
            }
        }
        if (err < args->in_convergence_criterion)
        {
            args->out_last_error = err;
            args->out_last_iteration = n_iteration;
            return JMTX_RESULT_SUCCESS;
        }
    }
    args->out_last_error = err;
    args->out_last_iteration = n_iteration;
    return JMTX_RESULT_NOT_CONVERGED;
}

JMTX_INDEX_T JMTX_NAME_TYPED(gmresm_round_cds)(
    const JMTX_NAME_TYPED(matrix_cds) * mtx, const JMTX_INDEX_T n, const JMTX_INDEX_T m, const JMTX_REAL_T y_mag,
    const JMTX_REAL_T tol, const JMTX_SCALAR_T residual[JMTX_ARRAY_ATTRIB(const restrict static n)],
    JMTX_SCALAR_T x[JMTX_ARRAY_ATTRIB(const restrict static n)], JMTX_NAME_TYPED(matrix_brm) * r,
    JMTX_SCALAR_T ck[JMTX_ARRAY_ATTRIB(const restrict m)], JMTX_SCALAR_T sk[JMTX_ARRAY_ATTRIB(const restrict m)],
    JMTX_SCALAR_T g[JMTX_ARRAY_ATTRIB(const restrict m)], JMTX_SCALAR_T alpha[JMTX_ARRAY_ATTRIB(const restrict m)],
    JMTX_SCALAR_T h[JMTX_ARRAY_ATTRIB(const restrict m)], JMTX_SCALAR_T p_mat[JMTX_ARRAY_ATTRIB(const restrict m *n)])
{
    JMTX_SCALAR_T *p = p_mat;
    JMTX_INDEX_T n_iteration = 0;
    JMTX_REAL_T err, r_mag = 0;
    for (JMTX_INDEX_T i = 0; i < n; ++i)
    {
        r_mag += JMTX_DOT(residual[i], residual[i]);
    }
    r_mag = JMTX_REAL_ROOT(r_mag);
    err = r_mag / y_mag;
    if (err < tol)
    {
        return 0;
    }

    for (JMTX_INDEX_T i = 0; i < n; ++i)
    {
        p[i] = residual[i] / r_mag;
    }
    g[0] = r_mag;
    JMTX_INDEX_T k;
    for (k = 1; k < m; ++k)
    {
        //  Generate new basis vector
        JMTX_NAME_TYPED(matrix_cds_vector_multiply)(mtx, p, p + n);
        p += n;
        //  Make the basis orthogonal to other basis
        for (JMTX_INDEX_T l = 0; l < k; ++l)
        {
            h[l] = 0;
            const JMTX_SCALAR_T *old_p = p_mat + n * l;
            for (JMTX_INDEX_T i = 0; i < n; ++i)
            {
                h[l] += JMTX_DOT(old_p[i], p[i]);
            }
            for (JMTX_INDEX_T i = 0; i < n; ++i)
            {
                p[i] -= JMTX_DOT(h[l], old_p[i]);
            }
        }

        //  Find magnitude
        JMTX_REAL_T mag_p = 0;
        for (JMTX_INDEX_T i = 0; i < n; ++i)
        {
            mag_p += JMTX_DOT(p[i], p[i]);
        }
        const JMTX_SCALAR_T mag_p2 = mag_p;
        mag_p = JMTX_REAL_ROOT(mag_p);

        //  Normalize basis vector
        for (JMTX_INDEX_T i = 0; i < n; ++i)
        {
            p[i] /= mag_p;
        }

        //  Apply previous Givens rotations to the new column of R
        for (JMTX_INDEX_T l = 0; l < k - 1; ++l)
        {
            const JMTX_SCALAR_T tmp = ck[l] * h[l] + sk[l] * h[l + 1];
            h[l + 1] = -sk[l] * h[l] + ck[l] * h[l + 1];
            h[l] = tmp;
        }

        //  Compute the new givens rotation
        const JMTX_REAL_T rho = JMTX_REAL_ROOT(mag_p2 + h[k - 1] * h[k - 1]);
        const JMTX_SCALAR_T c_new = h[k - 1] / rho;
        const JMTX_SCALAR_T s_new = mag_p / rho;
        ck[k - 1] = c_new;
        sk[k - 1] = s_new;
        h[k - 1] = c_new * h[k - 1] + s_new * mag_p;

        JMTX_NAME_TYPED(matrix_brm_set_col)(r, k - 1, h);

        g[k] = -s_new * g[k - 1];
        g[k - 1] = c_new * g[k - 1];

        r_mag = JMTX_ABS(g[k]);
        err = r_mag / y_mag;
        n_iteration += 1;
        if (err < tol || k + 1 == m)
        {
            break;
        }
    }

    //  Solve the least squares problem using back substitution
    for (JMTX_FAST_INT_T row = 0; row < k; ++row)
    {
        const JMTX_FAST_INT_T i = k - 1 - row;
        JMTX_SCALAR_T *elements;
        JMTX_NAME_TYPED(matrix_brm_get_row)(r, i, &elements);
        JMTX_SCALAR_T sum = 0;
        for (JMTX_FAST_INT_T j = 1; j < row + 1; ++j)
        {
            sum += JMTX_DOT(alpha[i + j], elements[j]);
        }
        alpha[i] = (g[i] - sum) / elements[0];
    }

    //  Compute improvement to x
    for (JMTX_INDEX_T i = 0; i < n; ++i)
    {
        for (JMTX_INDEX_T j = 0; j < k; ++j)
        {
            x[i] += alpha[j] * p_mat[j * n + i];
        }
    }

    return n_iteration;
}

JMTX_INDEX_T JMTX_NAME_TYPED(gmresm_round_crs)(
    const JMTX_NAME_TYPED(matrix_crs) * mtx, const JMTX_INDEX_T n, const JMTX_INDEX_T m, const JMTX_REAL_T y_mag,
    const JMTX_REAL_T tol, const JMTX_SCALAR_T residual[JMTX_ARRAY_ATTRIB(const restrict static n)],
    JMTX_SCALAR_T x[JMTX_ARRAY_ATTRIB(const restrict static n)], JMTX_NAME_TYPED(matrix_brm) * r,
    JMTX_SCALAR_T ck[JMTX_ARRAY_ATTRIB(const restrict m)], JMTX_SCALAR_T sk[JMTX_ARRAY_ATTRIB(const restrict m)],
    JMTX_SCALAR_T g[JMTX_ARRAY_ATTRIB(const restrict m)], JMTX_SCALAR_T alpha[JMTX_ARRAY_ATTRIB(const restrict m)],
    JMTX_SCALAR_T h[JMTX_ARRAY_ATTRIB(const restrict m)], JMTX_SCALAR_T p_mat[JMTX_ARRAY_ATTRIB(const restrict m *n)])
{
    JMTX_SCALAR_T *p = p_mat;
    JMTX_INDEX_T n_iteration = 0;
    JMTX_REAL_T err, r_mag = 0;
    for (JMTX_INDEX_T i = 0; i < n; ++i)
    {
        r_mag += JMTX_DOT(residual[i], residual[i]);
    }
    r_mag = JMTX_REAL_ROOT(r_mag);
    err = r_mag / y_mag;
    if (err < tol)
    {
        return 0;
    }

    for (JMTX_INDEX_T i = 0; i < n; ++i)
    {
        p[i] = residual[i] / r_mag;
    }
    g[0] = r_mag;
    JMTX_INDEX_T k;
    for (k = 1; k < m; ++k)
    {
        //  Generate new basis vector
        JMTX_NAME_TYPED(matrix_crs_vector_multiply)(mtx, p, p + n);
        p += n;
        //  Make the basis orthogonal to other basis
        for (JMTX_INDEX_T l = 0; l < k; ++l)
        {
            h[l] = 0;
            const JMTX_SCALAR_T *old_p = p_mat + n * l;
            for (JMTX_INDEX_T i = 0; i < n; ++i)
            {
                h[l] += JMTX_DOT(old_p[i], p[i]);
            }
            for (JMTX_INDEX_T i = 0; i < n; ++i)
            {
                p[i] -= JMTX_DOT(h[l], old_p[i]);
            }
        }

        //  Find magnitude
        JMTX_REAL_T mag_p = 0;
        for (JMTX_INDEX_T i = 0; i < n; ++i)
        {
            mag_p += JMTX_DOT(p[i], p[i]);
        }
        const JMTX_SCALAR_T mag_p2 = mag_p;
        mag_p = JMTX_REAL_ROOT(mag_p);

        //  Normalize basis vector
        for (JMTX_INDEX_T i = 0; i < n; ++i)
        {
            p[i] /= mag_p;
        }

        //  Apply previous Givens rotations to the new column of R
        for (JMTX_INDEX_T l = 0; l < k - 1; ++l)
        {
            const JMTX_SCALAR_T tmp = ck[l] * h[l] + sk[l] * h[l + 1];
            h[l + 1] = -sk[l] * h[l] + ck[l] * h[l + 1];
            h[l] = tmp;
        }

        //  Compute the new givens rotation
        const JMTX_REAL_T rho = JMTX_REAL_ROOT(mag_p2 + h[k - 1] * h[k - 1]);
        const JMTX_SCALAR_T c_new = h[k - 1] / rho;
        const JMTX_SCALAR_T s_new = mag_p / rho;
        ck[k - 1] = c_new;
        sk[k - 1] = s_new;
        h[k - 1] = c_new * h[k - 1] + s_new * mag_p;

        JMTX_NAME_TYPED(matrix_brm_set_col)(r, k - 1, h);

        g[k] = -s_new * g[k - 1];
        g[k - 1] = c_new * g[k - 1];

        r_mag = JMTX_ABS(g[k]);
        err = r_mag / y_mag;
        n_iteration += 1;
        if (err < tol || k + 1 == m)
        {
            break;
        }
    }

    //  Solve the least squares problem using back substitution
    for (JMTX_FAST_INT_T row = 0; row < k; ++row)
    {
        const JMTX_FAST_INT_T i = k - 1 - row;
        JMTX_SCALAR_T *elements;
        JMTX_NAME_TYPED(matrix_brm_get_row)(r, i, &elements);
        JMTX_SCALAR_T sum = 0;
        for (JMTX_FAST_INT_T j = 1; j < row + 1; ++j)
        {
            sum += JMTX_DOT(alpha[i + j], elements[j]);
        }
        alpha[i] = (g[i] - sum) / elements[0];
    }

    //  Compute improvement to x
    for (JMTX_INDEX_T i = 0; i < n; ++i)
    {
        for (JMTX_INDEX_T j = 0; j < k; ++j)
        {
            x[i] += alpha[j] * p_mat[j * n + i];
        }
    }

    return n_iteration;
}

/**
 * Applies Generalized Minimum Residual method with a restart interval of M (known as GMRES(M)). Builds up a set of m
 * orthonormal basis for the Krylov subspace, then solves a least squares problem to minimize the residual using these
 * basis to solve a problem A x = y.
 *
 *
 * @param mtx system matrix A
 * @param y the solution of the system A x = y
 * @param x the solution vector which contains the initial guess of the solution
 * @param m the GMRES restart interval
 * @param r an m by m upper triangular matrix (lbw = 0, ubw = m - 1) that is to be used in solving the least squares
 * problem
 * @param aux_vec1 auxiliary memory for a vector of m elements
 * @param aux_vec2 auxiliary memory for a vector of m elements
 * @param aux_vec3 auxiliary memory for a vector of m elements
 * @param aux_vec4 auxiliary memory for a vector of m elements
 * @param aux_vec5 auxiliary memory for a vector of m elements
 * @param aux_vecs auxiliary memory for m vectors of the same size as x and y (n by m)
 * @param args::in_convergence_criterion tolerance to determine if the solution is close enough
 * @param args::in_max_iterations number of iterations to stop at
 * @param args::out_last_error receives the value of the error criterion at the final iteration
 * @param args::out_last_iteration receives the number of the final iteration
 * @param args::opt_error_evolution (optional) pointer to an array of length max_iterations, that receives the error
 * value of each iteration
 * @return JMTX_RESULT_SUCCESS if solution converged, JMTX_RESULT_NOT_CONVERGED if solution did not converge in the
 * given number of iterations, other error codes for other errors
 */
jmtx_result JMTX_NAME_TYPED(solve_iterative_gmresm_cds)(
    const JMTX_NAME_TYPED(matrix_cds) * mtx, const JMTX_SCALAR_T *restrict y, JMTX_SCALAR_T *restrict x, JMTX_INDEX_T m,
    JMTX_NAME_TYPED(matrix_brm) * r, JMTX_SCALAR_T aux_vec1[JMTX_ARRAY_ATTRIB(restrict m)],
    JMTX_SCALAR_T aux_vec2[JMTX_ARRAY_ATTRIB(restrict m)], JMTX_SCALAR_T aux_vec3[JMTX_ARRAY_ATTRIB(restrict m)],
    JMTX_SCALAR_T aux_vec4[JMTX_ARRAY_ATTRIB(restrict m)], JMTX_SCALAR_T aux_vec5[JMTX_ARRAY_ATTRIB(restrict m)],
    JMTX_SCALAR_T *restrict aux_vecs, JMTX_NAME_TYPED(solver_arguments) * args)
{
    JMTX_REAL_T err = 0, y_mag = 0, r_mag = 0;
    JMTX_INDEX_T n_iteration = 0;
    const JMTX_INDEX_T n = mtx->base.rows;
    for (JMTX_INDEX_T i = 0; i < n; ++i)
    {
        y_mag += JMTX_DOT(y[i], y[i]);
    }
    y_mag = JMTX_REAL_ROOT(y_mag);

    const JMTX_INDEX_T round_count = args->in_max_iterations / m;
    for (JMTX_INDEX_T round = 0; round < round_count; ++round)
    {
        JMTX_SCALAR_T *const q = aux_vecs;
        JMTX_SCALAR_T *const ck = aux_vec1;
        JMTX_SCALAR_T *const sk = aux_vec2;
        JMTX_SCALAR_T *const g = aux_vec3;
        JMTX_SCALAR_T *const alpha = aux_vec4;
        JMTX_SCALAR_T *const h = aux_vec5;
        JMTX_SCALAR_T *p = q;
        JMTX_NAME_TYPED(matrix_cds_vector_multiply)(mtx, x, p);
        for (JMTX_INDEX_T i = 0; i < n; ++i)
        {
            const JMTX_SCALAR_T res = y[i] - p[i];
            r_mag += JMTX_DOT(res, res);
            p[i] = res;
        }
        r_mag = JMTX_REAL_ROOT(r_mag);
        err = r_mag / y_mag;
        if (err < args->in_convergence_criterion)
        {
            args->out_last_iteration = n_iteration;
            args->out_last_error = err;
            return JMTX_RESULT_SUCCESS;
        }
        for (JMTX_INDEX_T i = 0; i < n; ++i)
        {
            p[i] /= r_mag;
        }
        g[0] = r_mag;
        JMTX_INDEX_T k;
        for (k = 1; k < m; ++k)
        {
            //  Generate new basis vector
            JMTX_NAME_TYPED(matrix_cds_vector_multiply)(mtx, p, p + n);
            p += n;
            //  Make the basis orthogonal to other basis
            for (JMTX_INDEX_T l = 0; l < k; ++l)
            {
                h[l] = 0;
                const JMTX_SCALAR_T *old_p = q + n * l;
                for (JMTX_INDEX_T i = 0; i < n; ++i)
                {
                    h[l] += JMTX_DOT(old_p[i], p[i]);
                }
                for (JMTX_INDEX_T i = 0; i < n; ++i)
                {
                    p[i] -= JMTX_DOT(h[l], old_p[i]);
                }
            }

            //  Find magnitude
            JMTX_REAL_T mag_p = 0;
            for (JMTX_INDEX_T i = 0; i < n; ++i)
            {
                mag_p += JMTX_DOT(p[i], p[i]);
            }
            const JMTX_SCALAR_T mag_p2 = mag_p;
            mag_p = JMTX_REAL_ROOT(mag_p);

            //  Normalize basis vector
            for (JMTX_INDEX_T i = 0; i < n; ++i)
            {
                p[i] /= mag_p;
            }

            //  Apply previous Givens rotations to the new column of R
            for (JMTX_INDEX_T l = 0; l < k - 1; ++l)
            {
                const JMTX_SCALAR_T tmp = ck[l] * h[l] + sk[l] * h[l + 1];
                h[l + 1] = -sk[l] * h[l] + ck[l] * h[l + 1];
                h[l] = tmp;
            }

            //  Compute the new givens rotation
            const JMTX_REAL_T rho = JMTX_REAL_ROOT(mag_p2 + h[k - 1] * h[k - 1]);
            const JMTX_SCALAR_T c_new = h[k - 1] / rho;
            const JMTX_SCALAR_T s_new = mag_p / rho;
            ck[k - 1] = c_new;
            sk[k - 1] = s_new;
            h[k - 1] = c_new * h[k - 1] + s_new * mag_p;

            JMTX_NAME_TYPED(matrix_brm_set_col)(r, k - 1, h);

            g[k] = -s_new * g[k - 1];
            g[k - 1] = c_new * g[k - 1];

            r_mag = JMTX_ABS(g[k]);
            err = r_mag / y_mag;
            if (args->opt_error_evolution)
            {
                args->opt_error_evolution[n_iteration] = err;
            }
            n_iteration += 1;
            if (err < args->in_convergence_criterion || k + 1 == m)
            {
                break;
            }
        }

        //  Solve the least squares problem using back substitution
        for (JMTX_FAST_INT_T row = 0; row < k; ++row)
        {
            const JMTX_FAST_INT_T i = k - 1 - row;
            JMTX_SCALAR_T *elements;
            JMTX_NAME_TYPED(matrix_brm_get_row)(r, i, &elements);
            JMTX_SCALAR_T sum = 0;
            for (JMTX_FAST_INT_T j = 1; j < row + 1; ++j)
            {
                sum += JMTX_DOT(alpha[i + j], elements[j]);
            }
            alpha[i] = (g[i] - sum) / elements[0];
        }

        //  Compute improvement to x
        for (JMTX_INDEX_T i = 0; i < n; ++i)
        {
            for (JMTX_INDEX_T j = 0; j < k; ++j)
            {
                x[i] += alpha[j] * q[j * n + i];
            }
        }
        if (err < args->in_convergence_criterion)
        {
            args->out_last_error = err;
            args->out_last_iteration = n_iteration;
            return JMTX_RESULT_SUCCESS;
        }
    }
    args->out_last_error = err;
    args->out_last_iteration = n_iteration;
    return JMTX_RESULT_NOT_CONVERGED;
}

/**
 * Applies Generalized Minimum Residual method with a restart interval of M (known as GMRES(M)). Builds up a set of m
 * orthonormal basis for the Krylov subspace, then solves a least squares problem to minimize the residual using these
 * basis to solve a problem A x = y. Using a preconditioner means that the Krylov subspace is not based on the system
 * matrix A, but instead on the preconditioned matrix.
 *
 * Uses Right Preconditioning with the Jacobi iteration, meaning it uses the it actually solves a different system:
 *                                      A D⁻¹ z = y, then solves z = D x
 * This is done in hopes of A D⁻¹ having a lower condition number than A. Right preconditioning retains the exact
 * same residual as the original system A x = y, which may or may not be desired.
 *
 *
 * @param mtx system matrix A
 * @param y the solution of the system A x = y
 * @param x the solution vector which contains the initial guess of the solution
 * @param m the GMRES restart interval
 * @param r an m by m upper triangular matrix (lbw = 0, ubw = m - 1) that is to be used in solving the least squares
 * problem
 * @param aux_vec1 auxiliary memory for a vector of m elements
 * @param aux_vec2 auxiliary memory for a vector of m elements
 * @param aux_vec3 auxiliary memory for a vector of m elements
 * @param aux_vec4 auxiliary memory for a vector of m elements
 * @param aux_vec5 auxiliary memory for a vector of m elements
 * @param aux_vec6 auxiliary memory for a vector of of the same size as x and y
 * @param aux_vec7 auxiliary memory for a vector of of the same size as x and y
 * @param aux_vecs auxiliary memory for m vectors of the same size as x and y (n by m)
 * @param args::in_convergence_criterion tolerance to determine if the solution is close enough
 * @param args::in_max_iterations number of iterations to stop at
 * @param args::out_last_error receives the value of the error criterion at the final iteration
 * @param args::out_last_iteration receives the number of the final iteration
 * @param args::opt_error_evolution (optional) pointer to an array of length max_iterations, that receives the error
 * value of each iteration
 * @return JMTX_RESULT_SUCCESS if solution converged, JMTX_RESULT_NOT_CONVERGED if solution did not converge in the
 * given number of iterations, other error codes for other errors
 */
jmtx_result JMTX_NAME_TYPED(solve_iterative_gmresm_rpc_jacobi_cds)(
    const JMTX_NAME_TYPED(matrix_cds) * mtx, const JMTX_SCALAR_T *restrict y, JMTX_SCALAR_T *restrict x, JMTX_INDEX_T m,
    JMTX_NAME_TYPED(matrix_brm) * r, JMTX_SCALAR_T aux_vec1[JMTX_ARRAY_ATTRIB(restrict m)],
    JMTX_SCALAR_T aux_vec2[JMTX_ARRAY_ATTRIB(restrict m)], JMTX_SCALAR_T aux_vec3[JMTX_ARRAY_ATTRIB(restrict m)],
    JMTX_SCALAR_T aux_vec4[JMTX_ARRAY_ATTRIB(restrict m)], JMTX_SCALAR_T aux_vec5[JMTX_ARRAY_ATTRIB(restrict m)],
    JMTX_SCALAR_T *restrict aux_vec6, JMTX_SCALAR_T *restrict aux_vec7, JMTX_SCALAR_T *restrict aux_vecs,
    JMTX_NAME_TYPED(solver_arguments) * args)
{
    JMTX_REAL_T err = 0, y_mag = 0, r_mag = 0;
    JMTX_INDEX_T n_iteration = 0;
    JMTX_SCALAR_T *const d_inv = aux_vec6;
    const JMTX_INDEX_T n = mtx->base.rows;
    for (JMTX_INDEX_T i = 0; i < n; ++i)
    {
        d_inv[i] = 1 / mtx->main_diagonal[i];
        y_mag += JMTX_DOT(y[i], y[i]);
    }
    y_mag = JMTX_REAL_ROOT(y_mag);

    const JMTX_INDEX_T round_count = args->in_max_iterations / m;
    for (JMTX_INDEX_T round = 0; round < round_count; ++round)
    {
        JMTX_SCALAR_T *const q = aux_vecs;
        JMTX_SCALAR_T *const ck = aux_vec1;
        JMTX_SCALAR_T *const sk = aux_vec2;
        JMTX_SCALAR_T *const g = aux_vec3;
        JMTX_SCALAR_T *const alpha = aux_vec4;
        JMTX_SCALAR_T *const h = aux_vec5;
        JMTX_SCALAR_T *p = q;

        JMTX_NAME_TYPED(matrix_cds_vector_multiply)(mtx, x, p);
        for (JMTX_INDEX_T i = 0; i < n; ++i)
        {
            const JMTX_SCALAR_T res = y[i] - p[i];
            r_mag += JMTX_DOT(res, res);
            p[i] = res;
        }
        r_mag = JMTX_REAL_ROOT(r_mag);
        err = r_mag / y_mag;
        if (err < args->in_convergence_criterion)
        {
            args->out_last_iteration = n_iteration;
            args->out_last_error = err;
            return JMTX_RESULT_SUCCESS;
        }
        for (JMTX_INDEX_T i = 0; i < n; ++i)
        {
            p[i] /= r_mag;
        }
        g[0] = r_mag;
        JMTX_INDEX_T k;
        for (k = 1; k < m; ++k)
        {
            //  Generate new basis vector
            for (JMTX_INDEX_T i = 0; i < n; ++i)
            {
                aux_vec7[i] = p[i] * d_inv[i];
            }
            JMTX_NAME_TYPED(matrix_cds_vector_multiply)(mtx, aux_vec7, p + n);
            p += n;
            //  Make the basis orthogonal to other basis
            for (JMTX_INDEX_T l = 0; l < k; ++l)
            {
                h[l] = 0;
                const JMTX_SCALAR_T *old_p = q + n * l;
                for (JMTX_INDEX_T i = 0; i < n; ++i)
                {
                    h[l] += JMTX_DOT(old_p[i], p[i]);
                }
                for (JMTX_INDEX_T i = 0; i < n; ++i)
                {
                    p[i] -= JMTX_DOT(h[l], old_p[i]);
                }
            }

            //  Find magnitude
            JMTX_REAL_T mag_p = 0;
            for (JMTX_INDEX_T i = 0; i < n; ++i)
            {
                mag_p += JMTX_DOT(p[i], p[i]);
            }
            const JMTX_SCALAR_T mag_p2 = mag_p;
            mag_p = JMTX_REAL_ROOT(mag_p);

            //  Normalize basis vector
            for (JMTX_INDEX_T i = 0; i < n; ++i)
            {
                p[i] /= mag_p;
            }

            //  Apply previous Givens rotations to the new column of R
            for (JMTX_INDEX_T l = 0; l < k - 1; ++l)
            {
                const JMTX_SCALAR_T tmp = ck[l] * h[l] + sk[l] * h[l + 1];
                h[l + 1] = -sk[l] * h[l] + ck[l] * h[l + 1];
                h[l] = tmp;
            }

            //  Compute the new givens rotation
            const JMTX_REAL_T rho = JMTX_REAL_ROOT(mag_p2 + h[k - 1] * h[k - 1]);
            const JMTX_SCALAR_T c_new = h[k - 1] / rho;
            const JMTX_SCALAR_T s_new = mag_p / rho;
            ck[k - 1] = c_new;
            sk[k - 1] = s_new;
            h[k - 1] = c_new * h[k - 1] + s_new * mag_p;

            JMTX_NAME_TYPED(matrix_brm_set_col)(r, k - 1, h);

            g[k] = -s_new * g[k - 1];
            g[k - 1] = c_new * g[k - 1];

            r_mag = JMTX_ABS(g[k]);
            err = r_mag / y_mag;
            if (args->opt_error_evolution)
            {
                args->opt_error_evolution[n_iteration] = err;
            }
            n_iteration += 1;
            if (err < args->in_convergence_criterion || k + 1 == m)
            {
                break;
            }
        }

        //  Solve the least squares problem using back substitution
        for (JMTX_FAST_INT_T row = 0; row < k; ++row)
        {
            const JMTX_FAST_INT_T i = k - 1 - row;
            JMTX_SCALAR_T *elements;
            JMTX_NAME_TYPED(matrix_brm_get_row)(r, i, &elements);
            JMTX_SCALAR_T sum = 0;
            for (JMTX_FAST_INT_T j = 1; j < row + 1; ++j)
            {
                sum += JMTX_DOT(alpha[i + j], elements[j]);
            }
            alpha[i] = (g[i] - sum) / elements[0];
        }

        //  Compute improvement to x
        //      Apply preconditioner on x first
        for (JMTX_INDEX_T i = 0; i < n; ++i)
        {
            x[i] *= d_inv[i];
        }
        for (JMTX_INDEX_T i = 0; i < n; ++i)
        {
            for (JMTX_INDEX_T j = 0; j < k; ++j)
            {
                x[i] += alpha[j] * q[j * n + i];
            }
        }
        //      Undo the preconditioner
        for (JMTX_INDEX_T i = 0; i < n; ++i)
        {
            x[i] *= mtx->main_diagonal[i];
        }

        if (err < args->in_convergence_criterion)
        {
            args->out_last_error = err;
            args->out_last_iteration = n_iteration;
            return JMTX_RESULT_SUCCESS;
        }
    }
    args->out_last_error = err;
    args->out_last_iteration = n_iteration;
    return JMTX_RESULT_NOT_CONVERGED;
}

/**
 * Applies Generalized Minimum Residual method with a restart interval of M (known as GMRES(M)). Builds up a set of m
 * orthonormal basis for the Krylov subspace, then solves a least squares problem to minimize the residual using these
 * basis to solve a problem A x = y. Using a preconditioner means that the Krylov subspace is not based on the system
 * matrix A, but instead on the preconditioned matrix.
 *
 * Uses Left Preconditioning with the Jacobi iteration, meaning it uses the it actually solves a different system:
 *                                            D⁻¹ A x = D⁻¹ y
 * This is done in hopes of D⁻¹ A having a lower condition number than A. Left preconditioning has a the consequence of
 * not actually using/minimizing the real residual but instead the residual of the preconditioned system. This may or
 * may not be desired.
 *
 *
 * @param mtx system matrix A
 * @param y the solution of the system A x = y
 * @param x the solution vector which contains the initial guess of the solution
 * @param m the GMRES restart interval
 * @param r an m by m upper triangular matrix (lbw = 0, ubw = m - 1) that is to be used in solving the least squares
 * problem
 * @param aux_vec1 auxiliary memory for a vector of m elements
 * @param aux_vec2 auxiliary memory for a vector of m elements
 * @param aux_vec3 auxiliary memory for a vector of m elements
 * @param aux_vec4 auxiliary memory for a vector of m elements
 * @param aux_vec5 auxiliary memory for a vector of m elements
 * @param aux_vec6 auxiliary memory for a vector of of the same size as x and y
 * @param aux_vecs auxiliary memory for m vectors of the same size as x and y (n by m)
 * @param args::in_convergence_criterion tolerance to determine if the solution is close enough
 * @param args::in_max_iterations number of iterations to stop at
 * @param args::out_last_error receives the value of the error criterion at the final iteration
 * @param args::out_last_iteration receives the number of the final iteration
 * @param args::opt_error_evolution (optional) pointer to an array of length max_iterations, that receives the error
 * value of each iteration
 * @return JMTX_RESULT_SUCCESS if solution converged, JMTX_RESULT_NOT_CONVERGED if solution did not converge in the
 * given number of iterations, other error codes for other errors
 */
jmtx_result JMTX_NAME_TYPED(solve_iterative_gmresm_lpc_jacobi_cds)(
    const JMTX_NAME_TYPED(matrix_cds) * mtx, const JMTX_SCALAR_T *restrict y, JMTX_SCALAR_T *restrict x, JMTX_INDEX_T m,
    JMTX_NAME_TYPED(matrix_brm) * r, JMTX_SCALAR_T aux_vec1[JMTX_ARRAY_ATTRIB(restrict m)],
    JMTX_SCALAR_T aux_vec2[JMTX_ARRAY_ATTRIB(restrict m)], JMTX_SCALAR_T aux_vec3[JMTX_ARRAY_ATTRIB(restrict m)],
    JMTX_SCALAR_T aux_vec4[JMTX_ARRAY_ATTRIB(restrict m)], JMTX_SCALAR_T aux_vec5[JMTX_ARRAY_ATTRIB(restrict m)],
    JMTX_SCALAR_T *restrict aux_vec6, JMTX_SCALAR_T *restrict aux_vecs, JMTX_NAME_TYPED(solver_arguments) * args)
{
    JMTX_REAL_T err = 0, y_mag = 0, r_mag = 0;
    JMTX_INDEX_T n_iteration = 0;
    JMTX_SCALAR_T *const d_inv = aux_vec6;
    const JMTX_INDEX_T n = mtx->base.rows;
    for (JMTX_INDEX_T i = 0; i < n; ++i)
    {
        d_inv[i] = 1 / mtx->main_diagonal[i];
        const JMTX_SCALAR_T scaled_y = y[i] * d_inv[i];
        y_mag += scaled_y * scaled_y;
    }
    y_mag = JMTX_REAL_ROOT(y_mag);

    const JMTX_INDEX_T round_count = args->in_max_iterations / m;
    for (JMTX_INDEX_T round = 0; round < round_count; ++round)
    {
        JMTX_SCALAR_T *const q = aux_vecs;
        JMTX_SCALAR_T *const ck = aux_vec1;
        JMTX_SCALAR_T *const sk = aux_vec2;
        JMTX_SCALAR_T *const g = aux_vec3;
        JMTX_SCALAR_T *const alpha = aux_vec4;
        JMTX_SCALAR_T *const h = aux_vec5;
        JMTX_SCALAR_T *p = q;

        JMTX_NAME_TYPED(matrix_cds_vector_multiply)(mtx, x, p);
        for (JMTX_INDEX_T i = 0; i < n; ++i)
        {
            x[i] *= d_inv[i];
        }
        for (JMTX_INDEX_T i = 0; i < n; ++i)
        {
            const JMTX_SCALAR_T res = y[i] * d_inv[i] - p[i];
            r_mag += JMTX_DOT(res, res);
            p[i] = res;
        }
        r_mag = JMTX_REAL_ROOT(r_mag);
        err = r_mag / y_mag;
        if (err < args->in_convergence_criterion)
        {
            args->out_last_iteration = n_iteration;
            args->out_last_error = err;
            return JMTX_RESULT_SUCCESS;
        }
        for (JMTX_INDEX_T i = 0; i < n; ++i)
        {
            p[i] /= r_mag;
        }
        g[0] = r_mag;
        JMTX_INDEX_T k;
        for (k = 1; k < m; ++k)
        {
            //  Generate new basis vector
            JMTX_NAME_TYPED(matrix_cds_vector_multiply)(mtx, p, p + n);
            p += n;
            for (JMTX_INDEX_T i = 0; i < n; ++i)
            {
                p[i] *= d_inv[i];
            }
            //  Make the basis orthogonal to other basis
            for (JMTX_INDEX_T l = 0; l < k; ++l)
            {
                h[l] = 0;
                const JMTX_SCALAR_T *old_p = q + n * l;
                for (JMTX_INDEX_T i = 0; i < n; ++i)
                {
                    h[l] += JMTX_DOT(old_p[i], p[i]);
                }
                for (JMTX_INDEX_T i = 0; i < n; ++i)
                {
                    p[i] -= JMTX_DOT(h[l], old_p[i]);
                }
            }

            //  Find magnitude
            JMTX_REAL_T mag_p = 0;
            for (JMTX_INDEX_T i = 0; i < n; ++i)
            {
                mag_p += JMTX_DOT(p[i], p[i]);
            }
            const JMTX_SCALAR_T mag_p2 = mag_p;
            mag_p = JMTX_REAL_ROOT(mag_p);

            //  Normalize basis vector
            for (JMTX_INDEX_T i = 0; i < n; ++i)
            {
                p[i] /= mag_p;
            }

            //  Apply previous Givens rotations to the new column of R
            for (JMTX_INDEX_T l = 0; l < k - 1; ++l)
            {
                const JMTX_SCALAR_T tmp = ck[l] * h[l] + sk[l] * h[l + 1];
                h[l + 1] = -sk[l] * h[l] + ck[l] * h[l + 1];
                h[l] = tmp;
            }

            //  Compute the new givens rotation
            const JMTX_REAL_T rho = JMTX_REAL_ROOT(mag_p2 + h[k - 1] * h[k - 1]);
            const JMTX_SCALAR_T c_new = h[k - 1] / rho;
            const JMTX_SCALAR_T s_new = mag_p / rho;
            ck[k - 1] = c_new;
            sk[k - 1] = s_new;
            h[k - 1] = c_new * h[k - 1] + s_new * mag_p;

            JMTX_NAME_TYPED(matrix_brm_set_col)(r, k - 1, h);

            g[k] = -s_new * g[k - 1];
            g[k - 1] = c_new * g[k - 1];

            r_mag = JMTX_ABS(g[k]);
            err = r_mag / y_mag;
            if (args->opt_error_evolution)
            {
                args->opt_error_evolution[n_iteration] = err;
            }
            n_iteration += 1;
            if (err < args->in_convergence_criterion || k + 1 == m)
            {
                break;
            }
        }

        //  Solve the least squares problem using back substitution
        for (JMTX_FAST_INT_T row = 0; row < k; ++row)
        {
            const JMTX_FAST_INT_T i = k - 1 - row;
            JMTX_SCALAR_T *elements;
            JMTX_NAME_TYPED(matrix_brm_get_row)(r, i, &elements);
            JMTX_SCALAR_T sum = 0;
            for (JMTX_FAST_INT_T j = 1; j < row + 1; ++j)
            {
                sum += JMTX_DOT(alpha[i + j], elements[j]);
            }
            alpha[i] = (g[i] - sum) / elements[0];
        }

        //  Compute improvement to x
        for (JMTX_INDEX_T i = 0; i < n; ++i)
        {
            for (JMTX_INDEX_T j = 0; j < k; ++j)
            {
                x[i] += alpha[j] * q[j * n + i];
            }
        }

        if (err < args->in_convergence_criterion)
        {
            args->out_last_error = err;
            args->out_last_iteration = n_iteration;
            return JMTX_RESULT_SUCCESS;
        }
    }
    args->out_last_error = err;
    args->out_last_iteration = n_iteration;
    return JMTX_RESULT_NOT_CONVERGED;
}
