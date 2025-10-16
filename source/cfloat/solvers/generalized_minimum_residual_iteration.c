//
// Created by jan on 1.1.2024.
//

#include "../../../include/jmtx/cfloat/solvers/generalized_minimum_residual_iteration.h"
#include "../matrices/band_row_major_internal.h"
#include "../matrices/sparse_diagonal_compressed_internal.h"
#include "../matrices/sparse_row_compressed_internal.h"
#include <assert.h>
#include <complex.h>
#include <math.h>

jmtx_result jmtxc_solve_iterative_gmresm_crs(const jmtxc_matrix_crs *mtx, const _Complex float *restrict y,
                                             _Complex float *restrict x, uint32_t m, jmtxc_matrix_brm *r,
                                             _Complex float aux_vec1[JMTX_ARRAY_ATTRIB(restrict m)],
                                             _Complex float aux_vec2[JMTX_ARRAY_ATTRIB(restrict m)],
                                             _Complex float aux_vec3[JMTX_ARRAY_ATTRIB(restrict m)],
                                             _Complex float aux_vec4[JMTX_ARRAY_ATTRIB(restrict m)],
                                             _Complex float aux_vec5[JMTX_ARRAY_ATTRIB(restrict m)],
                                             _Complex float *restrict aux_vecs, jmtxf_solver_arguments *args)
{
    float err = 0, y_mag = 0, r_mag = 0;
    uint32_t n_iteration = 0;
    const uint32_t n = mtx->base.rows;
    for (uint32_t i = 0; i < n; ++i)
    {
        y_mag += conjf(y[i]) * y[i];
    }
    y_mag = sqrtf(y_mag);

    const uint32_t round_count = args->in_max_iterations / m;
    for (uint32_t round = 0; round < round_count; ++round)
    {
        _Complex float *const q = aux_vecs;
        _Complex float *const ck = aux_vec1;
        _Complex float *const sk = aux_vec2;
        _Complex float *const g = aux_vec3;
        _Complex float *const alpha = aux_vec4;
        _Complex float *const h = aux_vec5;
        _Complex float *p = q;
        for (uint32_t i = 0; i < n; ++i)
        {
            const _Complex float res = y[i] - jmtxc_matrix_crs_vector_multiply_row(mtx, x, i);
            r_mag += conjf(res) * res;
            p[i] = res;
        }
        r_mag = sqrtf(r_mag);
        err = r_mag / y_mag;
        if (err < args->in_convergence_criterion)
        {
            args->out_last_iteration = n_iteration;
            args->out_last_error = err;
            return JMTX_RESULT_SUCCESS;
        }
        for (uint32_t i = 0; i < n; ++i)
        {
            p[i] /= r_mag;
        }
        g[0] = r_mag;
        uint32_t k;
        for (k = 1; k < m; ++k)
        {
            //  Generate new basis vector
            jmtxc_matrix_crs_vector_multiply(mtx, p, p + n);
            p += n;
            //  Make the basis orthogonal to other basis
            for (uint32_t l = 0; l < k; ++l)
            {
                h[l] = 0;
                const _Complex float *old_p = q + n * l;
                for (uint32_t i = 0; i < n; ++i)
                {
                    h[l] += conjf(old_p[i]) * p[i];
                }
                for (uint32_t i = 0; i < n; ++i)
                {
                    p[i] -= h[l] * old_p[i];
                }
            }

            //  Find magnitude
            float mag_p = 0;
            for (uint32_t i = 0; i < n; ++i)
            {
                mag_p += conjf(p[i]) * p[i];
            }
            const float mag_p2 = mag_p;
            mag_p = sqrtf(mag_p);

            //  Normalize basis vector
            for (uint32_t i = 0; i < n; ++i)
            {
                p[i] /= mag_p;
            }

            //  Apply previous Givens rotations to the new column of R
            for (uint32_t l = 0; l < k - 1; ++l)
            {
                const _Complex float tmp = ck[l] * h[l] + sk[l] * h[l + 1];
                h[l + 1] = -sk[l] * h[l] + ck[l] * h[l + 1];
                h[l] = tmp;
            }

            //  Compute the new givens rotation
            const float rho = sqrtf(mag_p2 + conjf(h[k - 1]) * h[k - 1]);
            const _Complex float c_new = h[k - 1] / rho;
            const _Complex float s_new = mag_p / rho;
            ck[k - 1] = c_new;
            sk[k - 1] = s_new;
            h[k - 1] = c_new * h[k - 1] + s_new * mag_p;

            jmtxc_matrix_brm_set_col(r, k - 1, h);

            g[k] = -s_new * g[k - 1];
            g[k - 1] = c_new * g[k - 1];

            r_mag = cabsf(g[k]);
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
        for (uint_fast32_t row = 0; row < k; ++row)
        {
            const uint_fast32_t i = k - 1 - row;
            _Complex float *elements;
            jmtxc_matrix_brm_get_row(r, i, &elements);
            _Complex float sum = 0;
            for (uint_fast32_t j = 1; j < row + 1; ++j)
            {
                sum += alpha[i + j] * elements[j];
            }
            alpha[i] = (g[i] - sum) / elements[0];
        }

        //  Compute improvement to x
        for (uint32_t i = 0; i < n; ++i)
        {
            for (uint32_t j = 0; j < k; ++j)
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
 * basis to solve a problem A x = y.
 *
 *
 * @param mtx system matrix A
 * @param n size of the problem
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
jmtx_result jmtxcs_solve_iterative_gmresm_crs(
    const jmtxc_matrix_crs *mtx, uint32_t n, const _Complex float y[JMTX_ARRAY_ATTRIB(restrict static n)],
    _Complex float x[JMTX_ARRAY_ATTRIB(restrict static n)], uint32_t m, jmtxc_matrix_brm *r,
    _Complex float aux_vec1[JMTX_ARRAY_ATTRIB(restrict m)], _Complex float aux_vec2[JMTX_ARRAY_ATTRIB(restrict m)],
    _Complex float aux_vec3[JMTX_ARRAY_ATTRIB(restrict m)], _Complex float aux_vec4[JMTX_ARRAY_ATTRIB(restrict m)],
    _Complex float aux_vec5[JMTX_ARRAY_ATTRIB(restrict m)], _Complex float aux_vecs[JMTX_ARRAY_ATTRIB(restrict m * n)],
    jmtxf_solver_arguments *args)
{
    if (mtx->base.type != JMTXC_TYPE_CRS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (mtx->base.rows != n || mtx->base.cols != n)
    {
        return JMTX_RESULT_BAD_MATRIX;
    }
    if (r->base.type != JMTXC_TYPE_BRM)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (r->base.rows != m || r->base.cols != m)
    {
        return JMTX_RESULT_BAD_MATRIX;
    }
    if (r->upper_bandwidth != m - 1 || r->lower_bandwidth != 0)
    {
        return JMTX_RESULT_BAD_MATRIX;
    }
    if (m == 0)
    {
        return JMTX_RESULT_BAD_PARAM;
    }
    return jmtxc_solve_iterative_gmresm_crs(mtx, y, x, m, r, aux_vec1, aux_vec2, aux_vec3, aux_vec4, aux_vec5, aux_vecs,
                                            args);
}

uint32_t jmtxc_gmresm_round_cds(
    const jmtxc_matrix_cds *mtx, const uint32_t n, const uint32_t m, const _Complex float y_mag, const float tol,
    const _Complex float residual[JMTX_ARRAY_ATTRIB(const restrict static n)],
    _Complex float x[JMTX_ARRAY_ATTRIB(const restrict static n)], jmtxc_matrix_brm *r,
    _Complex float ck[JMTX_ARRAY_ATTRIB(const restrict m)], _Complex float sk[JMTX_ARRAY_ATTRIB(const restrict m)],
    _Complex float g[JMTX_ARRAY_ATTRIB(const restrict m)], _Complex float alpha[JMTX_ARRAY_ATTRIB(const restrict m)],
    _Complex float h[JMTX_ARRAY_ATTRIB(const restrict m)], _Complex float p_mat[JMTX_ARRAY_ATTRIB(const restrict m *n)])
{
    _Complex float *p = p_mat;
    uint32_t n_iteration = 0;
    float err, r_mag = 0;
    for (uint32_t i = 0; i < n; ++i)
    {
        r_mag += conjf(residual[i]) * residual[i];
    }
    r_mag = sqrtf(r_mag);
    err = r_mag / y_mag;
    if (err < tol)
    {
        return 0;
    }

    for (uint32_t i = 0; i < n; ++i)
    {
        p[i] = residual[i] / r_mag;
    }
    g[0] = r_mag;
    uint32_t k;
    for (k = 1; k < m; ++k)
    {
        //  Generate new basis vector
        jmtxc_matrix_cds_vector_multiply(mtx, p, p + n);
        p += n;
        //  Make the basis orthogonal to other basis
        for (uint32_t l = 0; l < k; ++l)
        {
            h[l] = 0;
            const _Complex float *old_p = p_mat + n * l;
            for (uint32_t i = 0; i < n; ++i)
            {
                h[l] += conjf(old_p[i]) * p[i];
            }
            for (uint32_t i = 0; i < n; ++i)
            {
                p[i] -= h[l] * old_p[i];
            }
        }

        //  Find magnitude
        float mag_p = 0;
        for (uint32_t i = 0; i < n; ++i)
        {
            mag_p += conjf(p[i]) * p[i];
        }
        const float mag_p2 = mag_p;
        mag_p = sqrtf(mag_p);

        //  Normalize basis vector
        for (uint32_t i = 0; i < n; ++i)
        {
            p[i] /= mag_p;
        }

        //  Apply previous Givens rotations to the new column of R
        for (uint32_t l = 0; l < k - 1; ++l)
        {
            const _Complex float tmp = ck[l] * h[l] + sk[l] * h[l + 1];
            h[l + 1] = -sk[l] * h[l] + ck[l] * h[l + 1];
            h[l] = tmp;
        }

        //  Compute the new givens rotation
        const float rho = sqrtf(mag_p2 + conjf(h[k - 1]) * h[k - 1]);
        const _Complex float c_new = h[k - 1] / rho;
        const _Complex float s_new = mag_p / rho;
        ck[k - 1] = c_new;
        sk[k - 1] = s_new;
        h[k - 1] = c_new * h[k - 1] + s_new * mag_p;

        jmtxc_matrix_brm_set_col(r, k - 1, h);

        g[k] = -s_new * g[k - 1];
        g[k - 1] = c_new * g[k - 1];

        r_mag = cabsf(g[k]);
        err = r_mag / y_mag;
        n_iteration += 1;
        if (err < tol || k + 1 == m)
        {
            break;
        }
    }

    //  Solve the least squares problem using back substitution
    for (uint_fast32_t row = 0; row < k; ++row)
    {
        const uint_fast32_t i = k - 1 - row;
        _Complex float *elements;
        jmtxc_matrix_brm_get_row(r, i, &elements);
        _Complex float sum = 0;
        for (uint_fast32_t j = 1; j < row + 1; ++j)
        {
            sum += alpha[i + j] * elements[j];
        }
        alpha[i] = (g[i] - sum) / elements[0];
    }

    //  Compute improvement to x
    for (uint32_t i = 0; i < n; ++i)
    {
        for (uint32_t j = 0; j < k; ++j)
        {
            x[i] += alpha[j] * p_mat[j * n + i];
        }
    }

    return n_iteration;
}

uint32_t jmtxc_gmresm_round_crs(
    const jmtxc_matrix_crs *mtx, const uint32_t n, const uint32_t m, const _Complex float y_mag, const float tol,
    const _Complex float residual[JMTX_ARRAY_ATTRIB(const restrict static n)],
    _Complex float x[JMTX_ARRAY_ATTRIB(const restrict static n)], jmtxc_matrix_brm *r,
    _Complex float ck[JMTX_ARRAY_ATTRIB(const restrict m)], _Complex float sk[JMTX_ARRAY_ATTRIB(const restrict m)],
    _Complex float g[JMTX_ARRAY_ATTRIB(const restrict m)], _Complex float alpha[JMTX_ARRAY_ATTRIB(const restrict m)],
    _Complex float h[JMTX_ARRAY_ATTRIB(const restrict m)], _Complex float p_mat[JMTX_ARRAY_ATTRIB(const restrict m *n)])
{
    _Complex float *p = p_mat;
    uint32_t n_iteration = 0;
    float err, r_mag = 0;
    for (uint32_t i = 0; i < n; ++i)
    {
        r_mag += conjf(residual[i]) * residual[i];
    }
    r_mag = sqrtf(r_mag);
    err = r_mag / y_mag;
    if (err < tol)
    {
        return 0;
    }

    for (uint32_t i = 0; i < n; ++i)
    {
        p[i] = residual[i] / r_mag;
    }
    g[0] = r_mag;
    uint32_t k;
    for (k = 1; k < m; ++k)
    {
        //  Generate new basis vector
        jmtxc_matrix_crs_vector_multiply(mtx, p, p + n);
        p += n;
        //  Make the basis orthogonal to other basis
        for (uint32_t l = 0; l < k; ++l)
        {
            h[l] = 0;
            const _Complex float *old_p = p_mat + n * l;
            for (uint32_t i = 0; i < n; ++i)
            {
                h[l] += conjf(old_p[i]) * p[i];
            }
            for (uint32_t i = 0; i < n; ++i)
            {
                p[i] -= h[l] * old_p[i];
            }
        }

        //  Find magnitude
        float mag_p = 0;
        for (uint32_t i = 0; i < n; ++i)
        {
            mag_p += conjf(p[i]) * p[i];
        }
        const float mag_p2 = mag_p;
        mag_p = sqrtf(mag_p);

        //  Normalize basis vector
        for (uint32_t i = 0; i < n; ++i)
        {
            p[i] /= mag_p;
        }

        //  Apply previous Givens rotations to the new column of R
        for (uint32_t l = 0; l < k - 1; ++l)
        {
            const _Complex float tmp = ck[l] * h[l] + sk[l] * h[l + 1];
            h[l + 1] = -sk[l] * h[l] + ck[l] * h[l + 1];
            h[l] = tmp;
        }

        //  Compute the new givens rotation
        const float rho = sqrtf(mag_p2 + conjf(h[k - 1]) * h[k - 1]);
        const _Complex float c_new = h[k - 1] / rho;
        const _Complex float s_new = mag_p / rho;
        ck[k - 1] = c_new;
        sk[k - 1] = s_new;
        h[k - 1] = c_new * h[k - 1] + s_new * mag_p;

        jmtxc_matrix_brm_set_col(r, k - 1, h);

        g[k] = -s_new * g[k - 1];
        g[k - 1] = c_new * g[k - 1];

        r_mag = cabsf(g[k]);
        err = r_mag / y_mag;
        n_iteration += 1;
        if (err < tol || k + 1 == m)
        {
            break;
        }
    }

    //  Solve the least squares problem using back substitution
    for (uint_fast32_t row = 0; row < k; ++row)
    {
        const uint_fast32_t i = k - 1 - row;
        _Complex float *elements;
        jmtxc_matrix_brm_get_row(r, i, &elements);
        _Complex float sum = 0;
        for (uint_fast32_t j = 1; j < row + 1; ++j)
        {
            sum += alpha[i + j] * elements[j];
        }
        alpha[i] = (g[i] - sum) / elements[0];
    }

    //  Compute improvement to x
    for (uint32_t i = 0; i < n; ++i)
    {
        for (uint32_t j = 0; j < k; ++j)
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
jmtx_result jmtxc_solve_iterative_gmresm_cds(const jmtxc_matrix_cds *mtx, const _Complex float *restrict y,
                                             _Complex float *restrict x, uint32_t m, jmtxc_matrix_brm *r,
                                             _Complex float aux_vec1[JMTX_ARRAY_ATTRIB(restrict m)],
                                             _Complex float aux_vec2[JMTX_ARRAY_ATTRIB(restrict m)],
                                             _Complex float aux_vec3[JMTX_ARRAY_ATTRIB(restrict m)],
                                             _Complex float aux_vec4[JMTX_ARRAY_ATTRIB(restrict m)],
                                             _Complex float aux_vec5[JMTX_ARRAY_ATTRIB(restrict m)],
                                             _Complex float *restrict aux_vecs, jmtxf_solver_arguments *args)
{
    float err = 0, y_mag = 0, r_mag = 0;
    uint32_t n_iteration = 0;
    const uint32_t n = mtx->base.rows;
    for (uint32_t i = 0; i < n; ++i)
    {
        y_mag += conjf(y[i]) * y[i];
    }
    y_mag = sqrtf(y_mag);

    const uint32_t round_count = args->in_max_iterations / m;
    for (uint32_t round = 0; round < round_count; ++round)
    {
        _Complex float *const q = aux_vecs;
        _Complex float *const ck = aux_vec1;
        _Complex float *const sk = aux_vec2;
        _Complex float *const g = aux_vec3;
        _Complex float *const alpha = aux_vec4;
        _Complex float *const h = aux_vec5;
        _Complex float *p = q;
        jmtxc_matrix_cds_vector_multiply(mtx, x, p);
        for (uint32_t i = 0; i < n; ++i)
        {
            const _Complex float res = y[i] - p[i];
            r_mag += conjf(res) * res;
            p[i] = res;
        }
        r_mag = sqrtf(r_mag);
        err = r_mag / y_mag;
        if (err < args->in_convergence_criterion)
        {
            args->out_last_iteration = n_iteration;
            args->out_last_error = err;
            return JMTX_RESULT_SUCCESS;
        }
        for (uint32_t i = 0; i < n; ++i)
        {
            p[i] /= r_mag;
        }
        g[0] = r_mag;
        uint32_t k;
        for (k = 1; k < m; ++k)
        {
            //  Generate new basis vector
            jmtxc_matrix_cds_vector_multiply(mtx, p, p + n);
            p += n;
            //  Make the basis orthogonal to other basis
            for (uint32_t l = 0; l < k; ++l)
            {
                h[l] = 0;
                const _Complex float *old_p = q + n * l;
                for (uint32_t i = 0; i < n; ++i)
                {
                    h[l] += conjf(old_p[i]) * p[i];
                }
                for (uint32_t i = 0; i < n; ++i)
                {
                    p[i] -= h[l] * old_p[i];
                }
            }

            //  Find magnitude
            float mag_p = 0;
            for (uint32_t i = 0; i < n; ++i)
            {
                mag_p += conjf(p[i]) * p[i];
            }
            const float mag_p2 = mag_p;
            mag_p = sqrtf(mag_p);

            //  Normalize basis vector
            for (uint32_t i = 0; i < n; ++i)
            {
                p[i] /= mag_p;
            }

            //  Apply previous Givens rotations to the new column of R
            for (uint32_t l = 0; l < k - 1; ++l)
            {
                const _Complex float tmp = ck[l] * h[l] + sk[l] * h[l + 1];
                h[l + 1] = -sk[l] * h[l] + ck[l] * h[l + 1];
                h[l] = tmp;
            }

            //  Compute the new givens rotation
            const float rho = sqrtf(mag_p2 + conjf(h[k - 1]) * h[k - 1]);
            const _Complex float c_new = h[k - 1] / rho;
            const _Complex float s_new = mag_p / rho;
            ck[k - 1] = c_new;
            sk[k - 1] = s_new;
            h[k - 1] = c_new * h[k - 1] + s_new * mag_p;

            jmtxc_matrix_brm_set_col(r, k - 1, h);

            g[k] = -s_new * g[k - 1];
            g[k - 1] = c_new * g[k - 1];

            r_mag = cabsf(g[k]);
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
        for (uint_fast32_t row = 0; row < k; ++row)
        {
            const uint_fast32_t i = k - 1 - row;
            _Complex float *elements;
            jmtxc_matrix_brm_get_row(r, i, &elements);
            _Complex float sum = 0;
            for (uint_fast32_t j = 1; j < row + 1; ++j)
            {
                sum += alpha[i + j] * elements[j];
            }
            alpha[i] = (g[i] - sum) / elements[0];
        }

        //  Compute improvement to x
        for (uint32_t i = 0; i < n; ++i)
        {
            for (uint32_t j = 0; j < k; ++j)
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
 * basis to solve a problem A x = y.
 *
 *
 * @param mtx system matrix A
 * @param n size of the problem
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
jmtx_result jmtxcs_solve_iterative_gmresm_cds(
    const jmtxc_matrix_cds *mtx, uint32_t n, const _Complex float y[JMTX_ARRAY_ATTRIB(static restrict n)],
    _Complex float x[JMTX_ARRAY_ATTRIB(static restrict n)], uint32_t m, jmtxc_matrix_brm *r,
    _Complex float aux_vec1[JMTX_ARRAY_ATTRIB(restrict m)], _Complex float aux_vec2[JMTX_ARRAY_ATTRIB(restrict m)],
    _Complex float aux_vec3[JMTX_ARRAY_ATTRIB(restrict m)], _Complex float aux_vec4[JMTX_ARRAY_ATTRIB(restrict m)],
    _Complex float aux_vec5[JMTX_ARRAY_ATTRIB(restrict m)], _Complex float aux_vecs[JMTX_ARRAY_ATTRIB(restrict m * n)],
    jmtxf_solver_arguments *args)
{
    if (mtx->base.type != JMTXC_TYPE_CDS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (mtx->base.rows != n || mtx->base.cols != n || mtx->main_diagonal)
    {
        return JMTX_RESULT_BAD_MATRIX;
    }
    if (r->base.type != JMTXC_TYPE_BRM)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (r->base.rows != m || r->base.cols != m)
    {
        return JMTX_RESULT_BAD_MATRIX;
    }
    if (r->upper_bandwidth != m - 1 || r->lower_bandwidth != 0)
    {
        return JMTX_RESULT_BAD_MATRIX;
    }
    if (m == 0)
    {
        return JMTX_RESULT_BAD_PARAM;
    }
    return jmtxc_solve_iterative_gmresm_cds(mtx, y, x, m, r, aux_vec1, aux_vec2, aux_vec3, aux_vec4, aux_vec5, aux_vecs,
                                            args);
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
jmtx_result jmtxc_solve_iterative_gmresm_rpc_jacobi_cds(
    const jmtxc_matrix_cds *mtx, const _Complex float *restrict y, _Complex float *restrict x, uint32_t m,
    jmtxc_matrix_brm *r, _Complex float aux_vec1[JMTX_ARRAY_ATTRIB(restrict m)],
    _Complex float aux_vec2[JMTX_ARRAY_ATTRIB(restrict m)], _Complex float aux_vec3[JMTX_ARRAY_ATTRIB(restrict m)],
    _Complex float aux_vec4[JMTX_ARRAY_ATTRIB(restrict m)], _Complex float aux_vec5[JMTX_ARRAY_ATTRIB(restrict m)],
    _Complex float *restrict aux_vec6, _Complex float *restrict aux_vec7, _Complex float *restrict aux_vecs,
    jmtxf_solver_arguments *args)
{
    float err = 0, y_mag = 0, r_mag = 0;
    uint32_t n_iteration = 0;
    _Complex float *const d_inv = aux_vec6;
    const uint32_t n = mtx->base.rows;
    for (uint32_t i = 0; i < n; ++i)
    {
        d_inv[i] = 1 / mtx->main_diagonal[i];
        y_mag += conjf(y[i]) * y[i];
    }
    y_mag = sqrtf(y_mag);

    const uint32_t round_count = args->in_max_iterations / m;
    for (uint32_t round = 0; round < round_count; ++round)
    {
        _Complex float *const q = aux_vecs;
        _Complex float *const ck = aux_vec1;
        _Complex float *const sk = aux_vec2;
        _Complex float *const g = aux_vec3;
        _Complex float *const alpha = aux_vec4;
        _Complex float *const h = aux_vec5;
        _Complex float *p = q;

        jmtxc_matrix_cds_vector_multiply(mtx, x, p);
        for (uint32_t i = 0; i < n; ++i)
        {
            const _Complex float res = y[i] - p[i];
            r_mag += conjf(res) * res;
            p[i] = res;
        }
        r_mag = sqrtf(r_mag);
        err = r_mag / y_mag;
        if (err < args->in_convergence_criterion)
        {
            args->out_last_iteration = n_iteration;
            args->out_last_error = err;
            return JMTX_RESULT_SUCCESS;
        }
        for (uint32_t i = 0; i < n; ++i)
        {
            p[i] /= r_mag;
        }
        g[0] = r_mag;
        uint32_t k;
        for (k = 1; k < m; ++k)
        {
            //  Generate new basis vector
            for (uint32_t i = 0; i < n; ++i)
            {
                aux_vec7[i] = p[i] * d_inv[i];
            }
            jmtxc_matrix_cds_vector_multiply(mtx, aux_vec7, p + n);
            p += n;
            //  Make the basis orthogonal to other basis
            for (uint32_t l = 0; l < k; ++l)
            {
                h[l] = 0;
                const _Complex float *old_p = q + n * l;
                for (uint32_t i = 0; i < n; ++i)
                {
                    h[l] += conjf(old_p[i]) * p[i];
                }
                for (uint32_t i = 0; i < n; ++i)
                {
                    p[i] -= h[l] * old_p[i];
                }
            }

            //  Find magnitude
            float mag_p = 0;
            for (uint32_t i = 0; i < n; ++i)
            {
                mag_p += conjf(p[i]) * p[i];
            }
            const float mag_p2 = mag_p;
            mag_p = sqrtf(mag_p);

            //  Normalize basis vector
            for (uint32_t i = 0; i < n; ++i)
            {
                p[i] /= mag_p;
            }

            //  Apply previous Givens rotations to the new column of R
            for (uint32_t l = 0; l < k - 1; ++l)
            {
                const _Complex float tmp = ck[l] * h[l] + sk[l] * h[l + 1];
                h[l + 1] = -sk[l] * h[l] + ck[l] * h[l + 1];
                h[l] = tmp;
            }

            //  Compute the new givens rotation
            const float rho = sqrtf(mag_p2 + conjf(h[k - 1]) * h[k - 1]);
            const _Complex float c_new = h[k - 1] / rho;
            const _Complex float s_new = mag_p / rho;
            ck[k - 1] = c_new;
            sk[k - 1] = s_new;
            h[k - 1] = c_new * h[k - 1] + s_new * mag_p;

            jmtxc_matrix_brm_set_col(r, k - 1, h);

            g[k] = -s_new * g[k - 1];
            g[k - 1] = c_new * g[k - 1];

            r_mag = cabsf(g[k]);
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
        for (uint_fast32_t row = 0; row < k; ++row)
        {
            const uint_fast32_t i = k - 1 - row;
            _Complex float *elements;
            jmtxc_matrix_brm_get_row(r, i, &elements);
            _Complex float sum = 0;
            for (uint_fast32_t j = 1; j < row + 1; ++j)
            {
                sum += alpha[i + j] * elements[j];
            }
            alpha[i] = (g[i] - sum) / elements[0];
        }

        //  Compute improvement to x
        //      Apply preconditioner on x first
        for (uint32_t i = 0; i < n; ++i)
        {
            x[i] *= d_inv[i];
        }
        for (uint32_t i = 0; i < n; ++i)
        {
            for (uint32_t j = 0; j < k; ++j)
            {
                x[i] += alpha[j] * q[j * n + i];
            }
        }
        //      Undo the preconditioner
        for (uint32_t i = 0; i < n; ++i)
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
jmtx_result jmtxc_solve_iterative_gmresm_lpc_jacobi_cds(
    const jmtxc_matrix_cds *mtx, const _Complex float *restrict y, _Complex float *restrict x, uint32_t m,
    jmtxc_matrix_brm *r, _Complex float aux_vec1[JMTX_ARRAY_ATTRIB(restrict m)],
    _Complex float aux_vec2[JMTX_ARRAY_ATTRIB(restrict m)], _Complex float aux_vec3[JMTX_ARRAY_ATTRIB(restrict m)],
    _Complex float aux_vec4[JMTX_ARRAY_ATTRIB(restrict m)], _Complex float aux_vec5[JMTX_ARRAY_ATTRIB(restrict m)],
    _Complex float *restrict aux_vec6, _Complex float *restrict aux_vecs, jmtxf_solver_arguments *args)
{
    float err = 0, y_mag = 0, r_mag = 0;
    uint32_t n_iteration = 0;
    _Complex float *const d_inv = aux_vec6;
    const uint32_t n = mtx->base.rows;
    for (uint32_t i = 0; i < n; ++i)
    {
        d_inv[i] = 1 / mtx->main_diagonal[i];
        const _Complex float scaled_y = y[i] * d_inv[i];
        y_mag += conjf(scaled_y) * scaled_y;
    }
    y_mag = sqrtf(y_mag);

    const uint32_t round_count = args->in_max_iterations / m;
    for (uint32_t round = 0; round < round_count; ++round)
    {
        _Complex float *const q = aux_vecs;
        _Complex float *const ck = aux_vec1;
        _Complex float *const sk = aux_vec2;
        _Complex float *const g = aux_vec3;
        _Complex float *const alpha = aux_vec4;
        _Complex float *const h = aux_vec5;
        _Complex float *p = q;

        jmtxc_matrix_cds_vector_multiply(mtx, x, p);
        for (uint32_t i = 0; i < n; ++i)
        {
            x[i] *= d_inv[i];
        }
        for (uint32_t i = 0; i < n; ++i)
        {
            const _Complex float res = y[i] * d_inv[i] - p[i];
            r_mag += conjf(res) * res;
            p[i] = res;
        }
        r_mag = sqrtf(r_mag);
        err = r_mag / y_mag;
        if (err < args->in_convergence_criterion)
        {
            args->out_last_iteration = n_iteration;
            args->out_last_error = err;
            return JMTX_RESULT_SUCCESS;
        }
        for (uint32_t i = 0; i < n; ++i)
        {
            p[i] /= r_mag;
        }
        g[0] = r_mag;
        uint32_t k;
        for (k = 1; k < m; ++k)
        {
            //  Generate new basis vector
            jmtxc_matrix_cds_vector_multiply(mtx, p, p + n);
            p += n;
            for (uint32_t i = 0; i < n; ++i)
            {
                p[i] *= d_inv[i];
            }
            //  Make the basis orthogonal to other basis
            for (uint32_t l = 0; l < k; ++l)
            {
                h[l] = 0;
                const _Complex float *old_p = q + n * l;
                for (uint32_t i = 0; i < n; ++i)
                {
                    h[l] += conjf(old_p[i]) * p[i];
                }
                for (uint32_t i = 0; i < n; ++i)
                {
                    p[i] -= h[l] * old_p[i];
                }
            }

            //  Find magnitude
            float mag_p = 0;
            for (uint32_t i = 0; i < n; ++i)
            {
                mag_p += conjf(p[i]) * p[i];
            }
            const float mag_p2 = mag_p;
            mag_p = sqrtf(mag_p);

            //  Normalize basis vector
            for (uint32_t i = 0; i < n; ++i)
            {
                p[i] /= mag_p;
            }

            //  Apply previous Givens rotations to the new column of R
            for (uint32_t l = 0; l < k - 1; ++l)
            {
                const _Complex float tmp = ck[l] * h[l] + sk[l] * h[l + 1];
                h[l + 1] = -sk[l] * h[l] + ck[l] * h[l + 1];
                h[l] = tmp;
            }

            //  Compute the new givens rotation
            const float rho = sqrtf(mag_p2 + conjf(h[k - 1]) * h[k - 1]);
            const _Complex float c_new = h[k - 1] / rho;
            const _Complex float s_new = mag_p / rho;
            ck[k - 1] = c_new;
            sk[k - 1] = s_new;
            h[k - 1] = c_new * h[k - 1] + s_new * mag_p;

            jmtxc_matrix_brm_set_col(r, k - 1, h);

            g[k] = -s_new * g[k - 1];
            g[k - 1] = c_new * g[k - 1];

            r_mag = cabsf(g[k]);
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
        for (uint_fast32_t row = 0; row < k; ++row)
        {
            const uint_fast32_t i = k - 1 - row;
            _Complex float *elements;
            jmtxc_matrix_brm_get_row(r, i, &elements);
            _Complex float sum = 0;
            for (uint_fast32_t j = 1; j < row + 1; ++j)
            {
                sum += alpha[i + j] * elements[j];
            }
            alpha[i] = (g[i] - sum) / elements[0];
        }

        //  Compute improvement to x
        for (uint32_t i = 0; i < n; ++i)
        {
            for (uint32_t j = 0; j < k; ++j)
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
