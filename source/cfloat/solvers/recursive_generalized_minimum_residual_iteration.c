//
// Created by jan on 1.1.2024.
//
#include <complex.h>
#include <math.h>

#include "../../../include/jmtx/cfloat/solvers/recursive_generalized_minimum_residual_iteration.h"
#include "../matrices/band_row_major_internal.h"
#include "../matrices/sparse_diagonal_compressed_internal.h"

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
jmtx_result jmtxc_solve_iterative_gmresr_cds(
    const jmtxc_matrix_cds *mtx, const _Complex float *restrict y, _Complex float *restrict x, uint32_t m, uint32_t l,
    jmtxc_matrix_brm *r_mtx, _Complex float aux_vec1[JMTX_ARRAY_ATTRIB(restrict m)],
    _Complex float aux_vec2[JMTX_ARRAY_ATTRIB(restrict m)], _Complex float aux_vec3[JMTX_ARRAY_ATTRIB(restrict m)],
    _Complex float aux_vec4[JMTX_ARRAY_ATTRIB(restrict m)], _Complex float aux_vec5[JMTX_ARRAY_ATTRIB(restrict m)],
    _Complex float *restrict aux_vec6, _Complex float *restrict aux_vecs1, _Complex float *restrict aux_vecs2,
    _Complex float *restrict aux_vecs3, jmtxf_solver_arguments *args)
{
    float err, r_mag, y_mag;
    const uint32_t n = mtx->base.rows;
    _Complex float *const r = aux_vec6;
    jmtxc_matrix_cds_vector_multiply(mtx, x, r);
    y_mag = 0;
    r_mag = 0;
    for (uint32_t i = 0; i < n; ++i)
    {
        const _Complex float res = y[i] - r[i];
        y_mag += conjf(y[i]) * y[i];
        r_mag += conjf(res) * res;
        r[i] = res;
    }
    y_mag = sqrtf(y_mag);
    r_mag = sqrtf(r_mag);
    err = r_mag / y_mag;
    if (err < args->in_convergence_criterion || args->in_max_iterations == 0)
    {
        args->out_last_iteration = 0;
        args->out_last_error = err;
        return JMTX_RESULT_SUCCESS;
    }

    uint32_t n_iter = 0, vec_idx = 0, vec_cnt = 0;
    _Complex float *const s_vectors = aux_vecs2;
    _Complex float *const p_vectors = aux_vecs3;
    //  Begin iterations
    for (;;)
    {
        _Complex float *const s = s_vectors + n * vec_idx;
        _Complex float *const p = p_vectors + n * vec_idx;
        for (uint32_t i = 0; i < n; ++i)
        {
            p[i] = 0;
        }
        //  Use GMRES(m) to generate a correction to the current error
        (void)jmtxc_gmresm_round_cds(mtx, n, m, y_mag, args->in_convergence_criterion, r, p, r_mtx, aux_vec1, aux_vec2,
                                     aux_vec3, aux_vec4, aux_vec5, aux_vecs1);
        //  Generate conjugate vector
        jmtxc_matrix_cds_vector_multiply(mtx, p, s);

        //  Make conjugate vector orthonormal to previous l search vectors
        for (uint32_t i = 0; i < vec_cnt; ++i)
        {
            if (i == vec_idx)
            {
                continue;
            }
            _Complex float dp = 0;
            const _Complex float *other_s = s_vectors + i * n;
            const _Complex float *other_p = p_vectors + i * n;
            for (uint32_t j = 0; j < n; ++j)
            {
                dp += conjf(other_s[j]) * s[j];
            }
            for (uint32_t j = 0; j < n; ++j)
            {
                s[j] -= dp * other_s[j];
                p[j] -= dp * other_p[j];
            }
        }
        float s_mag = 0;
        for (uint32_t j = 0; j < n; ++j)
        {
            s_mag += conjf(s[j]) * s[j];
        }
        s_mag = sqrtf(s_mag);
        for (uint32_t j = 0; j < n; ++j)
        {
            s[j] /= s_mag;
            p[j] /= s_mag;
        }

        _Complex float alpha = 0;
        for (uint32_t j = 0; j < n; ++j)
        {
            alpha += conjf(r[j]) * s[j];
        }

        r_mag = 0;
        for (uint32_t j = 0; j < n; ++j)
        {
            r[j] -= alpha * s[j];
            r_mag += conjf(r[j]) * r[j];
            x[j] += alpha * p[j];
        }
        r_mag = sqrtf(r_mag);
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
 * @param n the size of the system
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
 * @param aux_vecs2 auxiliary memory for l vectors of the same size as x and y (n by l)
 * @param aux_vecs3 auxiliary memory for l vectors of the same size as x and y (n by l)
 * @param args::in_convergence_criterion tolerance to determine if the solution is close enough
 * @param args::in_max_iterations number of iterations to stop at
 * @param args::out_last_error receives the value of the error criterion at the final iteration
 * @param args::out_last_iteration receives the number of the final iteration
 * @param args::opt_error_evolution (optional) pointer to an array of length max_iterations, that receives the error
 * value of each iteration
 * @return JMTX_RESULT_SUCCESS if solution converged, JMTX_RESULT_NOT_CONVERGED if solution did not converge in the
 * given number of iterations, other error codes for other errors
 */
jmtx_result jmtxcs_solve_iterative_gmresr_cds(
    const jmtxc_matrix_cds *mtx, uint32_t n, const _Complex float y[JMTX_ARRAY_ATTRIB(static restrict n)],
    _Complex float x[JMTX_ARRAY_ATTRIB(static restrict n)], uint32_t m, uint32_t l, jmtxc_matrix_brm *r_mtx,
    _Complex float aux_vec1[JMTX_ARRAY_ATTRIB(restrict m)], _Complex float aux_vec2[JMTX_ARRAY_ATTRIB(restrict m)],
    _Complex float aux_vec3[JMTX_ARRAY_ATTRIB(restrict m)], _Complex float aux_vec4[JMTX_ARRAY_ATTRIB(restrict m)],
    _Complex float aux_vec5[JMTX_ARRAY_ATTRIB(restrict m)], _Complex float aux_vec6[JMTX_ARRAY_ATTRIB(restrict n)],
    _Complex float aux_vecs1[JMTX_ARRAY_ATTRIB(restrict m * n)],
    _Complex float aux_vecs2[JMTX_ARRAY_ATTRIB(restrict l * n)],
    _Complex float aux_vecs3[JMTX_ARRAY_ATTRIB(restrict l * n)], jmtxf_solver_arguments *args)
{
    if (mtx->base.type != JMTXC_TYPE_CDS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (mtx->base.rows != n || mtx->base.cols != n || mtx->main_diagonal)
    {
        return JMTX_RESULT_BAD_MATRIX;
    }
    if (r_mtx->base.type != JMTXC_TYPE_BRM)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (r_mtx->base.rows != m || r_mtx->base.cols != m)
    {
        return JMTX_RESULT_BAD_MATRIX;
    }
    if (r_mtx->upper_bandwidth != m - 1 || r_mtx->lower_bandwidth != 0)
    {
        return JMTX_RESULT_BAD_MATRIX;
    }
    if (m == 0)
    {
        return JMTX_RESULT_BAD_PARAM;
    }
    return jmtxc_solve_iterative_gmresr_cds(mtx, y, x, m, l, r_mtx, aux_vec1, aux_vec2, aux_vec3, aux_vec4, aux_vec5,
                                            aux_vec6, aux_vecs1, aux_vecs2, aux_vecs3, args);
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
jmtx_result jmtxc_solve_iterative_gmresr_crs(
    const jmtxc_matrix_crs *mtx, const _Complex float *restrict y, _Complex float *restrict x, uint32_t m, uint32_t l,
    jmtxc_matrix_brm *r_mtx, _Complex float aux_vec1[JMTX_ARRAY_ATTRIB(restrict m)],
    _Complex float aux_vec2[JMTX_ARRAY_ATTRIB(restrict m)], _Complex float aux_vec3[JMTX_ARRAY_ATTRIB(restrict m)],
    _Complex float aux_vec4[JMTX_ARRAY_ATTRIB(restrict m)], _Complex float aux_vec5[JMTX_ARRAY_ATTRIB(restrict m)],
    _Complex float *restrict aux_vec6, _Complex float *restrict aux_vecs1, _Complex float *restrict aux_vecs2,
    _Complex float *restrict aux_vecs3, jmtxf_solver_arguments *args)
{
    float err, r_mag, y_mag;
    const uint32_t n = mtx->base.rows;
    _Complex float *const r = aux_vec6;
    jmtxc_matrix_crs_vector_multiply(mtx, x, r);
    y_mag = 0;
    r_mag = 0;
    for (uint32_t i = 0; i < n; ++i)
    {
        const _Complex float res = y[i] - r[i];
        y_mag += conjf(y[i]) * y[i];
        r_mag += conjf(res) * res;
        r[i] = res;
    }
    y_mag = sqrtf(y_mag);
    r_mag = sqrtf(r_mag);
    err = r_mag / y_mag;
    if (err < args->in_convergence_criterion || args->in_max_iterations == 0)
    {
        args->out_last_iteration = 0;
        args->out_last_error = err;
        return JMTX_RESULT_SUCCESS;
    }

    uint32_t n_iter = 0, vec_idx = 0, vec_cnt = 0;
    _Complex float *const s_vectors = aux_vecs2;
    _Complex float *const p_vectors = aux_vecs3;
    //  Begin iterations
    for (;;)
    {
        _Complex float *const s = s_vectors + n * vec_idx;
        _Complex float *const p = p_vectors + n * vec_idx;
        for (uint32_t i = 0; i < n; ++i)
        {
            p[i] = 0;
        }
        //  Use GMRES(m) to generate a correction to the current error
        (void)jmtxc_gmresm_round_crs(mtx, n, m, y_mag, args->in_convergence_criterion, r, p, r_mtx, aux_vec1, aux_vec2,
                                     aux_vec3, aux_vec4, aux_vec5, aux_vecs1);
        //  Generate conjugate vector
        jmtxc_matrix_crs_vector_multiply(mtx, p, s);

        //  Make conjugate vector orthonormal to previous l search vectors
        for (uint32_t i = 0; i < vec_cnt; ++i)
        {
            if (i == vec_idx)
            {
                continue;
            }
            _Complex float dp = 0;
            const _Complex float *other_s = s_vectors + i * n;
            const _Complex float *other_p = p_vectors + i * n;
            for (uint32_t j = 0; j < n; ++j)
            {
                dp += conjf(other_s[j]) * s[j];
            }
            for (uint32_t j = 0; j < n; ++j)
            {
                s[j] -= dp * other_s[j];
                p[j] -= dp * other_p[j];
            }
        }
        float s_mag = 0;
        for (uint32_t j = 0; j < n; ++j)
        {
            s_mag += conjf(s[j]) * s[j];
        }
        s_mag = sqrtf(s_mag);
        for (uint32_t j = 0; j < n; ++j)
        {
            s[j] /= s_mag;
            p[j] /= s_mag;
        }

        _Complex float alpha = 0;
        for (uint32_t j = 0; j < n; ++j)
        {
            alpha += conjf(r[j]) * s[j];
        }

        r_mag = 0;
        for (uint32_t j = 0; j < n; ++j)
        {
            r[j] -= alpha * s[j];
            r_mag += conjf(r[j]) * r[j];
            x[j] += alpha * p[j];
        }
        r_mag = sqrtf(r_mag);
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
 * @param n the size of the system
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
 * @param aux_vecs2 auxiliary memory for l vectors of the same size as x and y (n by l)
 * @param aux_vecs3 auxiliary memory for l vectors of the same size as x and y (n by l)
 * @param args::in_convergence_criterion tolerance to determine if the solution is close enough
 * @param args::in_max_iterations number of iterations to stop at
 * @param args::out_last_error receives the value of the error criterion at the final iteration
 * @param args::out_last_iteration receives the number of the final iteration
 * @param args::opt_error_evolution (optional) pointer to an array of length max_iterations, that receives the error
 * value of each iteration
 * @return JMTX_RESULT_SUCCESS if solution converged, JMTX_RESULT_NOT_CONVERGED if solution did not converge in the
 * given number of iterations, other error codes for other errors
 */
jmtx_result jmtxcs_solve_iterative_gmresr_crs(
    const jmtxc_matrix_crs *mtx, uint32_t n, const _Complex float y[JMTX_ARRAY_ATTRIB(static restrict n)],
    _Complex float x[JMTX_ARRAY_ATTRIB(static restrict n)], uint32_t m, uint32_t l, jmtxc_matrix_brm *r_mtx,
    _Complex float aux_vec1[JMTX_ARRAY_ATTRIB(restrict m)], _Complex float aux_vec2[JMTX_ARRAY_ATTRIB(restrict m)],
    _Complex float aux_vec3[JMTX_ARRAY_ATTRIB(restrict m)], _Complex float aux_vec4[JMTX_ARRAY_ATTRIB(restrict m)],
    _Complex float aux_vec5[JMTX_ARRAY_ATTRIB(restrict m)], _Complex float aux_vec6[JMTX_ARRAY_ATTRIB(restrict n)],
    _Complex float aux_vecs1[JMTX_ARRAY_ATTRIB(restrict m * n)],
    _Complex float aux_vecs2[JMTX_ARRAY_ATTRIB(restrict l * n)],
    _Complex float aux_vecs3[JMTX_ARRAY_ATTRIB(restrict l * n)], jmtxf_solver_arguments *args)
{
    if (mtx->base.type != JMTXC_TYPE_CRS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (mtx->base.rows != n || mtx->base.cols != n)
    {
        return JMTX_RESULT_BAD_MATRIX;
    }
    if (r_mtx->base.type != JMTXC_TYPE_BRM)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (r_mtx->base.rows != m || r_mtx->base.cols != m)
    {
        return JMTX_RESULT_BAD_MATRIX;
    }
    if (r_mtx->upper_bandwidth != m - 1 || r_mtx->lower_bandwidth != 0)
    {
        return JMTX_RESULT_BAD_MATRIX;
    }
    if (m == 0)
    {
        return JMTX_RESULT_BAD_PARAM;
    }
    return jmtxc_solve_iterative_gmresr_crs(mtx, y, x, m, l, r_mtx, aux_vec1, aux_vec2, aux_vec3, aux_vec4, aux_vec5,
                                            aux_vec6, aux_vecs1, aux_vecs2, aux_vecs3, args);
}
