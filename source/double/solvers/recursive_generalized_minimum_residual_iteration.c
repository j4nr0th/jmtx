//
// Created by jan on 1.1.2024.
//
#include <math.h>

#include "../../../include/jmtx/double/solvers/recursive_generalized_minimum_residual_iteration.h"
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
 * @param args::opt_error_evolution (optional) pointer to an array of length max_iterations, that receives the error value of each
 * iteration
 * @return JMTX_RESULT_SUCCESS if solution converged, JMTX_RESULT_NOT_CONVERGED if solution did not converge in the
 * given number of iterations, other error codes for other errors
 */
jmtx_result jmtxd_solve_iterative_gmresr_cds(const jmtxd_matrix_cds* mtx, const double* restrict y, double* restrict x,
                                             uint32_t m, uint32_t l, jmtxd_matrix_brm* r_mtx,
                                             double aux_vec1[restrict m], double aux_vec2[restrict m],
                                             double aux_vec3[restrict m], double aux_vec4[restrict m],
                                             double aux_vec5[restrict m], double* restrict aux_vec6,
                                             double* restrict aux_vecs1, double* restrict aux_vecs2,
                                             double* restrict aux_vecs3, jmtxd_solver_arguments* args)
{
    double err, r_mag, y_mag;
    const uint32_t n = mtx->base.rows;
    double* const r = aux_vec6;
    jmtxd_matrix_cds_vector_multiply(mtx, x, r);
    y_mag = 0;
    r_mag = 0;
    for (uint32_t i = 0; i < n; ++i)
    {
        const double res = y[i] - r[i];
        y_mag += y[i] * y[i];
        r_mag += res * res;
        r[i] = res;
    }
    y_mag = sqrt(y_mag);
    r_mag = sqrt(r_mag);
    err = r_mag / y_mag;
    if (err < args->in_convergence_criterion || args->in_max_iterations == 0)
    {
        args->out_last_iteration = 0;
        args->out_last_error = err;
        return JMTX_RESULT_SUCCESS;
    }

    uint32_t n_iter = 0, vec_idx = 0, vec_cnt = 0;
    double* const s_vectors = aux_vecs2;
    double* const p_vectors = aux_vecs3;
    //  Begin iterations
    for (;;)
    {
        double* const s = s_vectors + n * vec_idx;
        double* const p = p_vectors + n * vec_idx;
        for (uint32_t i = 0; i < n; ++i)
        {
            p[i] = 0;
        }
        //  Use GMRES(m) to generate a correction to the current error
        (void)jmtxd_gmresm_round_cds(mtx, n, m, y_mag, args->in_convergence_criterion, r,
                               p, r_mtx, aux_vec1, aux_vec2, aux_vec3,
                               aux_vec4, aux_vec5, aux_vecs1);
        //  Generate conjugate vector
        jmtxd_matrix_cds_vector_multiply(mtx, p, s);

        //  Make conjugate vector orthonormal to previous l search vectors
        for (uint32_t i = 0; i < vec_cnt; ++i)
        {
            if (i == vec_idx)
            {
                continue;
            }
            double dp = 0;
            const double* other_s = s_vectors + i * n;
            const double* other_p = p_vectors + i * n;
            for (uint32_t j = 0; j < n; ++j)
            {
                dp += other_s[j] * s[j];
            }
            for (uint32_t j = 0; j < n; ++j)
            {
                s[j] -= dp * other_s[j];
                p[j] -= dp * other_p[j];
            }
        }
        double s_mag = 0;
        for (uint32_t j = 0; j < n; ++j)
        {
            s_mag += s[j] * s[j];
        }
        s_mag = sqrt(s_mag);
        for (uint32_t j = 0; j < n; ++j)
        {
            s[j] /= s_mag;
            p[j] /= s_mag;
        }

        double alpha = 0;
        for (uint32_t j = 0; j < n; ++j)
        {
            alpha += r[j] * s[j];
        }

        r_mag = 0;
        for (uint32_t j = 0; j < n; ++j)
        {
            r[j] -= alpha * s[j];
            r_mag += r[j] * r[j];
            x[j] += alpha * p[j];
        }
        r_mag = sqrt(r_mag);
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
