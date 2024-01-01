//
// Created by jan on 1.1.2024.
//

#include <math.h>
#include <assert.h>
#include "../../../include/jmtx/float/solvers/generalized_minimum_residual_iteration.h"
#include "../matrices/sparse_row_compressed_internal.h"
#include "../matrices/sparse_diagonal_compressed_internal.h"

jmtx_result jmtx_solve_iterative_gmresm_crs(const jmtx_matrix_crs* mtx, const float* restrict y, float* restrict x,
                                             uint32_t m, jmtx_matrix_brm* r, float aux_vec1[restrict m],
                                             float aux_vec2[restrict m], float aux_vec3[restrict m],
                                             float aux_vec4[restrict m], float aux_vec5[restrict m],
                                             float* restrict aux_vecs, jmtx_solver_arguments* args)
{
    float err = 0, y_mag = 0, r_mag = 0;
    uint32_t n_iteration = 0;
    const uint32_t n = mtx->base.rows;
    for (uint32_t i = 0; i < n; ++i)
    {
        y_mag += y[i] * y[i];
    }
    y_mag = sqrtf(y_mag);

    const uint32_t round_count = args->in_max_iterations / m;
    for (uint32_t round = 0; round < round_count; ++round)
    {
        float* const q = aux_vecs;
        float* const ck = aux_vec1;
        float* const sk = aux_vec2;
        float* const g = aux_vec3;
        float* const alpha = aux_vec4;
        float* const h = aux_vec5;
        float* p = q;
        for (uint32_t i = 0; i < n; ++i)
        {
            const float res = y[i] - jmtx_matrix_crs_vector_multiply_row(mtx, x, i);
            r_mag += res * res;
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
            jmtx_matrix_crs_vector_multiply(mtx, p, p + n);
            p += n;
            //  Make the basis orthogonal to other basis
            for (uint32_t l = 0; l < k; ++l)
            {
                h[l] = 0;
                const float* old_p = q + n * l;
                for (uint32_t i = 0; i < n; ++i)
                {
                    h[l] += old_p[i] * p[i];
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
                mag_p += p[i] * p[i];
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
                const float tmp = ck[l] * h[l] + sk[l] * h[l + 1];
                h[l + 1] = -sk[l] * h[l] + ck[l] * h[l + 1];
                h[l] = tmp;
            }

            //  Compute the new givens rotation
            const float rho = sqrtf(mag_p2 + h[k - 1] * h[k - 1]);
            const float c_new = h[k - 1] / rho;
            const float s_new = mag_p / rho;
            ck[k - 1] = c_new;
            sk[k - 1] = s_new;
            h[k - 1] = c_new * h[k - 1] + s_new * mag_p;

            jmtx_matrix_brm_set_col(r, k - 1, h);

            g[k] = -s_new * g[k - 1];
            g[k - 1] = c_new * g[k - 1];

            r_mag = fabsf(g[k]);
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
            float* elements;
            jmtx_matrix_brm_get_row(r,  i, &elements);
            float sum = 0;
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

uint32_t jmtx_gmresm_round_cds(const jmtx_matrix_cds* mtx, const uint32_t n, const uint32_t m, const float y_mag,
                                const float tol, const float residual[const restrict static n],
                                float x[const restrict static n], jmtx_matrix_brm* r, float ck[const restrict m],
                                float sk[const restrict m], float g[const restrict m], float alpha[const restrict m],
                                float h[const restrict m], float p_mat[const restrict m * n])
{
    float* p = p_mat;
    uint32_t n_iteration = 0;
    float err, r_mag = 0;
    for (uint32_t i = 0; i < n; ++i)
    {
        r_mag += residual[i] * residual[i];
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
        jmtx_matrix_cds_vector_multiply(mtx, p, p + n);
        p += n;
        //  Make the basis orthogonal to other basis
        for (uint32_t l = 0; l < k; ++l)
        {
            h[l] = 0;
            const float* old_p = p_mat + n * l;
            for (uint32_t i = 0; i < n; ++i)
            {
                h[l] += old_p[i] * p[i];
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
            mag_p += p[i] * p[i];
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
            const float tmp = ck[l] * h[l] + sk[l] * h[l + 1];
            h[l + 1] = -sk[l] * h[l] + ck[l] * h[l + 1];
            h[l] = tmp;
        }

        //  Compute the new givens rotation
        const float rho = sqrtf(mag_p2 + h[k - 1] * h[k - 1]);
        const float c_new = h[k - 1] / rho;
        const float s_new = mag_p / rho;
        ck[k - 1] = c_new;
        sk[k - 1] = s_new;
        h[k - 1] = c_new * h[k - 1] + s_new * mag_p;

        jmtx_matrix_brm_set_col(r, k - 1, h);

        g[k] = -s_new * g[k - 1];
        g[k - 1] = c_new * g[k - 1];

        r_mag = fabsf(g[k]);
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
        float* elements;
        jmtx_matrix_brm_get_row(r,  i, &elements);
        float sum = 0;
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
 * @param args::opt_error_evolution (optional) pointer to an array of length max_iterations, that receives the error value of each
 * iteration
 * @return JMTX_RESULT_SUCCESS if solution converged, JMTX_RESULT_NOT_CONVERGED if solution did not converge in the
 * given number of iterations, other error codes for other errors
 */
jmtx_result jmtx_solve_iterative_gmresm_cds(const jmtx_matrix_cds* mtx, const float* restrict y, float* restrict x,
                                             uint32_t m, jmtx_matrix_brm* r, float aux_vec1[restrict m],
                                             float aux_vec2[restrict m], float aux_vec3[restrict m],
                                             float aux_vec4[restrict m], float aux_vec5[restrict m],
                                             float* restrict aux_vecs, jmtx_solver_arguments* args){
    float err = 0, y_mag = 0, r_mag = 0;
    uint32_t n_iteration = 0;
    const uint32_t n = mtx->base.rows;
    for (uint32_t i = 0; i < n; ++i)
    {
        y_mag += y[i] * y[i];
    }
    y_mag = sqrtf(y_mag);

    const uint32_t round_count = args->in_max_iterations / m;
    for (uint32_t round = 0; round < round_count; ++round)
    {
        float* const q = aux_vecs;
        float* const ck = aux_vec1;
        float* const sk = aux_vec2;
        float* const g = aux_vec3;
        float* const alpha = aux_vec4;
        float* const h = aux_vec5;
        float* p = q;
        jmtx_matrix_cds_vector_multiply(mtx, x, p);
        for (uint32_t i = 0; i < n; ++i)
        {
            const float res = y[i] - p[i];
            r_mag += res * res;
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
            jmtx_matrix_cds_vector_multiply(mtx, p, p + n);
            p += n;
            //  Make the basis orthogonal to other basis
            for (uint32_t l = 0; l < k; ++l)
            {
                h[l] = 0;
                const float* old_p = q + n * l;
                for (uint32_t i = 0; i < n; ++i)
                {
                    h[l] += old_p[i] * p[i];
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
                mag_p += p[i] * p[i];
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
                const float tmp = ck[l] * h[l] + sk[l] * h[l + 1];
                h[l + 1] = -sk[l] * h[l] + ck[l] * h[l + 1];
                h[l] = tmp;
            }

            //  Compute the new givens rotation
            const float rho = sqrtf(mag_p2 + h[k - 1] * h[k - 1]);
            const float c_new = h[k - 1] / rho;
            const float s_new = mag_p / rho;
            ck[k - 1] = c_new;
            sk[k - 1] = s_new;
            h[k - 1] = c_new * h[k - 1] + s_new * mag_p;

            jmtx_matrix_brm_set_col(r, k - 1, h);

            g[k] = -s_new * g[k - 1];
            g[k - 1] = c_new * g[k - 1];

            r_mag = fabsf(g[k]);
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
            float* elements;
            jmtx_matrix_brm_get_row(r,  i, &elements);
            float sum = 0;
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
 * @param args::opt_error_evolution (optional) pointer to an array of length max_iterations, that receives the error value of each
 * iteration
 * @return JMTX_RESULT_SUCCESS if solution converged, JMTX_RESULT_NOT_CONVERGED if solution did not converge in the
 * given number of iterations, other error codes for other errors
 */
jmtx_result jmtx_solve_iterative_gmresm_rpc_jacobi_cds(const jmtx_matrix_cds* mtx, const float* restrict y,
                                                        float* restrict x, uint32_t m, jmtx_matrix_brm* r,
                                                        float aux_vec1[restrict m], float aux_vec2[restrict m],
                                                        float aux_vec3[restrict m], float aux_vec4[restrict m],
                                                        float aux_vec5[restrict m], float* restrict aux_vec6,
                                                        float* restrict aux_vec7, float* restrict aux_vecs,
                                                        jmtx_solver_arguments* args)
{
    float err = 0, y_mag = 0, r_mag = 0;
    uint32_t n_iteration = 0;
    float* const d_inv = aux_vec6;
    const uint32_t n = mtx->base.rows;
    for (uint32_t i = 0; i < n; ++i)
    {
        d_inv[i] = 1 / mtx->main_diagonal[i];
        y_mag += y[i] * y[i];
    }
    y_mag = sqrtf(y_mag);

    const uint32_t round_count = args->in_max_iterations / m;
    for (uint32_t round = 0; round < round_count; ++round)
    {
        float* const q = aux_vecs;
        float* const ck = aux_vec1;
        float* const sk = aux_vec2;
        float* const g = aux_vec3;
        float* const alpha = aux_vec4;
        float* const h = aux_vec5;
        float* p = q;

        jmtx_matrix_cds_vector_multiply(mtx, x, p);
        for (uint32_t i = 0; i < n; ++i)
        {
            const float res = y[i] - p[i];
            r_mag += res * res;
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
            jmtx_matrix_cds_vector_multiply(mtx, aux_vec7, p + n);
            p += n;
            //  Make the basis orthogonal to other basis
            for (uint32_t l = 0; l < k; ++l)
            {
                h[l] = 0;
                const float* old_p = q + n * l;
                for (uint32_t i = 0; i < n; ++i)
                {
                    h[l] += old_p[i] * p[i];
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
                mag_p += p[i] * p[i];
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
                const float tmp = ck[l] * h[l] + sk[l] * h[l + 1];
                h[l + 1] = -sk[l] * h[l] + ck[l] * h[l + 1];
                h[l] = tmp;
            }

            //  Compute the new givens rotation
            const float rho = sqrtf(mag_p2 + h[k - 1] * h[k - 1]);
            const float c_new = h[k - 1] / rho;
            const float s_new = mag_p / rho;
            ck[k - 1] = c_new;
            sk[k - 1] = s_new;
            h[k - 1] = c_new * h[k - 1] + s_new * mag_p;

            jmtx_matrix_brm_set_col(r, k - 1, h);

            g[k] = -s_new * g[k - 1];
            g[k - 1] = c_new * g[k - 1];

            r_mag = fabsf(g[k]);
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
            float* elements;
            jmtx_matrix_brm_get_row(r,  i, &elements);
            float sum = 0;
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
 * @param args::opt_error_evolution (optional) pointer to an array of length max_iterations, that receives the error value of each
 * iteration
 * @return JMTX_RESULT_SUCCESS if solution converged, JMTX_RESULT_NOT_CONVERGED if solution did not converge in the
 * given number of iterations, other error codes for other errors
 */
jmtx_result jmtx_solve_iterative_gmresm_lpc_jacobi_cds(const jmtx_matrix_cds* mtx, const float* restrict y,
                                                        float* restrict x, uint32_t m, jmtx_matrix_brm* r,
                                                        float aux_vec1[restrict m], float aux_vec2[restrict m],
                                                        float aux_vec3[restrict m], float aux_vec4[restrict m],
                                                        float aux_vec5[restrict m], float* restrict aux_vec6,
                                                        float* restrict aux_vecs, jmtx_solver_arguments* args)
{
    float err = 0, y_mag = 0, r_mag = 0;
    uint32_t n_iteration = 0;
    float* const d_inv = aux_vec6;
    const uint32_t n = mtx->base.rows;
    for (uint32_t i = 0; i < n; ++i)
    {
        d_inv[i] = 1 / mtx->main_diagonal[i];
        const float scaled_y = y[i] * d_inv[i];
        y_mag += scaled_y * scaled_y;
    }
    y_mag = sqrtf(y_mag);

    const uint32_t round_count = args->in_max_iterations / m;
    for (uint32_t round = 0; round < round_count; ++round)
    {
        float* const q = aux_vecs;
        float* const ck = aux_vec1;
        float* const sk = aux_vec2;
        float* const g = aux_vec3;
        float* const alpha = aux_vec4;
        float* const h = aux_vec5;
        float* p = q;

        jmtx_matrix_cds_vector_multiply(mtx, x, p);
        for (uint32_t i = 0; i < n; ++i)
        {
            x[i] *= d_inv[i];
        }
        for (uint32_t i = 0; i < n; ++i)
        {
            const float res = y[i] * d_inv[i] - p[i];
            r_mag += res * res;
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
            jmtx_matrix_cds_vector_multiply(mtx, p, p + n);
            p += n;
            for (uint32_t i = 0; i < n; ++i)
            {
                p[i] *= d_inv[i];
            }
            //  Make the basis orthogonal to other basis
            for (uint32_t l = 0; l < k; ++l)
            {
                h[l] = 0;
                const float* old_p = q + n * l;
                for (uint32_t i = 0; i < n; ++i)
                {
                    h[l] += old_p[i] * p[i];
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
                mag_p += p[i] * p[i];
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
                const float tmp = ck[l] * h[l] + sk[l] * h[l + 1];
                h[l + 1] = -sk[l] * h[l] + ck[l] * h[l + 1];
                h[l] = tmp;
            }

            //  Compute the new givens rotation
            const float rho = sqrtf(mag_p2 + h[k - 1] * h[k - 1]);
            const float c_new = h[k - 1] / rho;
            const float s_new = mag_p / rho;
            ck[k - 1] = c_new;
            sk[k - 1] = s_new;
            h[k - 1] = c_new * h[k - 1] + s_new * mag_p;

            jmtx_matrix_brm_set_col(r, k - 1, h);

            g[k] = -s_new * g[k - 1];
            g[k - 1] = c_new * g[k - 1];

            r_mag = fabsf(g[k]);
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
            float* elements;
            jmtx_matrix_brm_get_row(r,  i, &elements);
            float sum = 0;
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

