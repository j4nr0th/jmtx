//
// Created by jan on 1.1.2024.
//

#include <math.h>
#include <assert.h>
#include "../matrices/sparse_row_compressed_internal.h"
#include "../matrices/sparse_diagonal_compressed_internal.h"
#include "../matrices/band_row_major_internal.h"
#include "../../../include/jmtx/double/solvers/generalized_minimum_residual_iteration.h"

jmtx_result jmtxd_solve_iterative_gmresm_crs(const jmtxd_matrix_crs* mtx, const double* restrict y, double* restrict x,
                                             uint32_t m, jmtxd_matrix_brm* r, double aux_vec1[restrict m],
                                             double aux_vec2[restrict m], double aux_vec3[restrict m],
                                             double aux_vec4[restrict m], double aux_vec5[restrict m],
                                             double* restrict aux_vecs, jmtxd_solver_arguments* args)
{
    double err = 0, y_mag = 0, r_mag = 0;
    uint32_t n_iteration = 0;
    const uint32_t n = mtx->base.rows;
    for (uint32_t i = 0; i < n; ++i)
    {
        y_mag += y[i] * y[i];
    }
    y_mag = sqrt(y_mag);

    const uint32_t round_count = args->in_max_iterations / m;
    for (uint32_t round = 0; round < round_count; ++round)
    {
        double* const q = aux_vecs;
        double* const ck = aux_vec1;
        double* const sk = aux_vec2;
        double* const g = aux_vec3;
        double* const alpha = aux_vec4;
        double* const h = aux_vec5;
        double* p = q;
        for (uint32_t i = 0; i < n; ++i)
        {
            const double res = y[i] - jmtxd_matrix_crs_vector_multiply_row(mtx, x, i);
            r_mag += res * res;
            p[i] = res;
        }
        r_mag = sqrt(r_mag);
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
            jmtxd_matrix_crs_vector_multiply(mtx, p, p + n);
            p += n;
            //  Make the basis orthogonal to other basis
            for (uint32_t l = 0; l < k; ++l)
            {
                h[l] = 0;
                const double* old_p = q + n * l;
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
            double mag_p = 0;
            for (uint32_t i = 0; i < n; ++i)
            {
                mag_p += p[i] * p[i];
            }
            const double mag_p2 = mag_p;
            mag_p = sqrt(mag_p);

            //  Normalize basis vector
            for (uint32_t i = 0; i < n; ++i)
            {
                p[i] /= mag_p;
            }

            //  Apply previous Givens rotations to the new column of R
            for (uint32_t l = 0; l < k - 1; ++l)
            {
                const double tmp = ck[l] * h[l] + sk[l] * h[l + 1];
                h[l + 1] = -sk[l] * h[l] + ck[l] * h[l + 1];
                h[l] = tmp;
            }

            //  Compute the new givens rotation
            const double rho = sqrt(mag_p2 + h[k - 1] * h[k - 1]);
            const double c_new = h[k - 1] / rho;
            const double s_new = mag_p / rho;
            ck[k - 1] = c_new;
            sk[k - 1] = s_new;
            h[k - 1] = c_new * h[k - 1] + s_new * mag_p;

            jmtxd_matrix_brm_set_col(r, k - 1, h);

            g[k] = -s_new * g[k - 1];
            g[k - 1] = c_new * g[k - 1];

            r_mag = fabs(g[k]);
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
            double* elements;
            jmtxd_matrix_brm_get_row(r,  i, &elements);
            double sum = 0;
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
 * @param args::opt_error_evolution (optional) pointer to an array of length max_iterations, that receives the error value of each
 * iteration
 * @return JMTX_RESULT_SUCCESS if solution converged, JMTX_RESULT_NOT_CONVERGED if solution did not converge in the
 * given number of iterations, other error codes for other errors
 */
jmtx_result jmtxds_solve_iterative_gmresm_crs(const jmtxd_matrix_crs* mtx, uint32_t n, const double y[restrict static n],
                                             double x[restrict static n], uint32_t m, jmtxd_matrix_brm* r,
                                             double aux_vec1[restrict m], double aux_vec2[restrict m],
                                             double aux_vec3[restrict m], double aux_vec4[restrict m],
                                             double aux_vec5[restrict m], double aux_vecs[restrict m * n],
                                             jmtxd_solver_arguments* args)
{
    if (mtx->base.type != JMTXD_TYPE_CRS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (mtx->base.rows != n || mtx->base.cols != n)
    {
        return JMTX_RESULT_BAD_MATRIX;
    }
    if (r->base.type != JMTXD_TYPE_BRM)
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
    return jmtxd_solve_iterative_gmresm_crs(mtx, y, x, m ,r, aux_vec1, aux_vec2, aux_vec3, aux_vec4, aux_vec5, aux_vecs,
                                           args);
}

uint32_t jmtxd_gmresm_round_cds(const jmtxd_matrix_cds* mtx, const uint32_t n, const uint32_t m, const double y_mag,
                                const double tol, const double residual[const restrict static n],
                                double x[const restrict static n], jmtxd_matrix_brm* r, double ck[const restrict m],
                                double sk[const restrict m], double g[const restrict m], double alpha[const restrict m],
                                double h[const restrict m], double p_mat[const restrict m * n])
{
    double* p = p_mat;
    uint32_t n_iteration = 0;
    double err, r_mag = 0;
    for (uint32_t i = 0; i < n; ++i)
    {
        r_mag += residual[i] * residual[i];
    }
    r_mag = sqrt(r_mag);
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
        jmtxd_matrix_cds_vector_multiply(mtx, p, p + n);
        p += n;
        //  Make the basis orthogonal to other basis
        for (uint32_t l = 0; l < k; ++l)
        {
            h[l] = 0;
            const double* old_p = p_mat + n * l;
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
        double mag_p = 0;
        for (uint32_t i = 0; i < n; ++i)
        {
            mag_p += p[i] * p[i];
        }
        const double mag_p2 = mag_p;
        mag_p = sqrt(mag_p);

        //  Normalize basis vector
        for (uint32_t i = 0; i < n; ++i)
        {
            p[i] /= mag_p;
        }

        //  Apply previous Givens rotations to the new column of R
        for (uint32_t l = 0; l < k - 1; ++l)
        {
            const double tmp = ck[l] * h[l] + sk[l] * h[l + 1];
            h[l + 1] = -sk[l] * h[l] + ck[l] * h[l + 1];
            h[l] = tmp;
        }

        //  Compute the new givens rotation
        const double rho = sqrt(mag_p2 + h[k - 1] * h[k - 1]);
        const double c_new = h[k - 1] / rho;
        const double s_new = mag_p / rho;
        ck[k - 1] = c_new;
        sk[k - 1] = s_new;
        h[k - 1] = c_new * h[k - 1] + s_new * mag_p;

        jmtxd_matrix_brm_set_col(r, k - 1, h);

        g[k] = -s_new * g[k - 1];
        g[k - 1] = c_new * g[k - 1];

        r_mag = fabs(g[k]);
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
        double* elements;
        jmtxd_matrix_brm_get_row(r,  i, &elements);
        double sum = 0;
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

uint32_t jmtxd_gmresm_round_crs(const jmtxd_matrix_crs* mtx, const uint32_t n, const uint32_t m, const double y_mag,
                                const double tol, const double residual[const restrict static n],
                                double x[const restrict static n], jmtxd_matrix_brm* r, double ck[const restrict m],
                                double sk[const restrict m], double g[const restrict m], double alpha[const restrict m],
                                double h[const restrict m], double p_mat[const restrict m * n])
{
    double* p = p_mat;
    uint32_t n_iteration = 0;
    double err, r_mag = 0;
    for (uint32_t i = 0; i < n; ++i)
    {
        r_mag += residual[i] * residual[i];
    }
    r_mag = sqrt(r_mag);
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
        jmtxd_matrix_crs_vector_multiply(mtx, p, p + n);
        p += n;
        //  Make the basis orthogonal to other basis
        for (uint32_t l = 0; l < k; ++l)
        {
            h[l] = 0;
            const double* old_p = p_mat + n * l;
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
        double mag_p = 0;
        for (uint32_t i = 0; i < n; ++i)
        {
            mag_p += p[i] * p[i];
        }
        const double mag_p2 = mag_p;
        mag_p = sqrt(mag_p);

        //  Normalize basis vector
        for (uint32_t i = 0; i < n; ++i)
        {
            p[i] /= mag_p;
        }

        //  Apply previous Givens rotations to the new column of R
        for (uint32_t l = 0; l < k - 1; ++l)
        {
            const double tmp = ck[l] * h[l] + sk[l] * h[l + 1];
            h[l + 1] = -sk[l] * h[l] + ck[l] * h[l + 1];
            h[l] = tmp;
        }

        //  Compute the new givens rotation
        const double rho = sqrt(mag_p2 + h[k - 1] * h[k - 1]);
        const double c_new = h[k - 1] / rho;
        const double s_new = mag_p / rho;
        ck[k - 1] = c_new;
        sk[k - 1] = s_new;
        h[k - 1] = c_new * h[k - 1] + s_new * mag_p;

        jmtxd_matrix_brm_set_col(r, k - 1, h);

        g[k] = -s_new * g[k - 1];
        g[k - 1] = c_new * g[k - 1];

        r_mag = fabs(g[k]);
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
        double* elements;
        jmtxd_matrix_brm_get_row(r,  i, &elements);
        double sum = 0;
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
jmtx_result jmtxd_solve_iterative_gmresm_cds(const jmtxd_matrix_cds* mtx, const double* restrict y, double* restrict x,
                                             uint32_t m, jmtxd_matrix_brm* r, double aux_vec1[restrict m],
                                             double aux_vec2[restrict m], double aux_vec3[restrict m],
                                             double aux_vec4[restrict m], double aux_vec5[restrict m],
                                             double* restrict aux_vecs, jmtxd_solver_arguments* args){
    double err = 0, y_mag = 0, r_mag = 0;
    uint32_t n_iteration = 0;
    const uint32_t n = mtx->base.rows;
    for (uint32_t i = 0; i < n; ++i)
    {
        y_mag += y[i] * y[i];
    }
    y_mag = sqrt(y_mag);

    const uint32_t round_count = args->in_max_iterations / m;
    for (uint32_t round = 0; round < round_count; ++round)
    {
        double* const q = aux_vecs;
        double* const ck = aux_vec1;
        double* const sk = aux_vec2;
        double* const g = aux_vec3;
        double* const alpha = aux_vec4;
        double* const h = aux_vec5;
        double* p = q;
        jmtxd_matrix_cds_vector_multiply(mtx, x, p);
        for (uint32_t i = 0; i < n; ++i)
        {
            const double res = y[i] - p[i];
            r_mag += res * res;
            p[i] = res;
        }
        r_mag = sqrt(r_mag);
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
            jmtxd_matrix_cds_vector_multiply(mtx, p, p + n);
            p += n;
            //  Make the basis orthogonal to other basis
            for (uint32_t l = 0; l < k; ++l)
            {
                h[l] = 0;
                const double* old_p = q + n * l;
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
            double mag_p = 0;
            for (uint32_t i = 0; i < n; ++i)
            {
                mag_p += p[i] * p[i];
            }
            const double mag_p2 = mag_p;
            mag_p = sqrt(mag_p);

            //  Normalize basis vector
            for (uint32_t i = 0; i < n; ++i)
            {
                p[i] /= mag_p;
            }

            //  Apply previous Givens rotations to the new column of R
            for (uint32_t l = 0; l < k - 1; ++l)
            {
                const double tmp = ck[l] * h[l] + sk[l] * h[l + 1];
                h[l + 1] = -sk[l] * h[l] + ck[l] * h[l + 1];
                h[l] = tmp;
            }

            //  Compute the new givens rotation
            const double rho = sqrt(mag_p2 + h[k - 1] * h[k - 1]);
            const double c_new = h[k - 1] / rho;
            const double s_new = mag_p / rho;
            ck[k - 1] = c_new;
            sk[k - 1] = s_new;
            h[k - 1] = c_new * h[k - 1] + s_new * mag_p;

            jmtxd_matrix_brm_set_col(r, k - 1, h);

            g[k] = -s_new * g[k - 1];
            g[k - 1] = c_new * g[k - 1];

            r_mag = fabs(g[k]);
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
            double* elements;
            jmtxd_matrix_brm_get_row(r,  i, &elements);
            double sum = 0;
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
 * @param args::opt_error_evolution (optional) pointer to an array of length max_iterations, that receives the error value of each
 * iteration
 * @return JMTX_RESULT_SUCCESS if solution converged, JMTX_RESULT_NOT_CONVERGED if solution did not converge in the
 * given number of iterations, other error codes for other errors
 */
jmtx_result jmtxds_solve_iterative_gmresm_cds(const jmtxd_matrix_cds* mtx, uint32_t n, const double y[static restrict n],
                                             double x[static restrict n], uint32_t m, jmtxd_matrix_brm* r,
                                             double aux_vec1[restrict m], double aux_vec2[restrict m],
                                             double aux_vec3[restrict m], double aux_vec4[restrict m],
                                             double aux_vec5[restrict m], double aux_vecs[restrict n * m],
                                             jmtxd_solver_arguments* args)
{
    if (mtx->base.type != JMTXD_TYPE_CDS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (mtx->base.rows != n || mtx->base.cols != n || mtx->main_diagonal)
    {
        return JMTX_RESULT_BAD_MATRIX;
    }
    if (r->base.type != JMTXD_TYPE_BRM)
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
    return jmtxd_solve_iterative_gmresm_cds(mtx, y, x, m ,r, aux_vec1, aux_vec2, aux_vec3, aux_vec4, aux_vec5, aux_vecs,
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
 * @param args::opt_error_evolution (optional) pointer to an array of length max_iterations, that receives the error value of each
 * iteration
 * @return JMTX_RESULT_SUCCESS if solution converged, JMTX_RESULT_NOT_CONVERGED if solution did not converge in the
 * given number of iterations, other error codes for other errors
 */
jmtx_result jmtxd_solve_iterative_gmresm_rpc_jacobi_cds(const jmtxd_matrix_cds* mtx, const double* restrict y,
                                                        double* restrict x, uint32_t m, jmtxd_matrix_brm* r,
                                                        double aux_vec1[restrict m], double aux_vec2[restrict m],
                                                        double aux_vec3[restrict m], double aux_vec4[restrict m],
                                                        double aux_vec5[restrict m], double* restrict aux_vec6,
                                                        double* restrict aux_vec7, double* restrict aux_vecs,
                                                        jmtxd_solver_arguments* args)
{
    double err = 0, y_mag = 0, r_mag = 0;
    uint32_t n_iteration = 0;
    double* const d_inv = aux_vec6;
    const uint32_t n = mtx->base.rows;
    for (uint32_t i = 0; i < n; ++i)
    {
        d_inv[i] = 1 / mtx->main_diagonal[i];
        y_mag += y[i] * y[i];
    }
    y_mag = sqrt(y_mag);

    const uint32_t round_count = args->in_max_iterations / m;
    for (uint32_t round = 0; round < round_count; ++round)
    {
        double* const q = aux_vecs;
        double* const ck = aux_vec1;
        double* const sk = aux_vec2;
        double* const g = aux_vec3;
        double* const alpha = aux_vec4;
        double* const h = aux_vec5;
        double* p = q;

        jmtxd_matrix_cds_vector_multiply(mtx, x, p);
        for (uint32_t i = 0; i < n; ++i)
        {
            const double res = y[i] - p[i];
            r_mag += res * res;
            p[i] = res;
        }
        r_mag = sqrt(r_mag);
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
            jmtxd_matrix_cds_vector_multiply(mtx, aux_vec7, p + n);
            p += n;
            //  Make the basis orthogonal to other basis
            for (uint32_t l = 0; l < k; ++l)
            {
                h[l] = 0;
                const double* old_p = q + n * l;
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
            double mag_p = 0;
            for (uint32_t i = 0; i < n; ++i)
            {
                mag_p += p[i] * p[i];
            }
            const double mag_p2 = mag_p;
            mag_p = sqrt(mag_p);

            //  Normalize basis vector
            for (uint32_t i = 0; i < n; ++i)
            {
                p[i] /= mag_p;
            }

            //  Apply previous Givens rotations to the new column of R
            for (uint32_t l = 0; l < k - 1; ++l)
            {
                const double tmp = ck[l] * h[l] + sk[l] * h[l + 1];
                h[l + 1] = -sk[l] * h[l] + ck[l] * h[l + 1];
                h[l] = tmp;
            }

            //  Compute the new givens rotation
            const double rho = sqrt(mag_p2 + h[k - 1] * h[k - 1]);
            const double c_new = h[k - 1] / rho;
            const double s_new = mag_p / rho;
            ck[k - 1] = c_new;
            sk[k - 1] = s_new;
            h[k - 1] = c_new * h[k - 1] + s_new * mag_p;

            jmtxd_matrix_brm_set_col(r, k - 1, h);

            g[k] = -s_new * g[k - 1];
            g[k - 1] = c_new * g[k - 1];

            r_mag = fabs(g[k]);
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
            double* elements;
            jmtxd_matrix_brm_get_row(r,  i, &elements);
            double sum = 0;
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
jmtx_result jmtxd_solve_iterative_gmresm_lpc_jacobi_cds(const jmtxd_matrix_cds* mtx, const double* restrict y,
                                                        double* restrict x, uint32_t m, jmtxd_matrix_brm* r,
                                                        double aux_vec1[restrict m], double aux_vec2[restrict m],
                                                        double aux_vec3[restrict m], double aux_vec4[restrict m],
                                                        double aux_vec5[restrict m], double* restrict aux_vec6,
                                                        double* restrict aux_vecs, jmtxd_solver_arguments* args)
{
    double err = 0, y_mag = 0, r_mag = 0;
    uint32_t n_iteration = 0;
    double* const d_inv = aux_vec6;
    const uint32_t n = mtx->base.rows;
    for (uint32_t i = 0; i < n; ++i)
    {
        d_inv[i] = 1 / mtx->main_diagonal[i];
        const double scaled_y = y[i] * d_inv[i];
        y_mag += scaled_y * scaled_y;
    }
    y_mag = sqrt(y_mag);

    const uint32_t round_count = args->in_max_iterations / m;
    for (uint32_t round = 0; round < round_count; ++round)
    {
        double* const q = aux_vecs;
        double* const ck = aux_vec1;
        double* const sk = aux_vec2;
        double* const g = aux_vec3;
        double* const alpha = aux_vec4;
        double* const h = aux_vec5;
        double* p = q;

        jmtxd_matrix_cds_vector_multiply(mtx, x, p);
        for (uint32_t i = 0; i < n; ++i)
        {
            x[i] *= d_inv[i];
        }
        for (uint32_t i = 0; i < n; ++i)
        {
            const double res = y[i] * d_inv[i] - p[i];
            r_mag += res * res;
            p[i] = res;
        }
        r_mag = sqrt(r_mag);
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
            jmtxd_matrix_cds_vector_multiply(mtx, p, p + n);
            p += n;
            for (uint32_t i = 0; i < n; ++i)
            {
                p[i] *= d_inv[i];
            }
            //  Make the basis orthogonal to other basis
            for (uint32_t l = 0; l < k; ++l)
            {
                h[l] = 0;
                const double* old_p = q + n * l;
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
            double mag_p = 0;
            for (uint32_t i = 0; i < n; ++i)
            {
                mag_p += p[i] * p[i];
            }
            const double mag_p2 = mag_p;
            mag_p = sqrt(mag_p);

            //  Normalize basis vector
            for (uint32_t i = 0; i < n; ++i)
            {
                p[i] /= mag_p;
            }

            //  Apply previous Givens rotations to the new column of R
            for (uint32_t l = 0; l < k - 1; ++l)
            {
                const double tmp = ck[l] * h[l] + sk[l] * h[l + 1];
                h[l + 1] = -sk[l] * h[l] + ck[l] * h[l + 1];
                h[l] = tmp;
            }

            //  Compute the new givens rotation
            const double rho = sqrt(mag_p2 + h[k - 1] * h[k - 1]);
            const double c_new = h[k - 1] / rho;
            const double s_new = mag_p / rho;
            ck[k - 1] = c_new;
            sk[k - 1] = s_new;
            h[k - 1] = c_new * h[k - 1] + s_new * mag_p;

            jmtxd_matrix_brm_set_col(r, k - 1, h);

            g[k] = -s_new * g[k - 1];
            g[k - 1] = c_new * g[k - 1];

            r_mag = fabs(g[k]);
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
            double* elements;
            jmtxd_matrix_brm_get_row(r,  i, &elements);
            double sum = 0;
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

