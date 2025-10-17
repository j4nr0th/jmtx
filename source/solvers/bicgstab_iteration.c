#include "bicgstab_iteration.h"
#include "../matrices/band_row_major.h"
#include "../matrices/sparse_diagonal_compressed.h"
#include "../matrices/sparse_row_compressed.h"
#include <math.h>
#include <omp.h>

#include <stdio.h>

#include "lu_solving.h"

/**
 *  Solves the linear problem A x = y for a general matrix A by using the relations used for Bi-CG, but does not
 *  explicitly solve the adjoint problem, instead computing values by computing results of polynomial relations for it.
 *  Stabilized method also computes these indirectly by using a polynomial with a lower condition number, giving better
 *  convergence behaviour.
 *
 *  This version of the function does not check if its inputs are valid and just assumes they are.
 *
 * @param mtx system matrix A
 * @param y solution to the system A x = y
 * @param x the initial guess of x, which receives the final solution
 * @param aux_vec1 auxiliary memory used by the algorithm which needs to be the same size as x and y
 * @param aux_vec2 auxiliary memory used by the algorithm which needs to be the same size as x and y
 * @param aux_vec3 auxiliary memory used by the algorithm which needs to be the same size as x and y
 * @param aux_vec4 auxiliary memory used by the algorithm which needs to be the same size as x and y
 * @param aux_vec5 auxiliary memory used by the algorithm which needs to be the same size as x and y
 * @param aux_vec6 auxiliary memory used by the algorithm which needs to be the same size as x and y
 * @param args::in_convergence_criterion tolerance to determine if the solution is close enough
 * @param args::in_max_iterations number of iterations to stop at
 * @param args::out_last_error receives the value of the error criterion at the final iteration
 * @param args::out_last_iteration receives the number of the final iteration
 * @param args::opt_error_evolution (optional) pointer to an array of length max_iterations, that receives the error
 * value of each iteration
 * @return JMTX_RESULT_SUCCESS if solution converged, JMTX_RESULT_NOT_CONVERGED if solution did not converge in the
 * given number of iterations
 */
jmtx_result JMTX_NAME_TYPED(solve_iterative_bicgstab_crs)(
    const JMTX_NAME_TYPED(matrix_crs) * mtx, const JMTX_SCALAR_T *restrict y, JMTX_SCALAR_T *restrict x,
    JMTX_SCALAR_T *restrict aux_vec1, JMTX_SCALAR_T *restrict aux_vec2, JMTX_SCALAR_T *restrict aux_vec3,
    JMTX_SCALAR_T *restrict aux_vec4, JMTX_SCALAR_T *restrict aux_vec5, JMTX_SCALAR_T *restrict aux_vec6,
    JMTX_NAME_TYPED(solver_arguments) * args)
{
    const JMTX_INDEX_T n = mtx->base.rows;

    JMTX_SCALAR_T rho = 1, alpha = 1, omega = 1;

    JMTX_SCALAR_T *const r = aux_vec1;
    JMTX_SCALAR_T *const rQ = aux_vec2;
    JMTX_SCALAR_T *const p = aux_vec3;
    JMTX_SCALAR_T *const Ap = aux_vec4;
    JMTX_SCALAR_T *const s = aux_vec5;
    JMTX_SCALAR_T *const As = aux_vec6;

    JMTX_REAL_T err = 0;

    JMTX_NAME_TYPED(matrix_crs_vector_multiply)(mtx, x, r);
    JMTX_REAL_T y_mag = 0;
    for (JMTX_INDEX_T i = 0; i < n; ++i)
    {
        r[i] = y[i] - r[i];
        rQ[i] = r[i];
        p[i] = r[i];
        //        Ap[i] = 0;
        y_mag += JMTX_DOT(y[i], y[i]);
        err += JMTX_DOT(r[i], r[i]);
    }
    y_mag = JMTX_REAL_ROOT(y_mag);
    err = JMTX_REAL_ROOT(err) / y_mag;
    if (err < args->in_convergence_criterion)
    {
        args->out_last_error = err;
        args->out_last_iteration = 0;
        return JMTX_RESULT_SUCCESS;
    }

    JMTX_INDEX_T iter_count = 0;
    for (;;)
    {
        JMTX_NAME_TYPED(matrix_crs_vector_multiply)(mtx, p, Ap);
        JMTX_SCALAR_T rQAp = 0;
        for (JMTX_INDEX_T i = 0; i < n; ++i)
        {
            rQAp += JMTX_DOT(rQ[i], Ap[i]);
        }
        alpha = rho / rQAp;
        for (JMTX_INDEX_T i = 0; i < n; ++i)
        {
            x[i] = x[i] + alpha * p[i];
        }
        JMTX_REAL_T sksk_dp = 0;
        for (JMTX_INDEX_T i = 0; i < n; ++i)
        {
            const JMTX_SCALAR_T si = r[i] - alpha * Ap[i];
            s[i] = si;
            sksk_dp += JMTX_DOT(si, si);
        }
        err = JMTX_REAL_ROOT(sksk_dp) / y_mag;
        if (err < args->in_convergence_criterion)
        {
            break;
        }
        JMTX_NAME_TYPED(matrix_crs_vector_multiply)(mtx, s, As);
        JMTX_REAL_T sAAs_dp = 0;
        JMTX_SCALAR_T sAs_dp = 0;
        for (JMTX_INDEX_T i = 0; i < n; ++i)
        {
            sAAs_dp += JMTX_DOT(As[i], As[i]);
            sAs_dp += JMTX_DOT(s[i], As[i]);
        }
        omega = sAs_dp / sAAs_dp;
        for (JMTX_INDEX_T i = 0; i < n; ++i)
        {
            x[i] = x[i] + omega * s[i];
        }
        JMTX_REAL_T rkrk_dp = 0;
        for (JMTX_INDEX_T i = 0; i < n; ++i)
        {
            const JMTX_SCALAR_T ri = s[i] - omega * As[i];
            r[i] = ri;
            rkrk_dp += JMTX_DOT(ri, ri);
        }
        err = JMTX_REAL_ROOT(rkrk_dp) / y_mag;
        if (iter_count == args->in_max_iterations)
        {
            break;
        }
        if (args->opt_error_evolution)
        {
            args->opt_error_evolution[iter_count] = err;
        }
        iter_count += 1;
        if (err < args->in_convergence_criterion)
        {
            break;
        }

        JMTX_SCALAR_T rQrk_dp = 0;
        for (JMTX_INDEX_T i = 0; i < n; ++i)
        {
            rQrk_dp += JMTX_DOT(rQ[i], r[i]);
        }
        const JMTX_SCALAR_T beta = rQrk_dp / rho * alpha / omega;
        rho = rQrk_dp;
        for (JMTX_INDEX_T i = 0; i < n; ++i)
        {
            p[i] = r[i] + beta * (p[i] - omega * Ap[i]);
        }
    }

    args->out_last_iteration = iter_count;
    args->out_last_error = err;
    if (!isfinite(err) || err > args->in_convergence_criterion)
    {
        return JMTX_RESULT_NOT_CONVERGED;
    }

    return JMTX_RESULT_SUCCESS;
}

/**
 *  Solves the linear problem A x = y for a general matrix A by using the relations used for Bi-CG, but does not
 *  explicitly solve the adjoint problem, instead computing values by computing results of polynomial relations for it.
 *  Stabilized method also computes these indirectly by using a polynomial with a lower condition number, giving better
 *  convergence behaviour.
 *
 *  This version of the function does not check if its inputs are valid and just assumes they are.
 *
 * @param mtx system matrix A
 * @param y solution to the system A x = y
 * @param x the initial guess of x, which receives the final solution
 * @param aux_vec1 auxiliary memory used by the algorithm which needs to be the same size as x and y
 * @param aux_vec2 auxiliary memory used by the algorithm which needs to be the same size as x and y
 * @param aux_vec3 auxiliary memory used by the algorithm which needs to be the same size as x and y
 * @param aux_vec4 auxiliary memory used by the algorithm which needs to be the same size as x and y
 * @param aux_vec5 auxiliary memory used by the algorithm which needs to be the same size as x and y
 * @param aux_vec6 auxiliary memory used by the algorithm which needs to be the same size as x and y
 * @param args::in_convergence_criterion tolerance to determine if the solution is close enough
 * @param args::in_max_iterations number of iterations to stop at
 * @param args::out_last_error receives the value of the error criterion at the final iteration
 * @param args::out_last_iteration receives the number of the final iteration
 * @param args::opt_error_evolution (optional) pointer to an array of length max_iterations, that receives the error
 * value of each iteration
 * @return JMTX_RESULT_SUCCESS if solution converged, JMTX_RESULT_NOT_CONVERGED if solution did not converge in the
 * given number of iterations
 */
jmtx_result JMTX_NAME_TYPED(solve_iterative_bicgstab_cds)(
    const JMTX_NAME_TYPED(matrix_cds) * mtx, const JMTX_SCALAR_T *restrict y, JMTX_SCALAR_T *restrict x,
    JMTX_SCALAR_T *restrict aux_vec1, JMTX_SCALAR_T *restrict aux_vec2, JMTX_SCALAR_T *restrict aux_vec3,
    JMTX_SCALAR_T *restrict aux_vec4, JMTX_SCALAR_T *restrict aux_vec5, JMTX_SCALAR_T *restrict aux_vec6,
    JMTX_NAME_TYPED(solver_arguments) * args)
{
    const JMTX_INDEX_T n = mtx->base.rows;

    JMTX_SCALAR_T rho = 1, alpha = 1, omega = 1;

    JMTX_SCALAR_T *const r = aux_vec1;
    JMTX_SCALAR_T *const rQ = aux_vec2;
    JMTX_SCALAR_T *const p = aux_vec3;
    JMTX_SCALAR_T *const Ap = aux_vec4;
    JMTX_SCALAR_T *const s = aux_vec5;
    JMTX_SCALAR_T *const As = aux_vec6;

    JMTX_REAL_T err = 0;

    JMTX_NAME_TYPED(matrix_cds_vector_multiply)(mtx, x, r);
    JMTX_REAL_T y_mag = 0;
    for (JMTX_INDEX_T i = 0; i < n; ++i)
    {
        r[i] = y[i] - r[i];
        rQ[i] = r[i];
        p[i] = r[i];
        //        Ap[i] = 0;
        y_mag += JMTX_DOT(y[i], y[i]);
        err += JMTX_DOT(r[i], r[i]);
    }
    y_mag = JMTX_REAL_ROOT(y_mag);
    err = JMTX_REAL_ROOT(err) / y_mag;
    if (err < args->in_convergence_criterion)
    {
        args->out_last_error = err;
        args->out_last_iteration = 0;
        return JMTX_RESULT_SUCCESS;
    }

    JMTX_INDEX_T iter_count = 0;
    for (;;)
    {
        JMTX_NAME_TYPED(matrix_cds_vector_multiply)(mtx, p, Ap);
        JMTX_SCALAR_T rQAp = 0;
        for (JMTX_INDEX_T i = 0; i < n; ++i)
        {
            rQAp += JMTX_DOT(rQ[i], Ap[i]);
        }
        alpha = rho / rQAp;
        for (JMTX_INDEX_T i = 0; i < n; ++i)
        {
            x[i] = x[i] + alpha * p[i];
        }
        JMTX_REAL_T sksk_dp = 0;
        for (JMTX_INDEX_T i = 0; i < n; ++i)
        {
            const JMTX_SCALAR_T si = r[i] - alpha * Ap[i];
            s[i] = si;
            sksk_dp += JMTX_DOT(si, si);
        }
        err = JMTX_REAL_ROOT(sksk_dp) / y_mag;
        if (err < args->in_convergence_criterion)
        {
            break;
        }
        JMTX_NAME_TYPED(matrix_cds_vector_multiply)(mtx, s, As);
        JMTX_SCALAR_T sAs_dp = 0;
        JMTX_REAL_T sAAs_dp = 0;
        for (JMTX_INDEX_T i = 0; i < n; ++i)
        {
            sAAs_dp += JMTX_DOT(As[i], As[i]);
            sAs_dp += JMTX_DOT(s[i], As[i]);
        }
        omega = sAs_dp / sAAs_dp;
        for (JMTX_INDEX_T i = 0; i < n; ++i)
        {
            x[i] = x[i] + omega * s[i];
        }
        JMTX_REAL_T rkrk_dp = 0;
        for (JMTX_INDEX_T i = 0; i < n; ++i)
        {
            const JMTX_SCALAR_T ri = s[i] - omega * As[i];
            r[i] = ri;
            rkrk_dp += JMTX_DOT(ri, ri);
        }
        err = JMTX_REAL_ROOT(rkrk_dp) / y_mag;
        if (iter_count == args->in_max_iterations)
        {
            break;
        }
        if (args->opt_error_evolution)
        {
            args->opt_error_evolution[iter_count] = err;
        }
        iter_count += 1;
        if (err < args->in_convergence_criterion)
        {
            break;
        }

        JMTX_SCALAR_T rQrk_dp = 0;
        for (JMTX_INDEX_T i = 0; i < n; ++i)
        {
            rQrk_dp += JMTX_DOT(rQ[i], r[i]);
        }
        const JMTX_SCALAR_T beta = rQrk_dp / rho * alpha / omega;
        rho = rQrk_dp;
        for (JMTX_INDEX_T i = 0; i < n; ++i)
        {
            p[i] = r[i] + beta * (p[i] - omega * Ap[i]);
        }
    }

    args->out_last_iteration = iter_count;
    args->out_last_error = err;
    if (!isfinite(err) || err > args->in_convergence_criterion)
    {
        return JMTX_RESULT_NOT_CONVERGED;
    }

    return JMTX_RESULT_SUCCESS;
}

/**
 *  Solves the linear problem A x = y for a general matrix A by using the relations used for Bi-CG, but does not
 *  explicitly solve the adjoint problem, instead computing values by computing results of polynomial relations for it.
 *  Stabilized method also computes these indirectly by using a polynomial with a lower condition number, giving better
 *  convergence behaviour.
 *
 *  This version uses incomplete LU decomposition (ILU) of the matrix, which then allows for better convergence
 *  properties. The decomposition must be given to the function.
 *
 *  This version of the function does not check if its inputs are valid and just assumes they are.
 *
 * @param mtx system matrix A
 * @param l lower triangular matrix
 * @param u upper triangular matrix
 * @param y solution to the system A x = y
 * @param x the initial guess of x, which receives the final solution
 * @param aux_vec1 auxiliary memory used by the algorithm which needs to be the same size as x and y
 * @param aux_vec2 auxiliary memory used by the algorithm which needs to be the same size as x and y
 * @param aux_vec3 auxiliary memory used by the algorithm which needs to be the same size as x and y
 * @param aux_vec4 auxiliary memory used by the algorithm which needs to be the same size as x and y
 * @param aux_vec5 auxiliary memory used by the algorithm which needs to be the same size as x and y
 * @param aux_vec6 auxiliary memory used by the algorithm which needs to be the same size as x and y
 * @param aux_vec7 auxiliary memory used by the algorithm which needs to be the same size as x and y
 * @param aux_vec8 auxiliary memory used by the algorithm which needs to be the same size as x and y
 * @param args::in_convergence_criterion tolerance to determine if the solution is close enough
 * @param args::in_max_iterations number of iterations to stop at
 * @param args::out_last_error receives the value of the error criterion at the final iteration
 * @param args::out_last_iteration receives the number of the final iteration
 * @param args::opt_error_evolution (optional) pointer to an array of length max_iterations, that receives the error
 * value of each iteration
 * @return JMTX_RESULT_SUCCESS if solution converged, JMTX_RESULT_NOT_CONVERGED if solution did not converge in the
 * given number of iterations
 */
jmtx_result JMTX_NAME_TYPED(solve_iterative_pilubicgstab_crs)(
    const JMTX_NAME_TYPED(matrix_crs) * mtx, const JMTX_NAME_TYPED(matrix_crs) * l,
    const JMTX_NAME_TYPED(matrix_crs) * u, const JMTX_SCALAR_T *restrict y, JMTX_SCALAR_T *restrict x,
    JMTX_SCALAR_T *restrict aux_vec1, JMTX_SCALAR_T *restrict aux_vec2, JMTX_SCALAR_T *restrict aux_vec3,
    JMTX_SCALAR_T *restrict aux_vec4, JMTX_SCALAR_T *restrict aux_vec5, JMTX_SCALAR_T *restrict aux_vec6,
    JMTX_SCALAR_T *restrict aux_vec7, JMTX_SCALAR_T *restrict aux_vec8, JMTX_NAME_TYPED(solver_arguments) * args)
{
    const JMTX_INDEX_T n = mtx->base.rows;

    JMTX_SCALAR_T rho = 1, alpha = 1, omega = 1;

    JMTX_SCALAR_T *const r = aux_vec1;
    JMTX_SCALAR_T *const rQ = aux_vec2;
    JMTX_SCALAR_T *const p = aux_vec3;
    JMTX_SCALAR_T *const Ap = aux_vec4;
    JMTX_SCALAR_T *const s = aux_vec5;
    JMTX_SCALAR_T *const As = aux_vec6;
    JMTX_SCALAR_T *const p_hat = aux_vec7;
    JMTX_SCALAR_T *const s_hat = aux_vec8;

    JMTX_REAL_T err = 0;

    JMTX_NAME_TYPED(matrix_crs_vector_multiply)(mtx, x, r);
    JMTX_REAL_T y_mag = 0;
    for (JMTX_INDEX_T i = 0; i < n; ++i)
    {
        r[i] = y[i] - r[i];
        rQ[i] = r[i];
        p[i] = r[i];
        //        Ap[i] = 0;
        y_mag += JMTX_DOT(y[i], y[i]);
        err += JMTX_DOT(r[i], r[i]);
    }
    y_mag = JMTX_REAL_ROOT(y_mag);
    err = JMTX_REAL_ROOT(err) / y_mag;
    if (err < args->in_convergence_criterion)
    {
        args->out_last_error = err;
        args->out_last_iteration = 0;
        return JMTX_RESULT_SUCCESS;
    }

    JMTX_INDEX_T iter_count = 0;
    for (;;)
    {
        JMTX_NAME_TYPED(solve_direct_lu_crs)(l, u, p, p_hat);
        JMTX_NAME_TYPED(matrix_crs_vector_multiply)(mtx, p_hat, Ap);
        JMTX_SCALAR_T rQAp = 0;
        for (JMTX_INDEX_T i = 0; i < n; ++i)
        {
            rQAp += JMTX_DOT(rQ[i], Ap[i]);
        }
        if (rQAp == 0)
        {
            break;
        }
        alpha = rho / rQAp;
        for (JMTX_INDEX_T i = 0; i < n; ++i)
        {
            x[i] = x[i] + alpha * p_hat[i];
        }
        JMTX_REAL_T sksk_dp = 0;
        for (JMTX_INDEX_T i = 0; i < n; ++i)
        {
            const JMTX_SCALAR_T si = r[i] - alpha * Ap[i];
            s[i] = si;
            sksk_dp += JMTX_DOT(si, si);
        }
        err = JMTX_REAL_ROOT(sksk_dp) / y_mag;
        if (err < args->in_convergence_criterion)
        {
            break;
        }
        JMTX_NAME_TYPED(solve_direct_lu_crs)(l, u, s, s_hat);
        JMTX_NAME_TYPED(matrix_crs_vector_multiply)(mtx, s_hat, As);
        JMTX_SCALAR_T sAs_dp = 0, sAAs_dp = 0;
        for (JMTX_INDEX_T i = 0; i < n; ++i)
        {
            sAAs_dp += JMTX_DOT(As[i], As[i]);
            sAs_dp += JMTX_DOT(s[i], As[i]);
        }
        if (sAs_dp == 0)
        {
            break;
        }
        omega = sAs_dp / sAAs_dp;
        for (JMTX_INDEX_T i = 0; i < n; ++i)
        {
            x[i] = x[i] + omega * s_hat[i];
        }
        JMTX_SCALAR_T rkrk_dp = 0;
        for (JMTX_INDEX_T i = 0; i < n; ++i)
        {
            const JMTX_SCALAR_T ri = s[i] - omega * As[i];
            r[i] = ri;
            rkrk_dp += JMTX_DOT(ri, ri);
        }
        err = JMTX_REAL_ROOT(rkrk_dp) / y_mag;
        if (iter_count == args->in_max_iterations)
        {
            break;
        }
        if (args->opt_error_evolution)
        {
            args->opt_error_evolution[iter_count] = err;
        }
        iter_count += 1;
        if (err < args->in_convergence_criterion)
        {
            break;
        }

        JMTX_SCALAR_T rQrk_dp = 0;
        for (JMTX_INDEX_T i = 0; i < n; ++i)
        {
            rQrk_dp += JMTX_DOT(rQ[i], r[i]);
        }
        const JMTX_SCALAR_T beta = rQrk_dp / rho * alpha / omega;
        rho = rQrk_dp;
        if (rQrk_dp == 0)
        {
            break;
        }
        for (JMTX_INDEX_T i = 0; i < n; ++i)
        {
            p[i] = r[i] + beta * (p[i] - omega * Ap[i]);
        }
    }

    args->out_last_iteration = iter_count;
    args->out_last_error = err;
    if (!isfinite(err) || err > args->in_convergence_criterion)
    {
        return JMTX_RESULT_NOT_CONVERGED;
    }

    return JMTX_RESULT_SUCCESS;
}

/**
 *  Solves the linear problem A x = y for a general matrix A by using the relations used for Bi-CG, but does not
 *  explicitly solve the adjoint problem, instead computing values by computing results of polynomial relations for it.
 *  Stabilized method also computes these indirectly by using a polynomial with a lower condition number, giving better
 *  convergence behaviour.
 *
 *  This version uses incomplete LU decomposition (ILU) of the matrix, which then allows for better convergence
 *  properties. The decomposition must be given to the function.
 *
 *  This version of the function does not check if its inputs are valid and just assumes they are.
 *
 *  This version uses OpenMP to solve the problem in parallel using multiple threads.
 *
 * @param mtx system matrix A
 * @param l lower triangular matrix
 * @param u upper triangular matrix
 * @param y solution to the system A x = y
 * @param x the initial guess of x, which receives the final solution
 * @param aux_vec1 auxiliary memory used by the algorithm which needs to be the same size as x and y
 * @param aux_vec2 auxiliary memory used by the algorithm which needs to be the same size as x and y
 * @param aux_vec3 auxiliary memory used by the algorithm which needs to be the same size as x and y
 * @param aux_vec4 auxiliary memory used by the algorithm which needs to be the same size as x and y
 * @param aux_vec5 auxiliary memory used by the algorithm which needs to be the same size as x and y
 * @param aux_vec6 auxiliary memory used by the algorithm which needs to be the same size as x and y
 * @param aux_vec7 auxiliary memory used by the algorithm which needs to be the same size as x and y
 * @param aux_vec8 auxiliary memory used by the algorithm which needs to be the same size as x and y
 * @param args::in_convergence_criterion tolerance to determine if the solution is close enough
 * @param args::in_max_iterations number of iterations to stop at
 * @param args::out_last_error receives the value of the error criterion at the final iteration
 * @param args::out_last_iteration receives the number of the final iteration
 * @param args::opt_error_evolution (optional) pointer to an array of length max_iterations, that receives the error
 * value of each iteration
 * @return JMTX_RESULT_SUCCESS if solution converged, JMTX_RESULT_NOT_CONVERGED if solution did not converge in the
 * given number of iterations
 */
jmtx_result JMTX_NAME_TYPED(solve_iterative_pilubicgstab_crs_parallel)(
    const JMTX_NAME_TYPED(matrix_crs) * mtx, const JMTX_NAME_TYPED(matrix_crs) * l,
    const JMTX_NAME_TYPED(matrix_crs) * u, const JMTX_SCALAR_T *restrict y, JMTX_SCALAR_T *restrict x,
    JMTX_SCALAR_T *restrict aux_vec1, JMTX_SCALAR_T *restrict aux_vec2, JMTX_SCALAR_T *restrict aux_vec3,
    JMTX_SCALAR_T *restrict aux_vec4, JMTX_SCALAR_T *restrict aux_vec5, JMTX_SCALAR_T *restrict aux_vec6,
    JMTX_SCALAR_T *restrict aux_vec7, JMTX_SCALAR_T *restrict aux_vec8, JMTX_NAME_TYPED(solver_arguments) * args)
{
    const JMTX_INDEX_T n = mtx->base.rows;

    JMTX_SCALAR_T *const r = aux_vec1;
    JMTX_SCALAR_T *const rQ = aux_vec2;
    JMTX_SCALAR_T *const p = aux_vec3;
    JMTX_SCALAR_T *const Ap = aux_vec4;
    JMTX_SCALAR_T *const s = aux_vec5;
    JMTX_SCALAR_T *const As = aux_vec6;
    JMTX_SCALAR_T *const p_hat = aux_vec7;
    JMTX_SCALAR_T *const s_hat = aux_vec8;

    JMTX_SCALAR_T rQAp = 0;
    JMTX_REAL_T sksk_dp = 0;
    JMTX_SCALAR_T sAs_dp = 0;
    JMTX_REAL_T sAAs_dp = 0;
    JMTX_REAL_T y_mag = 0;
    JMTX_SCALAR_T rkrk_dp = 0;
    JMTX_SCALAR_T rQrk_dp = 0;

    JMTX_INDEX_T iter_count = 0;

    JMTX_REAL_T err = 0;
#pragma omp parallel default(none) shared(x, mtx, y, u, l, args, rQrk_dp, err, y_mag, r, rQ, p, Ap, s, As, p_hat,      \
                                              s_hat, n, iter_count, rQAp, sksk_dp, sAs_dp, sAAs_dp, rkrk_dp)
    {
        JMTX_SCALAR_T rho = 1, alpha = 1, omega = 1, beta = 1;
        // JMTX_NAME_TYPED(matrix_crs_vector_multiply(mtx, x, r);
#pragma omp for reduction(+ : err, y_mag) schedule(static)
        for (JMTX_INDEX_T i = 0; i < n; ++i)
        {
            const JMTX_SCALAR_T yv = JMTX_NAME_TYPED(matrix_crs_vector_multiply_row)(mtx, x, i);
            r[i] = y[i] - yv;
            rQ[i] = r[i];
            p[i] = r[i];
            //        Ap[i] = 0;
            y_mag += JMTX_DOT(y[i], y[i]);
            err += JMTX_DOT(r[i], r[i]);
        }
#pragma omp single
        {
            y_mag = JMTX_REAL_ROOT(y_mag);
            err = JMTX_REAL_ROOT(err) / y_mag;
        }

        while (err > args->in_convergence_criterion)
        {
#pragma omp single
            {
                JMTX_NAME_TYPED(solve_direct_lu_crs)(l, u, p, p_hat);
                rQAp = 0;
                sksk_dp = 0;
            }

            // JMTX_NAME_TYPED(matrix_crs_vector_multiply(mtx, p_hat, Ap);

#pragma omp for reduction(+ : rQAp) schedule(static)
            for (JMTX_INDEX_T i = 0; i < n; ++i)
            {
                Ap[i] = JMTX_NAME_TYPED(matrix_crs_vector_multiply_row)(mtx, p_hat, i);
                rQAp += JMTX_DOT(rQ[i], Ap[i]);
            }

            if (rQAp == 0)
            {
                break;
            }
            // #pragma omp single

            alpha = rho / rQAp;

#pragma omp for reduction(+ : sksk_dp) schedule(static)
            for (JMTX_INDEX_T i = 0; i < n; ++i)
            {
                x[i] = x[i] + alpha * p_hat[i];
                s[i] = r[i] - alpha * Ap[i];
                sksk_dp += s[i] * s[i];
            }

            // #pragma omp single
            {
                err = JMTX_REAL_ROOT(sksk_dp) / y_mag;
            }
            if (err < args->in_convergence_criterion)
            {
                break;
            }
#pragma omp single
            {
                JMTX_NAME_TYPED(solve_direct_lu_crs)(l, u, s, s_hat);
                sAs_dp = 0, sAAs_dp = 0, rkrk_dp = 0, rQrk_dp = 0;
            }

            // JMTX_NAME_TYPED(matrix_crs_vector_multiply(mtx, s_hat, As);
#pragma omp for reduction(+ : sAs_dp, sAAs_dp) schedule(static)
            for (JMTX_INDEX_T i = 0; i < n; ++i)
            {
                As[i] = JMTX_NAME_TYPED(matrix_crs_vector_multiply_row)(mtx, s_hat, i);
                sAAs_dp += JMTX_DOT(As[i], As[i]);
                sAs_dp += JMTX_DOT(s[i], As[i]);
            }

            if (sAs_dp == 0)
            {
                break;
            }
            // #pragma omp single
            {
                omega = sAs_dp / sAAs_dp;
            }
#pragma omp for reduction(+ : rkrk_dp) schedule(static)
            for (JMTX_INDEX_T i = 0; i < n; ++i)
            {
                x[i] = x[i] + omega * s_hat[i];
                r[i] = s[i] - omega * As[i];
                rkrk_dp += JMTX_DOT(r[i], r[i]);
            }
            // #pragma omp single
            {
                err = JMTX_REAL_ROOT(rkrk_dp) / y_mag;
            }
            if (iter_count >= args->in_max_iterations)
            {
                break;
            }
#pragma omp barrier
#pragma omp master
            {
                if (args->opt_error_evolution)
                {
                    args->opt_error_evolution[iter_count] = err;
                }
                iter_count += 1;
            }
            if (err < args->in_convergence_criterion)
            {
                break;
            }

#pragma omp for reduction(+ : rQrk_dp) schedule(static)
            for (JMTX_INDEX_T i = 0; i < n; ++i)
            {
                rQrk_dp += JMTX_DOT(rQ[i], r[i]);
            }

            if (rQrk_dp == 0)
            {
                break;
            }
            // #pragma omp single
            {
                beta = rQrk_dp / rho * alpha / omega;
                rho = rQrk_dp;
            }
#pragma omp for schedule(static)
            for (JMTX_INDEX_T i = 0; i < n; ++i)
            {
                p[i] = r[i] + beta * (p[i] - omega * Ap[i]);
            }
        }
    }

    args->out_last_iteration = iter_count;
    args->out_last_error = err;
    if (!isfinite(err) || err > args->in_convergence_criterion)
    {
        return JMTX_RESULT_NOT_CONVERGED;
    }

    return JMTX_RESULT_SUCCESS;
}
