// Automatically generated from source/double/solvers/bicgstab_iteration.c on Sun Dec 17 15:54:41 2023
//
// Created by jan on 17.6.2022.
//

#include <math.h>
#include "../matrices/sparse_row_compressed_internal.h"
#include "../matrices/sparse_diagonal_compressed_internal.h"
#include "../matrices/band_row_major_internal.h"
#include "../../../include/jmtx/cdouble/solvers/bicgstab_iteration.h"
#include "../../../include/jmtx/cdouble/solvers/lu_solving.h"
#include <complex.h>


/**
 *  Solves the linear problem A x = y for a general matrix A by using the relations used for Bi-CG, but does not
 *  explicitly solve the adjoint problem, instead computing values by computing results of polynomial relations for it.
 *  Stabilized method also computes these indirectly by using a polynomial with a lower condition number, giving better
 *  convergence behaviour.
 *
 *  This version of the funciton does not check if its inputs are valid and just assumes they are.
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
jmtx_result jmtxz_solve_iterative_bicgstab_crs(
        const jmtxz_matrix_crs* mtx, const _Complex double* restrict y, _Complex double* restrict x, _Complex double* restrict aux_vec1,
        _Complex double* restrict aux_vec2, _Complex double* restrict aux_vec3, _Complex double* restrict aux_vec4, _Complex double* restrict aux_vec5,
        _Complex double* restrict aux_vec6, jmtxd_solver_arguments* args)
{
    const uint32_t n = mtx->base.rows;

    _Complex double rho = 1, alpha = 1, omega = 1;

    _Complex double* const r = aux_vec1;
    _Complex double* const rQt = aux_vec2;
    _Complex double* const p = aux_vec3;
    _Complex double* const Ap = aux_vec4;
    _Complex double* const s = aux_vec5;
    _Complex double* const As = aux_vec6;

    double err = 0;

    jmtxz_matrix_crs_vector_multiply(mtx, x, r);
    double y_mag = 0;
    for (uint32_t i = 0; i < n; ++i)
    {
        r[i] = y[i] - r[i];
        rQt[i] = conj(r[i]);
        p[i] = r[i];
//        Ap[i] = 0;
        y_mag += conj(y[i]) * y[i];
        err += rQt[i] * r[i];
    }
    y_mag = sqrt(y_mag);
    err = sqrt(err) / y_mag;
    if (err < args->in_convergence_criterion)
    {
        args->out_last_error = err;
        args->out_last_iteration = 0;
        return JMTX_RESULT_SUCCESS;
    }

    uint32_t iter_count = 0;
    for (;;)
    {
        jmtxz_matrix_crs_vector_multiply(mtx, p, Ap);
        _Complex double rQAp = 0;
        for (uint32_t i = 0; i < n; ++i)
        {
            rQAp += rQt[i] * Ap[i];
        }
        if (rQAp == 0)
        {
            break;
        }
        alpha = rho / rQAp;
        for (uint32_t i = 0; i < n; ++i)
        {
            x[i] = x[i] + alpha * p[i];
        }
        _Complex double sksk_dp = 0;
        for (uint32_t i = 0; i < n; ++i)
        {
            const _Complex double si = r[i] - alpha * Ap[i];
            s[i] = si;
            sksk_dp += conj(si) * si;
        }
        err = sqrt(sksk_dp) / y_mag;
        if (err < args->in_convergence_criterion)
        {
            break;
        }
        jmtxz_matrix_crs_vector_multiply(mtx, s, As);
        _Complex double sAs_dp = 0, sAAs_dp = 0;
        for (uint32_t i = 0; i < n; ++i)
        {
            sAAs_dp += conj(As[i]) * As[i];
            sAs_dp += conj(s[i]) * As[i];
        }
        if (sAs_dp == 0)
        {
            break;
        }
        omega = sAs_dp / sAAs_dp;
        for (uint32_t i = 0; i < n; ++i)
        {
            x[i] = x[i] + omega * s[i];
        }
        _Complex double rkrk_dp = 0;
        for (uint32_t i = 0; i < n; ++i)
        {
            const _Complex double ri = s[i] - omega * As[i];
            r[i] = ri;
            rkrk_dp += conj(ri) * ri;
        }
        err = sqrt(rkrk_dp) / y_mag;
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

        _Complex double rQrk_dp = 0;
        for (uint32_t i = 0; i < n; ++i)
        {
            rQrk_dp += rQt[i] * r[i];
        }
        const _Complex double beta = rQrk_dp / rho * alpha / omega;
        rho = rQrk_dp;
        if (rQrk_dp == 0)
        {
            break;
        }
        for (uint32_t i = 0; i < n; ++i)
        {
            p[i] = r[i] + beta * (p[i] - omega * Ap[i]);
        }
    }

    args->out_last_iteration = iter_count;
    args->out_last_error= err;
    if (!isfinite(err) || err > args->in_convergence_criterion)
    {
        return JMTX_RESULT_NOT_CONVERGED;
    }

    return JMTX_RESULT_SUCCESS;
}

static inline int check_vector_overlaps(const unsigned n, const size_t size, const void* ptrs[JMTX_ARRAY_ATTRIB(static const n)])
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
 *  Solves the linear problem A x = y for a general matrix A by using the relations used for Bi-CG, but does not
 *  explicitly solve the adjoint problem, instead computing values by computing results of polynomial relations for it.
 *  Stabilized method also computes these indirectly by using a polynomial with a lower condition number, giving better
 *  convergence behaviour.
 *
 *  This version of the funciton checks for appropriate matrix type and dimensions, as well as for memory not
 *  overlapping.
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
 * given number of iterations, other error codes in case of other errors
 */
jmtx_result jmtxzs_solve_iterative_bicgstab_crs(
        const jmtxz_matrix_crs* mtx, uint32_t n, const _Complex double y[JMTX_ARRAY_ATTRIB(restrict static n)], _Complex double x[JMTX_ARRAY_ATTRIB(restrict n)], _Complex double aux_vec1[JMTX_ARRAY_ATTRIB(restrict n)],
        _Complex double aux_vec2[JMTX_ARRAY_ATTRIB(restrict n)], _Complex double aux_vec3[JMTX_ARRAY_ATTRIB(restrict n)], _Complex double aux_vec4[JMTX_ARRAY_ATTRIB(restrict n)], _Complex double aux_vec5[JMTX_ARRAY_ATTRIB(restrict n)],
        _Complex double aux_vec6[JMTX_ARRAY_ATTRIB(restrict n)], jmtxd_solver_arguments* args)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTXZ_TYPE_CRS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (!x)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.rows != n || mtx->base.cols != n)
    {
        return JMTX_RESULT_BAD_MATRIX;
    }
    if (!isfinite(args->in_convergence_criterion))
    {
        return JMTX_RESULT_BAD_PARAM;
    }
    //  Check if any of the vectors overlap
    const void* memory_addresses[] = { x, y, aux_vec1, aux_vec2, aux_vec3, aux_vec4, aux_vec5, aux_vec6 };
    if (check_vector_overlaps(sizeof(memory_addresses) / sizeof(*memory_addresses), n * sizeof(*x), memory_addresses))
    {
        return JMTX_RESULT_BAD_PARAM;
    }

    return jmtxz_solve_iterative_bicgstab_crs(mtx, y, x, aux_vec1, aux_vec2, aux_vec3, aux_vec4, aux_vec5, aux_vec6, args);
}

/**
 *  Solves the linear problem A x = y for a general matrix A by using the relations used for Bi-CG, but does not
 *  explicitly solve the adjoint problem, instead computing values by computing results of polynomial relations for it.
 *  Stabilized method also computes these indirectly by using a polynomial with a lower condition number, giving better
 *  convergence behaviour.
 *
 *  This version of the funciton does not check if its inputs are valid and just assumes they are.
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
jmtx_result jmtxz_solve_iterative_bicgstab_cds(
        const jmtxz_matrix_cds* mtx, const _Complex double* restrict y, _Complex double* restrict x, _Complex double* restrict aux_vec1,
        _Complex double* restrict aux_vec2, _Complex double* restrict aux_vec3, _Complex double* restrict aux_vec4, _Complex double* restrict aux_vec5,
        _Complex double* restrict aux_vec6, jmtxd_solver_arguments* args)
{
    const uint32_t n = mtx->base.rows;

    _Complex double rho = 1, alpha = 1, omega = 1;

    _Complex double* const r = aux_vec1;
    _Complex double* const rQt = aux_vec2;
    _Complex double* const p = aux_vec3;
    _Complex double* const Ap = aux_vec4;
    _Complex double* const s = aux_vec5;
    _Complex double* const As = aux_vec6;

    double err = 0;

    jmtxz_matrix_cds_vector_multiply(mtx, x, r);
    double y_mag = 0;
    for (uint32_t i = 0; i < n; ++i)
    {
        r[i] = y[i] - r[i];
        rQt[i] = conj(r[i]);
        p[i] = r[i];
//        Ap[i] = 0;
        y_mag += conj(y[i]) * y[i];
        err += rQt[i] * r[i];
    }
    y_mag = sqrt(y_mag);
    err = sqrt(err) / y_mag;
    if (err < args->in_convergence_criterion)
    {
        args->out_last_error = err;
        args->out_last_iteration = 0;
        return JMTX_RESULT_SUCCESS;
    }

    uint32_t iter_count = 0;
    for (;;)
    {
        jmtxz_matrix_cds_vector_multiply(mtx, p, Ap);
        _Complex double rQAp = 0;
        for (uint32_t i = 0; i < n; ++i)
        {
            rQAp += rQt[i] * Ap[i];
        }
        alpha = rho / rQAp;
        for (uint32_t i = 0; i < n; ++i)
        {
            x[i] = x[i] + alpha * p[i];
        }
        _Complex double sksk_dp = 0;
        for (uint32_t i = 0; i < n; ++i)
        {
            const _Complex double si = r[i] - alpha * Ap[i];
            s[i] = si;
            sksk_dp += conj(si) * si;
        }
        err = sqrt(sksk_dp) / y_mag;
        if (err < args->in_convergence_criterion)
        {
            break;
        }
        jmtxz_matrix_cds_vector_multiply(mtx, s, As);
        _Complex double sAs_dp = 0, sAAs_dp = 0;
        for (uint32_t i = 0; i < n; ++i)
        {
            sAAs_dp += conj(As[i]) * As[i];
            sAs_dp += conj(s[i]) * As[i];
        }
        omega = sAs_dp / sAAs_dp;
        for (uint32_t i = 0; i < n; ++i)
        {
            x[i] = x[i] + omega * s[i];
        }
        _Complex double rkrk_dp = 0;
        for (uint32_t i = 0; i < n; ++i)
        {
            const _Complex double ri = s[i] - omega * As[i];
            r[i] = ri;
            rkrk_dp += conj(ri) * ri;
        }
        err = sqrt(rkrk_dp) / y_mag;
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

        _Complex double rQrk_dp = 0;
        for (uint32_t i = 0; i < n; ++i)
        {
            rQrk_dp += rQt[i] * r[i];
        }
        const _Complex double beta = rQrk_dp / rho * alpha / omega;
        rho = rQrk_dp;
        for (uint32_t i = 0; i < n; ++i)
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
 *  This version of the funciton checks for appropriate matrix type and dimensions, as well as for memory not
 *  overlapping.
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
 * given number of iterations, other error codes in case of other errors
 */
jmtx_result jmtxzs_solve_iterative_bicgstab_cds(
        const jmtxz_matrix_cds* mtx, uint32_t n, const _Complex double y[JMTX_ARRAY_ATTRIB(restrict static n)], _Complex double x[JMTX_ARRAY_ATTRIB(restrict n)], _Complex double aux_vec1[JMTX_ARRAY_ATTRIB(restrict n)],
        _Complex double aux_vec2[JMTX_ARRAY_ATTRIB(restrict n)], _Complex double aux_vec3[JMTX_ARRAY_ATTRIB(restrict n)], _Complex double aux_vec4[JMTX_ARRAY_ATTRIB(restrict n)], _Complex double aux_vec5[JMTX_ARRAY_ATTRIB(restrict n)],
        _Complex double aux_vec6[JMTX_ARRAY_ATTRIB(restrict n)], jmtxd_solver_arguments* args)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTXZ_TYPE_CDS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (!x)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.rows != n || mtx->base.cols != n)
    {
        return JMTX_RESULT_BAD_MATRIX;
    }
    if (!isfinite(args->in_convergence_criterion))
    {
        return JMTX_RESULT_BAD_PARAM;
    }
    //  Check if any of the vectors overlap
    const void* memory_addresses[] = { x, y, aux_vec1, aux_vec2, aux_vec3, aux_vec4, aux_vec5, aux_vec6 };
    if (check_vector_overlaps(sizeof(memory_addresses) / sizeof(*memory_addresses), n * sizeof(*x), memory_addresses))
    {
        return JMTX_RESULT_BAD_PARAM;
    }

    return jmtxz_solve_iterative_bicgstab_cds(mtx, y, x, aux_vec1, aux_vec2, aux_vec3, aux_vec4, aux_vec5, aux_vec6, args);
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
 *  This version of the funciton does not check if its inputs are valid and just assumes they are.
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
jmtx_result jmtxz_solve_iterative_pilubicgstab_crs(
        const jmtxz_matrix_crs* mtx, const jmtxz_matrix_crs* l, const jmtxz_matrix_crs* u, const _Complex double* restrict y,
        _Complex double* restrict x, _Complex double* restrict aux_vec1, _Complex double* restrict aux_vec2, _Complex double* restrict aux_vec3,
        _Complex double* restrict aux_vec4, _Complex double* restrict aux_vec5, _Complex double* restrict aux_vec6, _Complex double* restrict aux_vec7,
        _Complex double* restrict aux_vec8, jmtxd_solver_arguments* args)
{
    const uint32_t n = mtx->base.rows;

    _Complex double rho = 1, alpha = 1, omega = 1;

    _Complex double* const r = aux_vec1;
    _Complex double* const rQ = aux_vec2;
    _Complex double* const p = aux_vec3;
    _Complex double* const Ap = aux_vec4;
    _Complex double* const s = aux_vec5;
    _Complex double* const As = aux_vec6;
    _Complex double* const phat = aux_vec7;
    _Complex double* const shat = aux_vec8;

    double err = 0;

    jmtxz_matrix_crs_vector_multiply(mtx, x, r);
    double y_mag = 0;
    for (uint32_t i = 0; i < n; ++i)
    {
        r[i] = y[i] - r[i];
        rQ[i] = conj(r[i]);
        p[i] = r[i];
//        Ap[i] = 0;
        y_mag += conj(y[i]) * y[i];
        err += conj(r[i]) * r[i];
    }
    y_mag = sqrt(y_mag);
    err = sqrt(err) / y_mag;
    if (err < args->in_convergence_criterion)
    {
        args->out_last_error = err;
        args->out_last_iteration = 0;
        return JMTX_RESULT_SUCCESS;
    }

    uint32_t iter_count = 0;
    for (;;)
    {
        jmtxz_solve_direct_lu_crs(l, u, p, phat);
        jmtxz_matrix_crs_vector_multiply(mtx, phat, Ap);
        _Complex double rQAp = 0;
        for (uint32_t i = 0; i < n; ++i)
        {
            rQAp += rQ[i] * Ap[i];
        }
        alpha = rho / rQAp;
        for (uint32_t i = 0; i < n; ++i)
        {
            x[i] = x[i] + alpha * phat[i];
        }
        double sksk_dp = 0;
        for (uint32_t i = 0; i < n; ++i)
        {
            const _Complex double si = r[i] - alpha * Ap[i];
            s[i] = si;
            sksk_dp += conj(si) * si;
        }
        err = sqrt(sksk_dp) / y_mag;
        if (err < args->in_convergence_criterion)
        {
            break;
        }
        jmtxz_solve_direct_lu_crs(l, u, s, shat);
        jmtxz_matrix_crs_vector_multiply(mtx, shat, As);
        _Complex double sAs_dp = 0;
        double sAAs_dp = 0;
        for (uint32_t i = 0; i < n; ++i)
        {
            sAAs_dp += conj(As[i]) * As[i];
            sAs_dp  += conj(s[i]) * As[i];
        }
        omega = sAs_dp / sAAs_dp;
        for (uint32_t i = 0; i < n; ++i)
        {
            x[i] = x[i] + omega * shat[i];
        }
        double rkrk_dp = 0;
        for (uint32_t i = 0; i < n; ++i)
        {
            const _Complex double ri = s[i] - omega * As[i];
            r[i] = ri;
            rkrk_dp += conj(ri) * ri;
        }
        err = sqrt(rkrk_dp) / y_mag;
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

        _Complex double rQrk_dp = 0;
        for (uint32_t i = 0; i < n; ++i)
        {
            rQrk_dp += rQ[i] * r[i];
        }
        const double beta = rQrk_dp / rho * alpha / omega;
        rho = rQrk_dp;
        for (uint32_t i = 0; i < n; ++i)
        {
            p[i] = r[i] + beta * (p[i] - omega * Ap[i]);
        }
    }

    args->out_last_iteration = iter_count;
    args->out_last_error= err;
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
 *  This version of the funciton does not check if its inputs are valid and just assumes they are.
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
jmtx_result jmtxz_solve_iterative_pilubicgstab_crs_parallel(
        const jmtxz_matrix_crs* mtx, const jmtxz_matrix_crs* l, const jmtxz_matrix_crs* u, const _Complex double* restrict y,
        _Complex double* restrict x, _Complex double* restrict aux_vec1, _Complex double* restrict aux_vec2, _Complex double* restrict aux_vec3,
        _Complex double* restrict aux_vec4, _Complex double* restrict aux_vec5, _Complex double* restrict aux_vec6, _Complex double* restrict aux_vec7,
        _Complex double* restrict aux_vec8, jmtxd_solver_arguments* args)
{
    const uint32_t n = mtx->base.rows;

    _Complex double rho = 1, alpha = 1, omega = 1;

    _Complex double* const r = aux_vec1;
    _Complex double* const rQ = aux_vec2;
    _Complex double* const p = aux_vec3;
    _Complex double* const Ap = aux_vec4;
    _Complex double* const s = aux_vec5;
    _Complex double* const As = aux_vec6;
    _Complex double* const phat = aux_vec7;
    _Complex double* const shat = aux_vec8;

    _Complex double rQAp = 0;
    double sksk_dp = 0;
    double sAs_dp = 0, sAAs_dp = 0;
    double y_mag = 0;
    double rkrk_dp = 0;
    _Complex double rQrk_dp = 0;

    uint32_t iter_count = 0;

    double err = 0;
#pragma omp parallel shared(err, y_mag, r, rQ, p, Ap, s, As, phat, shat, rho, alpha, omega, n, iter_count, rQAp, sksk_dp, sAs_dp, sAAs_dp, rkrk_dp)
    {
        // jmtxz_matrix_crs_vector_multiply(mtx, x, r);
#pragma omp for reduction(+:err,y_mag)
        for (uint32_t i = 0; i < n; ++i)
        {
            const _Complex double yv = jmtxz_matrix_crs_vector_multiply_row(mtx, x, i);
            r[i] = y[i] - yv;
            rQ[i] = r[i];
            p[i] = r[i];
            //        Ap[i] = 0;
            y_mag += y[i] * conj(y[i]);
            err += conj(r[i]) * r[i];
        }
#pragma omp single
        {
            y_mag = sqrt(y_mag);
            err = sqrt(err) / y_mag;
        }

        for (;;)
        {
#pragma omp single
            {
                jmtxz_solve_direct_lu_crs(l, u, p, phat);
                rQAp = 0;
                sksk_dp = 0;
            }

            // jmtxz_matrix_crs_vector_multiply(mtx, phat, Ap);

#pragma omp for reduction(+:rQAp)
            for (uint32_t i = 0; i < n; ++i)
            {
                const _Complex double ap = jmtxz_matrix_crs_vector_multiply_row(mtx, phat, i);
                Ap[i] = ap;
                rQAp += conj(rQ[i]) * Ap[i];
            }

            if (rQAp == 0)
            {
                break;
            }
#pragma omp single
            alpha = rho / rQAp;
#pragma omp for
            for (uint32_t i = 0; i < n; ++i)
            {
                x[i] = x[i] + alpha * phat[i];
            }
#pragma omp for reduction(+:sksk_dp)
            for (uint32_t i = 0; i < n; ++i)
            {
                const _Complex double si = r[i] - alpha * Ap[i];
                s[i] = si;
                sksk_dp += conj(si) * si;
            }

#pragma omp single
            {
                err = sqrt(sksk_dp) / y_mag;
            }
            if (err < args->in_convergence_criterion)
            {
                break;
            }
#pragma omp single
            {
                jmtxz_solve_direct_lu_crs(l, u, s, shat);
                sAs_dp = 0, sAAs_dp = 0, rkrk_dp = 0, rQrk_dp = 0;
            }

            // jmtxz_matrix_crs_vector_multiply(mtx, shat, As);
#pragma omp for reduction(+:sAs_dp,sAAs_dp)
            for (uint32_t i = 0; i < n; ++i)
            {
                const _Complex double as = jmtxz_matrix_crs_vector_multiply_row(mtx, shat, i);
                As[i] = as;
                sAAs_dp += conj(As[i]) * As[i];
                sAs_dp += conj(s[i]) * As[i];
            }

            if (sAs_dp == 0)
            {
                break;
            }
#pragma omp single
            omega = sAs_dp / sAAs_dp;
#pragma omp for
            for (uint32_t i = 0; i < n; ++i)
            {
                x[i] = x[i] + omega * shat[i];
            }
#pragma omp for reduction(+:rkrk_dp)
            for (uint32_t i = 0; i < n; ++i)
            {
                const _Complex double ri = s[i] - omega * As[i];
                r[i] = ri;
                rkrk_dp += conj(ri) * ri;
            }
#pragma omp single
            {
                err = sqrt(rkrk_dp) / y_mag;
            }
            if (iter_count == args->in_max_iterations)
            {
                break;
            }
            if (args->opt_error_evolution)
            {
                args->opt_error_evolution[iter_count] = err;
            }
#pragma omp single
            iter_count += 1;
            if (err < args->in_convergence_criterion)
            {
                break;
            }

#pragma omp for reduction(+:rQrk_dp)
            for (uint32_t i = 0; i < n; ++i)
            {
                rQrk_dp += conj(rQ[i]) * r[i];
            }

            if (rQrk_dp == 0)
            {
                break;
            }

            const _Complex double beta = rQrk_dp / rho * alpha / omega;
            rho = rQrk_dp;
#pragma omp for
            for (uint32_t i = 0; i < n; ++i)
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


