// Automatically generated from source/float/solvers/bicgstab_iteration.c on Sun Dec 17 15:54:30 2023
//
// Created by jan on 17.6.2022.
//

#include <math.h>
#include "../matrices/sparse_row_compressed_internal.h"
#include "../matrices/sparse_diagonal_compressed_internal.h"
#include "../matrices/band_row_major_internal.h"
#include "../../../include/jmtx/cfloat/solvers/bicgstab_iteration.h"
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
jmtx_result jmtxc_bicgstab_crs(
        const jmtxc_matrix_crs* mtx, const _Complex float* restrict y, _Complex float* restrict x, _Complex float* restrict aux_vec1,
        _Complex float* restrict aux_vec2, _Complex float* restrict aux_vec3, _Complex float* restrict aux_vec4, _Complex float* restrict aux_vec5,
        _Complex float* restrict aux_vec6, jmtx_solver_arguments* args)
{
    const uint32_t n = mtx->base.rows;

    _Complex float rho = 1, alpha = 1, omega = 1;

    _Complex float* const r = aux_vec1;
    _Complex float* const rQt = aux_vec2;
    _Complex float* const p = aux_vec3;
    _Complex float* const Ap = aux_vec4;
    _Complex float* const s = aux_vec5;
    _Complex float* const As = aux_vec6;

    float err = 0;

    jmtxc_matrix_crs_vector_multiply(mtx, x, r);
    float y_mag = 0;
    for (uint32_t i = 0; i < n; ++i)
    {
        r[i] = y[i] - r[i];
        rQt[i] = conjf(r[i]);
        p[i] = r[i];
//        Ap[i] = 0;
        y_mag += conjf(y[i]) * y[i];
        err += rQt[i] * r[i];
    }
    y_mag = sqrtf(y_mag);
    err = sqrtf(err) / y_mag;
    if (err < args->in_convergence_criterion)
    {
        args->out_last_error = err;
        args->out_last_iteration = 0;
        return JMTX_RESULT_SUCCESS;
    }

    uint32_t iter_count = 0;
    for (;;)
    {
        jmtxc_matrix_crs_vector_multiply(mtx, p, Ap);
        _Complex float rQAp = 0;
        for (uint32_t i = 0; i < n; ++i)
        {
            rQAp += rQt[i] * Ap[i];
        }
        alpha = rho / rQAp;
        for (uint32_t i = 0; i < n; ++i)
        {
            x[i] = x[i] + alpha * p[i];
        }
        _Complex float sksk_dp = 0;
        for (uint32_t i = 0; i < n; ++i)
        {
            const _Complex float si = r[i] - alpha * Ap[i];
            s[i] = si;
            sksk_dp += conjf(si) * si;
        }
        err = sqrtf(sksk_dp) / y_mag;
        if (err < args->in_convergence_criterion)
        {
            break;
        }
        jmtxc_matrix_crs_vector_multiply(mtx, s, As);
        _Complex float sAs_dp = 0, sAAs_dp = 0;
        for (uint32_t i = 0; i < n; ++i)
        {
            sAAs_dp += conjf(As[i]) * As[i];
            sAs_dp += conjf(s[i]) * As[i];
        }
        omega = sAs_dp / sAAs_dp;
        for (uint32_t i = 0; i < n; ++i)
        {
            x[i] = x[i] + omega * s[i];
        }
        _Complex float rkrk_dp = 0;
        for (uint32_t i = 0; i < n; ++i)
        {
            const _Complex float ri = s[i] - omega * As[i];
            r[i] = ri;
            rkrk_dp += conjf(ri) * ri;
        }
        err = sqrtf(rkrk_dp) / y_mag;
        if (args->opt_error_evolution)
        {
            args->opt_error_evolution[iter_count] = err;
        }
        if (iter_count == args->in_max_iterations)
        {
            break;
        }
        iter_count += 1;
        if (err < args->in_convergence_criterion)
        {
            break;
        }

        _Complex float rQrk_dp = 0;
        for (uint32_t i = 0; i < n; ++i)
        {
            rQrk_dp += rQt[i] * r[i];
        }
        const _Complex float beta = rQrk_dp / rho * alpha / omega;
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

static inline int check_vector_overlaps(const unsigned n, const size_t size, const void* ptrs[static const n])
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
jmtx_result jmtxcs_bicgstab_crs(
        const jmtxc_matrix_crs* mtx, uint32_t n, const _Complex float y[restrict static n], _Complex float x[restrict n], _Complex float aux_vec1[restrict n],
        _Complex float aux_vec2[restrict n], _Complex float aux_vec3[restrict n], _Complex float aux_vec4[restrict n], _Complex float aux_vec5[restrict n],
        _Complex float aux_vec6[restrict n], jmtx_solver_arguments* args)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTXC_TYPE_CRS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (!y)
    {
        return JMTX_RESULT_NULL_PARAM;
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

    return jmtxc_bicgstab_crs(mtx, y, x, aux_vec1, aux_vec2, aux_vec3, aux_vec4, aux_vec5, aux_vec6, args);
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
jmtx_result jmtxc_bicgstab_cds(
        const jmtxc_matrix_cds* mtx, const _Complex float* restrict y, _Complex float* restrict x, _Complex float* restrict aux_vec1,
        _Complex float* restrict aux_vec2, _Complex float* restrict aux_vec3, _Complex float* restrict aux_vec4, _Complex float* restrict aux_vec5,
        _Complex float* restrict aux_vec6, jmtx_solver_arguments* args)
{
    const uint32_t n = mtx->base.rows;

    _Complex float rho = 1, alpha = 1, omega = 1;

    _Complex float* const r = aux_vec1;
    _Complex float* const rQt = aux_vec2;
    _Complex float* const p = aux_vec3;
    _Complex float* const Ap = aux_vec4;
    _Complex float* const s = aux_vec5;
    _Complex float* const As = aux_vec6;

    float err = 0;

    jmtxc_matrix_cds_vector_multiply(mtx, x, r);
    float y_mag = 0;
    for (uint32_t i = 0; i < n; ++i)
    {
        r[i] = y[i] - r[i];
        rQt[i] = conjf(r[i]);
        p[i] = r[i];
//        Ap[i] = 0;
        y_mag += conjf(y[i]) * y[i];
        err += rQt[i] * r[i];
    }
    y_mag = sqrtf(y_mag);
    err = sqrtf(err) / y_mag;
    if (err < args->in_convergence_criterion)
    {
        args->out_last_error = err;
        args->out_last_iteration = 0;
        return JMTX_RESULT_SUCCESS;
    }

    uint32_t iter_count = 0;
    for (;;)
    {
        jmtxc_matrix_cds_vector_multiply(mtx, p, Ap);
        _Complex float rQAp = 0;
        for (uint32_t i = 0; i < n; ++i)
        {
            rQAp += rQt[i] * Ap[i];
        }
        alpha = rho / rQAp;
        for (uint32_t i = 0; i < n; ++i)
        {
            x[i] = x[i] + alpha * p[i];
        }
        _Complex float sksk_dp = 0;
        for (uint32_t i = 0; i < n; ++i)
        {
            const _Complex float si = r[i] - alpha * Ap[i];
            s[i] = si;
            sksk_dp += conjf(si) * si;
        }
        err = sqrtf(sksk_dp) / y_mag;
        if (err < args->in_convergence_criterion)
        {
            break;
        }
        jmtxc_matrix_cds_vector_multiply(mtx, s, As);
        _Complex float sAs_dp = 0, sAAs_dp = 0;
        for (uint32_t i = 0; i < n; ++i)
        {
            sAAs_dp += conjf(As[i]) * As[i];
            sAs_dp += conjf(s[i]) * As[i];
        }
        omega = sAs_dp / sAAs_dp;
        for (uint32_t i = 0; i < n; ++i)
        {
            x[i] = x[i] + omega * s[i];
        }
        _Complex float rkrk_dp = 0;
        for (uint32_t i = 0; i < n; ++i)
        {
            const _Complex float ri = s[i] - omega * As[i];
            r[i] = ri;
            rkrk_dp += conjf(ri) * ri;
        }
        err = sqrtf(rkrk_dp) / y_mag;
        if (args->opt_error_evolution)
        {
            args->opt_error_evolution[iter_count] = err;
        }
        if (iter_count == args->in_max_iterations)
        {
            break;
        }
        iter_count += 1;
        if (err < args->in_convergence_criterion)
        {
            break;
        }

        _Complex float rQrk_dp = 0;
        for (uint32_t i = 0; i < n; ++i)
        {
            rQrk_dp += rQt[i] * r[i];
        }
        const _Complex float beta = rQrk_dp / rho * alpha / omega;
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
jmtx_result jmtxcs_bicgstab_cds(
        const jmtxc_matrix_cds* mtx, uint32_t n, const _Complex float y[restrict static n], _Complex float x[restrict n], _Complex float aux_vec1[restrict n],
        _Complex float aux_vec2[restrict n], _Complex float aux_vec3[restrict n], _Complex float aux_vec4[restrict n], _Complex float aux_vec5[restrict n],
        _Complex float aux_vec6[restrict n], jmtx_solver_arguments* args)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTXC_TYPE_CDS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (!y)
    {
        return JMTX_RESULT_NULL_PARAM;
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

    return jmtxc_bicgstab_cds(mtx, y, x, aux_vec1, aux_vec2, aux_vec3, aux_vec4, aux_vec5, aux_vec6, args);
}