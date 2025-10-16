//
// Created by jan on 17.6.2022.
//

#include "../matrices/band_row_major_internal.h"
#include "../matrices/sparse_diagonal_compressed_internal.h"
#include "../matrices/sparse_row_compressed_internal.h"
#include "../../../include/jmtx/float/solvers/bicgstab_iteration.h"
#include "../../../include/jmtx/float/solvers/lu_solving.h"
#include <math.h>

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
jmtx_result jmtx_solve_iterative_bicgstab_crs(const jmtxf_matrix_crs *mtx, const float *restrict y, float *restrict x,
                                              float *restrict aux_vec1, float *restrict aux_vec2,
                                              float *restrict aux_vec3, float *restrict aux_vec4,
                                              float *restrict aux_vec5, float *restrict aux_vec6,
                                              jmtxf_solver_arguments *args)
{
    const uint32_t n = mtx->base.rows;

    float rho = 1, alpha = 1, omega = 1;

    float *const r = aux_vec1;
    float *const rQ = aux_vec2;
    float *const p = aux_vec3;
    float *const Ap = aux_vec4;
    float *const s = aux_vec5;
    float *const As = aux_vec6;

    float err = 0;

    jmtxf_matrix_crs_vector_multiply(mtx, x, r);
    float y_mag = 0;
    for (uint32_t i = 0; i < n; ++i)
    {
        r[i] = y[i] - r[i];
        rQ[i] = r[i];
        p[i] = r[i];
        //        Ap[i] = 0;
        y_mag += y[i] * y[i];
        err += r[i] * r[i];
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
        jmtxf_matrix_crs_vector_multiply(mtx, p, Ap);
        float rQAp = 0;
        for (uint32_t i = 0; i < n; ++i)
        {
            rQAp += rQ[i] * Ap[i];
        }
        alpha = rho / rQAp;
        for (uint32_t i = 0; i < n; ++i)
        {
            x[i] = x[i] + alpha * p[i];
        }
        float sksk_dp = 0;
        for (uint32_t i = 0; i < n; ++i)
        {
            const float si = r[i] - alpha * Ap[i];
            s[i] = si;
            sksk_dp += si * si;
        }
        err = sqrtf(sksk_dp) / y_mag;
        if (err < args->in_convergence_criterion)
        {
            break;
        }
        jmtxf_matrix_crs_vector_multiply(mtx, s, As);
        float sAs_dp = 0, sAAs_dp = 0;
        for (uint32_t i = 0; i < n; ++i)
        {
            sAAs_dp += As[i] * As[i];
            sAs_dp += s[i] * As[i];
        }
        omega = sAs_dp / sAAs_dp;
        for (uint32_t i = 0; i < n; ++i)
        {
            x[i] = x[i] + omega * s[i];
        }
        float rkrk_dp = 0;
        for (uint32_t i = 0; i < n; ++i)
        {
            const float ri = s[i] - omega * As[i];
            r[i] = ri;
            rkrk_dp += ri * ri;
        }
        err = sqrtf(rkrk_dp) / y_mag;
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

        float rQrk_dp = 0;
        for (uint32_t i = 0; i < n; ++i)
        {
            rQrk_dp += rQ[i] * r[i];
        }
        const float beta = rQrk_dp / rho * alpha / omega;
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

static inline int check_vector_overlaps(const unsigned n, const size_t size,
                                        const void *ptrs[JMTX_ARRAY_ATTRIB(static const n)])
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
 *  This version of the function checks for appropriate matrix type and dimensions, as well as for memory not
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
jmtx_result jmtxs_solve_iterative_bicgstab_crs(
    const jmtxf_matrix_crs *mtx, uint32_t n, const float y[JMTX_ARRAY_ATTRIB(restrict static n)],
    float x[JMTX_ARRAY_ATTRIB(restrict n)], float aux_vec1[JMTX_ARRAY_ATTRIB(restrict n)],
    float aux_vec2[JMTX_ARRAY_ATTRIB(restrict n)], float aux_vec3[JMTX_ARRAY_ATTRIB(restrict n)],
    float aux_vec4[JMTX_ARRAY_ATTRIB(restrict n)], float aux_vec5[JMTX_ARRAY_ATTRIB(restrict n)],
    float aux_vec6[JMTX_ARRAY_ATTRIB(restrict n)], jmtxf_solver_arguments *args)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTX_TYPE_CRS)
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
    const void *memory_addresses[] = {x, y, aux_vec1, aux_vec2, aux_vec3, aux_vec4, aux_vec5, aux_vec6};
    if (check_vector_overlaps(sizeof(memory_addresses) / sizeof(*memory_addresses), n * sizeof(*x), memory_addresses))
    {
        return JMTX_RESULT_BAD_PARAM;
    }

    return jmtx_solve_iterative_bicgstab_crs(mtx, y, x, aux_vec1, aux_vec2, aux_vec3, aux_vec4, aux_vec5, aux_vec6,
                                             args);
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
jmtx_result jmtx_solve_iterative_bicgstab_cds(const jmtx_matrix_cds *mtx, const float *restrict y, float *restrict x,
                                              float *restrict aux_vec1, float *restrict aux_vec2,
                                              float *restrict aux_vec3, float *restrict aux_vec4,
                                              float *restrict aux_vec5, float *restrict aux_vec6,
                                              jmtxf_solver_arguments *args)
{
    const uint32_t n = mtx->base.rows;

    float rho = 1, alpha = 1, omega = 1;

    float *const r = aux_vec1;
    float *const rQ = aux_vec2;
    float *const p = aux_vec3;
    float *const Ap = aux_vec4;
    float *const s = aux_vec5;
    float *const As = aux_vec6;

    float err = 0;

    jmtx_matrix_cds_vector_multiply(mtx, x, r);
    float y_mag = 0;
    for (uint32_t i = 0; i < n; ++i)
    {
        r[i] = y[i] - r[i];
        rQ[i] = r[i];
        p[i] = r[i];
        //        Ap[i] = 0;
        y_mag += y[i] * y[i];
        err += r[i] * r[i];
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
        jmtx_matrix_cds_vector_multiply(mtx, p, Ap);
        float rQAp = 0;
        for (uint32_t i = 0; i < n; ++i)
        {
            rQAp += rQ[i] * Ap[i];
        }
        alpha = rho / rQAp;
        for (uint32_t i = 0; i < n; ++i)
        {
            x[i] = x[i] + alpha * p[i];
        }
        float sksk_dp = 0;
        for (uint32_t i = 0; i < n; ++i)
        {
            const float si = r[i] - alpha * Ap[i];
            s[i] = si;
            sksk_dp += si * si;
        }
        err = sqrtf(sksk_dp) / y_mag;
        if (err < args->in_convergence_criterion)
        {
            break;
        }
        jmtx_matrix_cds_vector_multiply(mtx, s, As);
        float sAs_dp = 0, sAAs_dp = 0;
        for (uint32_t i = 0; i < n; ++i)
        {
            sAAs_dp += As[i] * As[i];
            sAs_dp += s[i] * As[i];
        }
        omega = sAs_dp / sAAs_dp;
        for (uint32_t i = 0; i < n; ++i)
        {
            x[i] = x[i] + omega * s[i];
        }
        float rkrk_dp = 0;
        for (uint32_t i = 0; i < n; ++i)
        {
            const float ri = s[i] - omega * As[i];
            r[i] = ri;
            rkrk_dp += ri * ri;
        }
        err = sqrtf(rkrk_dp) / y_mag;
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

        float rQrk_dp = 0;
        for (uint32_t i = 0; i < n; ++i)
        {
            rQrk_dp += rQ[i] * r[i];
        }
        const float beta = rQrk_dp / rho * alpha / omega;
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
 *  This version of the function checks for appropriate matrix type and dimensions, as well as for memory not
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
jmtx_result jmtxs_solve_iterative_bicgstab_cds(
    const jmtx_matrix_cds *mtx, uint32_t n, const float y[JMTX_ARRAY_ATTRIB(restrict static n)],
    float x[JMTX_ARRAY_ATTRIB(restrict n)], float aux_vec1[JMTX_ARRAY_ATTRIB(restrict n)],
    float aux_vec2[JMTX_ARRAY_ATTRIB(restrict n)], float aux_vec3[JMTX_ARRAY_ATTRIB(restrict n)],
    float aux_vec4[JMTX_ARRAY_ATTRIB(restrict n)], float aux_vec5[JMTX_ARRAY_ATTRIB(restrict n)],
    float aux_vec6[JMTX_ARRAY_ATTRIB(restrict n)], jmtxf_solver_arguments *args)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTX_TYPE_CDS)
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
    const void *memory_addresses[] = {x, y, aux_vec1, aux_vec2, aux_vec3, aux_vec4, aux_vec5, aux_vec6};
    if (check_vector_overlaps(sizeof(memory_addresses) / sizeof(*memory_addresses), n * sizeof(*x), memory_addresses))
    {
        return JMTX_RESULT_BAD_PARAM;
    }

    return jmtx_solve_iterative_bicgstab_cds(mtx, y, x, aux_vec1, aux_vec2, aux_vec3, aux_vec4, aux_vec5, aux_vec6,
                                             args);
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
jmtx_result jmtx_solve_iterative_pilubicgstab_crs(const jmtxf_matrix_crs *mtx, const jmtxf_matrix_crs *l,
                                                  const jmtxf_matrix_crs *u, const float *restrict y, float *restrict x,
                                                  float *restrict aux_vec1, float *restrict aux_vec2,
                                                  float *restrict aux_vec3, float *restrict aux_vec4,
                                                  float *restrict aux_vec5, float *restrict aux_vec6,
                                                  float *restrict aux_vec7, float *restrict aux_vec8,
                                                  jmtxf_solver_arguments *args)
{
    const uint32_t n = mtx->base.rows;

    float rho = 1, alpha = 1, omega = 1;

    float *const r = aux_vec1;
    float *const rQ = aux_vec2;
    float *const p = aux_vec3;
    float *const Ap = aux_vec4;
    float *const s = aux_vec5;
    float *const As = aux_vec6;
    float *const p_hat = aux_vec7;
    float *const s_hat = aux_vec8;

    float err = 0;

    jmtxf_matrix_crs_vector_multiply(mtx, x, r);
    float y_mag = 0;
    for (uint32_t i = 0; i < n; ++i)
    {
        r[i] = y[i] - r[i];
        rQ[i] = r[i];
        p[i] = r[i];
        //        Ap[i] = 0;
        y_mag += y[i] * y[i];
        err += r[i] * r[i];
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
        jmtx_solve_direct_lu_crs(l, u, p, p_hat);
        jmtxf_matrix_crs_vector_multiply(mtx, p_hat, Ap);
        float rQAp = 0;
        for (uint32_t i = 0; i < n; ++i)
        {
            rQAp += rQ[i] * Ap[i];
        }
        if (rQAp == 0)
        {
            break;
        }
        alpha = rho / rQAp;
        for (uint32_t i = 0; i < n; ++i)
        {
            x[i] = x[i] + alpha * p_hat[i];
        }
        float sksk_dp = 0;
        for (uint32_t i = 0; i < n; ++i)
        {
            const float si = r[i] - alpha * Ap[i];
            s[i] = si;
            sksk_dp += si * si;
        }
        err = sqrtf(sksk_dp) / y_mag;
        if (err < args->in_convergence_criterion)
        {
            break;
        }
        jmtx_solve_direct_lu_crs(l, u, s, s_hat);
        jmtxf_matrix_crs_vector_multiply(mtx, s_hat, As);
        float sAs_dp = 0, sAAs_dp = 0;
        for (uint32_t i = 0; i < n; ++i)
        {
            sAAs_dp += As[i] * As[i];
            sAs_dp += s[i] * As[i];
        }
        if (sAs_dp == 0)
        {
            break;
        }
        omega = sAs_dp / sAAs_dp;
        for (uint32_t i = 0; i < n; ++i)
        {
            x[i] = x[i] + omega * s_hat[i];
        }
        float rkrk_dp = 0;
        for (uint32_t i = 0; i < n; ++i)
        {
            const float ri = s[i] - omega * As[i];
            r[i] = ri;
            rkrk_dp += ri * ri;
        }
        err = sqrtf(rkrk_dp) / y_mag;
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

        float rQrk_dp = 0;
        for (uint32_t i = 0; i < n; ++i)
        {
            rQrk_dp += rQ[i] * r[i];
        }
        const float beta = rQrk_dp / rho * alpha / omega;
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
jmtx_result jmtx_solve_iterative_pilubicgstab_crs_parallel(const jmtxf_matrix_crs *mtx, const jmtxf_matrix_crs *l,
                                                           const jmtxf_matrix_crs *u, const float *restrict y,
                                                           float *restrict x, float *restrict aux_vec1,
                                                           float *restrict aux_vec2, float *restrict aux_vec3,
                                                           float *restrict aux_vec4, float *restrict aux_vec5,
                                                           float *restrict aux_vec6, float *restrict aux_vec7,
                                                           float *restrict aux_vec8, jmtxf_solver_arguments *args)
{
    const uint32_t n = mtx->base.rows;

    float rho = 1, alpha = 1, omega = 1;

    float *const r = aux_vec1;
    float *const rQ = aux_vec2;
    float *const p = aux_vec3;
    float *const Ap = aux_vec4;
    float *const s = aux_vec5;
    float *const As = aux_vec6;
    float *const p_hat = aux_vec7;
    float *const s_hat = aux_vec8;

    float rQAp = 0;
    float sksk_dp = 0;
    float sAs_dp = 0, sAAs_dp = 0;
    float y_mag = 0;
    float rkrk_dp = 0;
    float rQrk_dp = 0;

    uint32_t iter_count = 0;

    float err = 0;
#pragma omp parallel shared(err, y_mag, r, rQ, p, Ap, s, As, p_hat, s_hat, rho, alpha, omega, n, iter_count, rQAp,     \
                                sksk_dp, sAs_dp, sAAs_dp, rkrk_dp)
    {
        // jmtxf_matrix_crs_vector_multiply(mtx, x, r);
#pragma omp for reduction(+ : err, y_mag)
        for (uint32_t i = 0; i < n; ++i)
        {
            const float yv = jmtxf_matrix_crs_vector_multiply_row(mtx, x, i);
            r[i] = y[i] - yv;
            rQ[i] = r[i];
            p[i] = r[i];
            //        Ap[i] = 0;
            y_mag += y[i] * y[i];
            err += r[i] * r[i];
        }
#pragma omp single
        {
            y_mag = sqrtf(y_mag);
            err = sqrtf(err) / y_mag;
        }

        for (;;)
        {
#pragma omp single
            {
                jmtx_solve_direct_lu_crs(l, u, p, p_hat);
                rQAp = 0;
                sksk_dp = 0;
            }

            // jmtxf_matrix_crs_vector_multiply(mtx, p_hat, Ap);

#pragma omp for reduction(+ : rQAp)
            for (uint32_t i = 0; i < n; ++i)
            {
                const float ap = jmtxf_matrix_crs_vector_multiply_row(mtx, p_hat, i);
                Ap[i] = ap;
                rQAp += rQ[i] * Ap[i];
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
                x[i] = x[i] + alpha * p_hat[i];
            }
#pragma omp for reduction(+ : sksk_dp)
            for (uint32_t i = 0; i < n; ++i)
            {
                const float si = r[i] - alpha * Ap[i];
                s[i] = si;
                sksk_dp += si * si;
            }

#pragma omp single
            {
                err = sqrtf(sksk_dp) / y_mag;
            }
            if (err < args->in_convergence_criterion)
            {
                break;
            }
#pragma omp single
            {
                jmtx_solve_direct_lu_crs(l, u, s, s_hat);
                sAs_dp = 0, sAAs_dp = 0, rkrk_dp = 0, rQrk_dp = 0;
            }

            // jmtxf_matrix_crs_vector_multiply(mtx, s_hat, As);
#pragma omp for reduction(+ : sAs_dp, sAAs_dp)
            for (uint32_t i = 0; i < n; ++i)
            {
                const float as = jmtxf_matrix_crs_vector_multiply_row(mtx, s_hat, i);
                As[i] = as;
                sAAs_dp += As[i] * As[i];
                sAs_dp += s[i] * As[i];
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
                x[i] = x[i] + omega * s_hat[i];
            }
#pragma omp for reduction(+ : rkrk_dp)
            for (uint32_t i = 0; i < n; ++i)
            {
                const float ri = s[i] - omega * As[i];
                r[i] = ri;
                rkrk_dp += ri * ri;
            }
#pragma omp single
            {
                err = sqrtf(rkrk_dp) / y_mag;
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

#pragma omp for reduction(+ : rQrk_dp)
            for (uint32_t i = 0; i < n; ++i)
            {
                rQrk_dp += rQ[i] * r[i];
            }

            if (rQrk_dp == 0)
            {
                break;
            }

            const float beta = rQrk_dp / rho * alpha / omega;
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
