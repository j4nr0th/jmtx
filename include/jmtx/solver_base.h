//
// Created by jan on 21.11.2023.
//

#ifndef JMTX_SOLVER_BASE_H
#define JMTX_SOLVER_BASE_H

#ifndef JMTX_COMMON_H
#    include "common.h"
#endif

/**
 * Parameters used for any kind of iterative solution algorithm.
 */
struct jmtxf_solver_arguments_t
{
    /**
     * Used to determine if the solution is close enough to the exact value. For linear solvers this is typically
     * implemented as ratio of magnitude of the residual vector divided by the magnitude of the forcing vector.
     */
    float in_convergence_criterion;
    /**
     * Maximum number of iterations to perform. Solver will stop once this number is reached and indicate by returning
     * JMTX_RESULT_NOT_CONVERGED if convergence criterion was not met.
     */
    uint32_t in_max_iterations;

    /**
     * Receives the value of the error measure computed in the last iteration. If the algorithm converged it must be
     * less than convergence criterion. For linear solvers this is typically implemented as ratio of magnitude of the
     * residual vector divided by the magnitude of the forcing vector.
     */
    float out_last_error;
    /**
     * Receives the number of the last iteration. Will be less or equal to the specified number of max iterations.
     */
    uint32_t out_last_iteration;

    /**
     * Pointer may be optionally specified to point to an array of the same length as the maximum number of iterations.
     * Each iteration, the error measure for computed for that iteration i will be written array[i].
     */
    float *opt_error_evolution;
};
typedef struct jmtxf_solver_arguments_t jmtxf_solver_arguments;

/**
 * Parameters used for any kind of iterative solution algorithm.
 */
struct jmtxd_solver_arguments_t
{
    /**
     * Used to determine if the solution is close enough to the exact value. For linear solvers this is typically
     * implemented as ratio of magnitude of the residual vector divided by the magnitude of the forcing vector.
     */
    double in_convergence_criterion;
    /**
     * Maximum number of iterations to perform. Solver will stop once this number is reached and indicate by returning
     * JMTX_RESULT_NOT_CONVERGED if convergence criterion was not met.
     */
    uint32_t in_max_iterations;

    /**
     * Receives the value of the error measure computed in the last iteration. If the algorithm converged it must be
     * less than convergence criterion. For linear solvers this is typically implemented as ratio of magnitude of the
     * residual vector divided by the magnitude of the forcing vector.
     */
    double out_last_error;
    /**
     * Receives the number of the last iteration. Will be less or equal to the specified number of max iterations.
     */
    uint32_t out_last_iteration;

    /**
     * Pointer may be optionally specified to point to an array of the same length as the maximum number of iterations.
     * Each iteration, the error measure for computed for that iteration i will be written array[i].
     */
    double *opt_error_evolution;
};
typedef struct jmtxd_solver_arguments_t jmtxd_solver_arguments;

#endif // JMTX_SOLVER_BASE_H
