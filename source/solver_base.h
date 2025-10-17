//
// Created by jan on 21.11.2023.
//

#ifndef JMTX_SOLVER_BASE_H
#define JMTX_SOLVER_BASE_H

#include "common.h"

/**
 * Parameters used for any kind of iterative solution algorithm.
 */
struct JMTX_NAME_TYPED(solver_arguments_t)
{
    /**
     * Used to determine if the solution is close enough to the exact value. For linear solvers this is typically
     * implemented as ratio of magnitude of the residual vector divided by the magnitude of the forcing vector.
     */
    JMTX_REAL_T in_convergence_criterion;
    /**
     * Maximum number of iterations to perform. Solver will stop once this number is reached and indicate by returning
     * JMTX_RESULT_NOT_CONVERGED if convergence criterion was not met.
     */
    JMTX_INDEX_T in_max_iterations;

    /**
     * Receives the value of the error measure computed in the last iteration. If the algorithm converged it must be
     * less than convergence criterion. For linear solvers this is typically implemented as ratio of magnitude of the
     * residual vector divided by the magnitude of the forcing vector.
     */
    JMTX_REAL_T out_last_error;
    /**
     * Receives the number of the last iteration. Will be less or equal to the specified number of max iterations.
     */
    JMTX_INDEX_T out_last_iteration;

    /**
     * Pointer may be optionally specified to point to an array of the same length as the maximum number of iterations.
     * Each iteration, the error measure for computed for that iteration i will be written array[i].
     */
    JMTX_REAL_T *opt_error_evolution;
};
typedef struct JMTX_NAME_TYPED(solver_arguments_t) JMTX_NAME_TYPED(solver_arguments);

#endif // JMTX_SOLVER_BASE_H
