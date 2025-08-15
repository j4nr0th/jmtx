#ifndef JMTX_MATRIX_BASE_H
#define JMTX_MATRIX_BASE_H
#define JMTX_MATRIX_BASE_H
#ifndef JMTX_COMMON_H
#    include "common.h"
#endif
#include <stddef.h>
#include <stdint.h>
#include <string.h>

enum jmtx_result_enum
{
    JMTX_RESULT_SUCCESS = 0, //  All is good

    JMTX_RESULT_NOT_CONVERGED,       //  Iteration didn't converge
    JMTX_RESULT_STAGNATED,           //  Iterations stagnated
    JMTX_RESULT_BAD_ALLOC,           //  memory allocator being a cunt
    JMTX_RESULT_BAD_THREAD,          //  Couldn't create a thread
    JMTX_RESULT_WRONG_TYPE,          //  Matrix has the wrong storage type (e.g.: ccs instead of crs)
    JMTX_RESULT_DIMS_MISMATCH,       //  Matrix or a vector don't match
    JMTX_RESULT_INDEX_OUT_OF_BOUNDS, //  Index was too large
    JMTX_RESULT_BAD_PARAM,           //  Parameter had a bad value
    JMTX_RESULT_NULL_PARAM,          //  Parameter was null

    JMTX_RESULT_BAD_MATRIX, //  Matrix is fucked

    JMTX_RESULT_UNARY_RETURN, //  Unary function returned non-zero

    JMTX_RESULT_NO_FILE, //  Could not open file for writing

    JMTX_RESULT_COUNT,
};
/**
 * Return type of functions that check parameters as well as functions which can fail when parameters are correct (like
 * due to memory allocation failing or solution not converging). Calling "jmtx_result_to_str" returns the string which
 * explains the meaning of each value.
 */
typedef enum jmtx_result_enum jmtx_result;

/**
 * Returns a string with more precise explanation of the passed jmtx_result value
 * @param res value from the jmtx_result enum for which to return the explanation for
 * @return statically allocated constant string which contains the explanation of the value of the passed value
 */
const char *jmtx_result_to_str(jmtx_result res);

/**
 * Struct which holds user supplied callbacks which are used by some matrix routines for
 * allocating/deallocating/reallocating internal memory. Except for the "state", all other members may not be non-null.
 */
typedef struct jmtx_allocator_callbacks_struct jmtx_allocator_callbacks;
struct jmtx_allocator_callbacks_struct
{
    /**
     * Callback which will be used to allocate memory
     * @param state member jmtx_allocator_callbacks::state which will always be passed to all callbacks in this struct
     * @param size minimum required size of the memory to allocate
     * @return pointer to a valid memory address, which can be used for at least size bytes, or NULL in case it was not
     * possible to allocate the memory
     */
    void *(*alloc)(void *state, uint64_t size);

    /**
     * Callback which will be used to deallocate memory
     * @param state member jmtx_allocator_callbacks::state which will always be passed to all callbacks in this struct
     * @param ptr pointer to the memory address previously allocated by the jmtx_allocator_callbacks::alloc function in
     * the same struct
     */
    void (*free)(void *state, void *ptr);

    /**
     * Callback which will be used to reallocate memory
     * @param state member jmtx_allocator_callbacks::state which will always be passed to all callbacks in this struct
     * @param ptr NULL or a pointer to a valid memory adderss previously allocated with jmtx_allocator_callbacks::alloc
     * @param new_size minimum required size of the memory after reallocation
     * @return pointer to a valid memory address, which can be used for at least new_size bytes, which will be the same
     * if ptr was non-null, or NULL in case it was not possible to reallocate the memory
     */
    void *(*realloc)(void *state, void *ptr, uint64_t new_size);

    /**
     * Pointer which is always passed as the first argument to all calls to the other functions specified in this struct
     */
    void *state;
};

#endif // !JMTX_MATRIX_BASE_H
