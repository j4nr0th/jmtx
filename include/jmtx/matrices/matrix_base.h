#ifndef JMTX_MATRIX_BASE_H
#define JMTX_MATRIX_BASE_H
#ifndef JMTX_COMMON_H
    #include "../common.h"
#endif
#include <stdint.h>
#include <stddef.h>
#include <string.h>



enum jmtx_result_enum
{
    JMTX_RESULT_SUCCESS         = 0,    //  All is good

    JMTX_RESULT_NOT_CONVERGED,      //  Iteration didn't converge
    JMTX_RESULT_STAGNATED,          //  Iterations stagnated
    JMTX_RESULT_BAD_ALLOC,          //  memory allocator being a cunt
    JMTX_RESULT_BAD_THREAD,         //  Couldn't create a thread
    JMTX_RESULT_WRONG_TYPE,         //  Matrix has the wrong storage type (e.g.: ccs instead of crs)
    JMTX_RESULT_DIMS_MISMATCH,      //  Matrix or a vector don't match
    JMTX_RESULT_INDEX_OUT_OF_BOUNDS,//  Index was too large
    JMTX_RESULT_BAD_PARAM,          //  Parameter had a bad value
    JMTX_RESULT_NULL_PARAM,         //  Parameter was null

    JMTX_RESULT_BAD_MATRIX,         //  Matrix is fucked

    JMTX_RESULT_UNARY_RETURN,       //  Unary function returned non-zero

    JMTX_RESULT_NO_FILE,            //  Could not open file for writing

    JMTX_RESULT_COUNT,
};
typedef enum jmtx_result_enum jmtx_result;

const char* jmtx_result_to_str(jmtx_result res);

typedef struct jmtx_allocator_callbacks_struct jmtx_allocator_callbacks;
struct jmtx_allocator_callbacks_struct
{
    void* (*alloc)(void* state, uint64_t size);
    void (*free)(void* state, void* ptr);
    void* (*realloc)(void* state, void* ptr, uint64_t new_size);
    void* state;
};

#endif // !JMTX_MATRIX_BASE_H
