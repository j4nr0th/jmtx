#include <assert.h>
#include <stdlib.h>
#include "matrix_base.h"

static const char* const jmtx_result_string_array[JMTX_RESULT_COUNT] =
        {
                [JMTX_RESULT_SUCCESS] = "Success",    //  All is good
                [JMTX_RESULT_NOT_CONVERGED] = "Iteration did not converge within the given time",      //  Iteration didn't converge
                [JMTX_RESULT_BAD_ALLOC] = "Memory allocation did not succeed",          //  memory allocator being a cunt
                [JMTX_RESULT_BAD_THREAD] = "Thread creation failed",         //  Couldn't create a thread
                [JMTX_RESULT_WRONG_TYPE] = "Wrong matrix type was passed",         //  Matrix has the wrong storage type (e.g.: ccs instead of crs)
                [JMTX_RESULT_DIMS_MISMATCH] = "Dimensions of matrices/vectors did not match",      //  Matrix or a vector don't match
                [JMTX_RESULT_INDEX_OUT_OF_BOUNDS] = "Index was too large",//  Index was too large
                [JMTX_RESULT_BAD_PARAM] = "Parameter had an invalid value",          //  Parameter had a bad value
                [JMTX_RESULT_NULL_PARAM] = "Parameter was null",         //  Parameter was null
        };

const char* jmtx_result_to_str(jmtx_result res)
{
    if (res >= JMTX_RESULT_SUCCESS && res < JMTX_RESULT_COUNT)
    {
        return jmtx_result_string_array[res];
    }
    return "Unknown";
}

static const char* const jmtx_matrix_type_string_array[JMTX_TYPE_COUNT] =
        {
                [JMTX_TYPE_NONE] = "Invalid", //  Invalid, here to force one of other values to be specified
                [JMTX_TYPE_CRS] = "Compressed-row sparse",      //  Compressed row sparse
                [JMTX_TYPE_CCS] = "Compressed-column sparse",      //  Compressed column sparse
                [JMTX_TYPE_DENSE] = "Dense/full",
        };

const char* jmtx_matrix_type_to_str(jmtx_matrix_type type)
{
    if (type >= JMTX_TYPE_NONE && type < JMTX_TYPE_COUNT)
    {
        return jmtx_matrix_type_string_array[type];
    }
    return "Unknown";
}

#ifndef NDEBUG
static void BEEF_FILL(void* restrict const ptr, const size_t size)
{
    for (uint32_t* p = ptr; ((uintptr_t)p) < ((uintptr_t)ptr) + size; ++p)
    {
        *p = 0xDeadBeef;
    }
    memset((uint8_t*)ptr + (size - (size & 3)), 0, size & 3);
}

#else
#define BEEF_FILL(ptr, size) (void)0
#endif

static const char* const funny_string = "Why are FEM engineers bad at deadlifting?\nTheir problem is in their weak form!\n";

static void* default_alloc(void* state, uint64_t size)
{
    assert(state == funny_string);
    void* const ptr = malloc(size);
    if (ptr) BEEF_FILL(ptr, size);
    return ptr;
}

static void default_free(void* state, void* ptr)
{
    assert(state == funny_string);
    free(ptr);
}

static void* default_realloc(void* state, void* ptr, uint64_t new_size)
{
    assert(state == funny_string);
    void* const p_out = realloc(ptr, new_size);
    if (p_out) BEEF_FILL(p_out, new_size);
    return p_out;
}

const jmtx_allocator_callbacks JMTX_DEFAULT_ALLOCATOR_CALLBACKS =
        {
        .alloc = default_alloc,
        .free = default_free,
        .realloc = default_realloc,
        .state = (void*)funny_string,
        };