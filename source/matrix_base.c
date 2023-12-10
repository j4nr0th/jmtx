#include <assert.h>
#include <stdlib.h>
#include "../include/jmtx/matrix_base.h"
#include "matrix_base_internal.h"

static const char* const jmtx_result_string_array[JMTX_RESULT_COUNT] =
        {
                [JMTX_RESULT_SUCCESS] = "Success",    //  All is good
                [JMTX_RESULT_NOT_CONVERGED] = "Iteration did not converge within the given time",      //  Iteration didn't converge
                [JMTX_RESULT_STAGNATED] = "Improvement between operations was too low", //  Iteration slowed down
                [JMTX_RESULT_BAD_ALLOC] = "Memory allocation did not succeed",          //  memory allocator being a cunt
                [JMTX_RESULT_BAD_THREAD] = "Thread creation failed",         //  Couldn't create a thread
                [JMTX_RESULT_WRONG_TYPE] = "Wrong matrix type was passed",         //  Matrix has the wrong storage type (e.g.: ccs instead of crs)
                [JMTX_RESULT_DIMS_MISMATCH] = "Dimensions of matrices/vectors did not match",      //  Matrix or a vector don't match
                [JMTX_RESULT_INDEX_OUT_OF_BOUNDS] = "Index was too large",//  Index was too large
                [JMTX_RESULT_BAD_PARAM] = "Parameter had an invalid value",          //  Parameter had a bad value
                [JMTX_RESULT_NULL_PARAM] = "Parameter was null",         //  Parameter was null
                [JMTX_RESULT_BAD_MATRIX] = "Matrix could not be handled by the iterative solver",   //  Matrix is fucked
                [JMTX_RESULT_UNARY_RETURN] = "Unary function returned non-zero",
                [JMTX_RESULT_NO_FILE] = "Could not open file for writing",
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
                [JMTX_TYPE_NONE] = "Invalid",                   //  Invalid, here to force one of other values to be specified
                //  Float
                [JMTX_TYPE_CRS] = "Compressed-row sparse",      //  Compressed row sparse
                [JMTX_TYPE_CCS] = "Compressed-column sparse",   //  Compressed column sparse
                [JMTX_TYPE_BRM] = "Band row-major",             //  Band row-major
                [JMTX_TYPE_CDS] = "Compressed-diagonal sparse", //  Compressed diagonal sparse
                //  Double
                [JMTXD_TYPE_CRS] = "Double precision compressed-row sparse",      //  Compressed row sparse
                [JMTXD_TYPE_CCS] = "Double precision compressed-column sparse",   //  Compressed column sparse
                [JMTXD_TYPE_BRM] = "Double precision band row-major",             //  Band row-major
                [JMTXD_TYPE_CDS] = "Double precision compressed-diagonal sparse", //  Compressed diagonal sparse
                //  Complex Float
                [JMTXC_TYPE_CRS] = "Complex Compressed-row sparse",      //  Compressed row sparse
                [JMTXC_TYPE_CCS] = "Complex Compressed-column sparse",   //  Compressed column sparse
                [JMTXC_TYPE_BRM] = "Complex Band row-major",             //  Band row-major
                [JMTXC_TYPE_CDS] = "Complex Compressed-diagonal sparse", //  Compressed diagonal sparse
                //  Double
                [JMTXZ_TYPE_CRS] = "Complex Double precision compressed-row sparse",      //  Compressed row sparse
                [JMTXZ_TYPE_CCS] = "Complex Double precision compressed-column sparse",   //  Compressed column sparse
                [JMTXZ_TYPE_BRM] = "Complex Double precision band row-major",             //  Band row-major
                [JMTXZ_TYPE_CDS] = "Complex Double precision compressed-diagonal sparse", //  Compressed diagonal sparse
        };

const char* jmtx_matrix_type_to_str(jmtx_matrix_type type)
{
    if (type >= JMTX_TYPE_NONE && type < JMTX_TYPE_COUNT)
    {
        return jmtx_matrix_type_string_array[type];
    }
    return "Unknown";
}

static const char* const funny_string = "Why are FEM engineers bad at deadlifting?\nTheir problem is in their weak form!\n";

static void* default_alloc(void* state, uint64_t size)
{
    assert(state == funny_string);
    void* const ptr = malloc(size);
#ifndef NDEBUG
    if (ptr)
    {
        memset(ptr, 0xCC, size);
    }
#endif
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
//    if (p_out) BEEF_FILL(p_out, new_size);    DO NOT BEEF THE NEW POINTER!
    return p_out;
}

const jmtx_allocator_callbacks JMTX_DEFAULT_ALLOCATOR_CALLBACKS =
        {
        .alloc = default_alloc,
        .free = default_free,
        .realloc = default_realloc,
        .state = (void*)funny_string,
        };

[[gnu::hot]]
uint32_t jmtx_internal_find_last_leq_value(uint32_t n_indices, const uint32_t p_indices[static n_indices], uint32_t value)
{
    uint32_t p = 0;
    uint32_t n = n_indices;
    uint32_t s = n - n / 2;
    const uint32_t* c = p_indices;
    while (s > 1)
    {
        if (c[p + s] > value)
        {
            n = s;
        }
        else
        {
            p += s;
            n -= s;
        }
        s = n - n / 2;
    }

    if (s != 0)
    {
        if (p == n_indices - 1)
        {
            return p;
        }
        assert(s == 1);
        if (c[p + s] > value)
        {
            return p;
        }
        else
        {
            return p + 1;
        }
    }


    return p;
}
