#ifndef JMTX_MATRIX_BASE_H
#define JMTX_MATRIX_BASE_H

#include <stdint.h>
#include <stddef.h>
#include <string.h>

#ifdef __GNUC__
#define JMTX_GCC_ONLY(x) x
#else
#define JMTX_GCC_ONLY(x)
#endif


enum jmtx_result_enum
{
    JMTX_RESULT_SUCCESS         = 0,    //  All is good

    JMTX_RESULT_NOT_CONVERGED,      //  Iteration didn't converge
    JMTX_RESULT_BAD_ALLOC,          //  memory allocator being a cunt
    JMTX_RESULT_BAD_THREAD,         //  Couldn't create a thread
    JMTX_RESULT_WRONG_TYPE,         //  Matrix has the wrong storage type (e.g.: ccs instead of crs)
    JMTX_RESULT_DIMS_MISMATCH,      //  Matrix or a vector don't match
    JMTX_RESULT_INDEX_OUT_OF_BOUNDS,//  Index was too large
    JMTX_RESULT_BAD_PARAM,          //  Parameter had a bad value
    JMTX_RESULT_NULL_PARAM,         //  Parameter was null

    JMTX_RESULT_BAD_MATRIX,         //  Matrix is fucked

    JMTX_RESULT_UNARY_RETURN,       //  Unary function returned non-zero

    JMTX_RESULT_COUNT,
};
typedef enum jmtx_result_enum jmtx_result;

const char* jmtx_result_to_str(jmtx_result res);



typedef uint32_t jmtx_index_t;

enum jmtx_matrix_type_enum
{
    JMTX_TYPE_NONE = 0, //  Invalid, here to force one of other values to be specified
    JMTX_TYPE_CRS,      //  Compressed row sparse
    JMTX_TYPE_CCS,      //  Compressed column sparse
    JMTX_TYPE_DRM,      //  Dense row-major
    JMTX_TYPE_DCM,      //  Dense column-major

    JMTX_TYPE_COUNT,    //  Here just as upper bound of what the enum should be
};
typedef enum jmtx_matrix_type_enum jmtx_matrix_type;

const char* jmtx_matrix_type_to_str(jmtx_matrix_type type);

typedef struct jmtx_allocator_callbacks_struct jmtx_allocator_callbacks;
struct jmtx_allocator_callbacks_struct
{
    void* (*alloc)(void* state, uint64_t size);
    void (*free)(void* state, void* ptr);
    void* (*realloc)(void* state, void* ptr, uint64_t new_size);
    void* state;
};

typedef struct jmtx_matrix_struct jmtx_matrix;
struct jmtx_matrix_struct
{
    jmtx_matrix_type type;
    uint32_t rows;
    uint32_t cols;
    jmtx_allocator_callbacks allocator_callbacks;
    jmtx_result (*get_element)(jmtx_matrix* mtx, jmtx_index_t row, jmtx_index_t col, float* p_out);
    int (*has_element)(jmtx_index_t row, jmtx_index_t col);
};

JMTX_GCC_ONLY(__attribute__((visibility("hidden"))))
uint32_t jmtx_internal_find_last_leq_value(uint32_t n_indices, const uint32_t* p_indices, uint32_t value);

extern const jmtx_allocator_callbacks JMTX_DEFAULT_ALLOCATOR_CALLBACKS;

#endif // !JMTX_MATRIX_BASE_H
