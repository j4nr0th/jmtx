#ifndef MTXLIB_MATRIX_BASE_H
#define MTXLIB_MATRIX_BASE_H

#include <stdint.h>
#include <stdlib.h>
#include <string.h>


typedef enum mtx_res_t_enum
{
    mtx_success         = 0,    //  All is good
    mtx_fail            = -1,   //  Generic error code
    mtx_not_converged   = -2,   //  Iteration didn't converge
    mtx_malloc_fail     = -3,   //  malloc/calloc/realloc being a cunt
    mtx_thread_fail     = -4,   //  Couldn't create a thread
#ifdef MTX_MATRIX_CHECKS        //  Checks for more detailed errors
    mtx_wrong_storage   = -5,   //  Matrix has the wrong storage type (e.g.: ccs instead of crs)
    mtx_wrong_size      = -6,   //  Matrix or a vector don't match
    mtx_out_of_range    = -7,   //  Index was too large
    mtx_bad_param       = -8,   //  Parameter had a bad value
#endif
} mtx_res_t;



typedef float scalar_t;

typedef enum enum_Matrix_Type
{
    mtx_type_crs, //  Compressed row sparse
    mtx_type_ccs, //  Compressed column sparse
} MatrixType;

#define MATRIX_STRUCT_BASE uint type;\
    uint rows;\
    uint columns;\
    uint n_elements;

#endif // !MTXLIB_MATRIX_BASE_H
