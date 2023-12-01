#ifndef JMTX_MATRIX_BASE_INTERNAL_H
#define JMTX_MATRIX_BASE_INTERNAL_H
#define JMTXD_MATRIX_BASE_INTERNAL_H

#ifndef JMTX_MATRIX_BASE_H
    #include "../include/jmtx/matrix_base.h"
#endif

enum jmtx_matrix_type_T
{
    JMTX_TYPE_NONE = 0, //  Invalid, here to force one of other values to be specified

    //  Single precision types
    JMTX_TYPE_CRS,      //  Compressed row sparse
    JMTX_TYPE_CCS,      //  Compressed column sparse
    JMTX_TYPE_BRM,      //  Band row-major
    JMTX_TYPE_CDS,      //  Compressed diagonal sparse

    //  Double precision
    JMTXD_TYPE_CRS,      //  Compressed row sparse
    JMTXD_TYPE_CCS,      //  Compressed column sparse
    JMTXD_TYPE_BRM,      //  Band row-major
    JMTXD_TYPE_CDS,      //  Compressed diagonal sparse


    JMTX_TYPE_COUNT,    //  Here just as upper bound of what the enum should be
};
typedef enum jmtx_matrix_type_T jmtx_matrix_type;


const char* jmtx_matrix_type_to_str(jmtx_matrix_type type);

typedef struct jmtx_matrix_base_T jmtx_matrix_base;
struct jmtx_matrix_base_T
{
    jmtx_matrix_type type;
    uint32_t rows;
    uint32_t cols;
    jmtx_allocator_callbacks allocator_callbacks;
};

JMTX_INTERNAL_FUNCTION
uint32_t jmtx_internal_find_last_leq_value(uint32_t n_indices, const uint32_t p_indices[static n_indices], uint32_t value);

extern const jmtx_allocator_callbacks JMTX_DEFAULT_ALLOCATOR_CALLBACKS;

#endif // !JMTX_MATRIX_BASE_INTERNAL_H
