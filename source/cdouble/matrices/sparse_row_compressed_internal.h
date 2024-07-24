// Automatically generated from source/cfloat/matrices/sparse_row_compressed_internal.h on Fri Dec  1 18:48:06 2023
// Automatically generated from source/cdouble/matrices/sparse_row_compressed_internal.h on Fri Dec  1 17:36:03 2023
//
// Created by jan on 21.7.2023.
//

#ifndef JMTXZ_SPARSE_ROW_COMPRESSED_INTERNAL_H
#define JMTXZ_SPARSE_ROW_COMPRESSED_INTERNAL_H
#ifndef JMTX_MATRIX_BASE_INTERNAL_H
#include "../../matrix_base_internal.h"
#endif
#ifndef JMTXZ_SPARSE_ROW_COMPRESSED_H
    #include "../../../include/jmtx/cdouble/matrices/sparse_row_compressed.h"
#endif

struct jmtxz_matrix_crs_struct
{
    jmtx_matrix_base base;
    //  end_of_row_offsets[i] has the number of values are there before the end of the row i
    uint32_t* restrict end_of_row_offsets;
    //  Column indices corresponding with the individual values
    uint32_t* restrict indices;
    //  Values of values
    _Complex double* restrict values;
    uint32_t n_entries;
    uint32_t capacity;
};
#endif //JMTXZ_SPARSE_ROW_COMPRESSED_INTERNAL_H
