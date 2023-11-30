// Automatically generated from source/float/matrices/sparse_row_compressed_internal.h on Thu Nov 30 20:33:08 2023
//
// Created by jan on 21.7.2023.
//

#ifndef JMTX_SPARSE_ROW_COMPRESSED_INTERNAL_H
#define JMTX_SPARSE_ROW_COMPRESSED_INTERNAL_H
#ifndef JMTX_MATRIX_BASE_INTERNAL_H
#include "../../matrix_base_internal.h"
#endif
#ifndef JMTX_SPARSE_ROW_COMPRESSED_H
    #include "../../../include/jmtx/double/matrices/sparse_row_compressed.h"
#endif

struct jmtxd_matrix_crs_struct
{
    jmtx_matrix_base base;
    //  end_of_row_offsets[i] has the number of values are there before the end of the row i
    uint32_t* end_of_row_offsets;
    //  Column indices corresponding with the individual values
    uint32_t* indices;
    //  Values of values
    double* values;
    uint32_t n_entries;
    uint32_t capacity;
};
#endif //JMTX_SPARSE_ROW_COMPRESSED_INTERNAL_H
