//
// Created by jan on 21.7.2023.
//

#ifndef JMTX_SPARSE_ROW_COMPRESSED_INTERNAL_H
#define JMTX_SPARSE_ROW_COMPRESSED_INTERNAL_H
#include "../../include/jmtx/matrices/sparse_row_compressed.h"
struct jmtx_matrix_crs_struct
{
    jmtx_matrix base;
    //  end_of_row_offsets[i] has the number of values are there before the end of the row i
    uint32_t* end_of_row_offsets;
    //  Column indices corresponding with the individual values
    uint32_t* indices;
    //  Values of values
    float* values;
    uint32_t n_entries;
    uint32_t capacity;
};
#endif //JMTX_SPARSE_ROW_COMPRESSED_INTERNAL_H
