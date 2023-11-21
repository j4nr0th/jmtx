//
// Created by jan on 21.7.2023.
//

#ifndef JMTX_SPARSE_COLUMN_COMPRESSED_INTERNAL_H
#define JMTX_SPARSE_COLUMN_COMPRESSED_INTERNAL_H
#include "../../include/jmtx/matrices/sparse_column_compressed.h"
struct jmtx_matrix_ccs_struct
{
    jmtx_matrix base;
    //  How many values exist in the columns left, so that column i is from index end_of_column_offsets[i] ot end_of_column_offsets[i + 1]
    uint32_t* end_of_column_offsets;
    //  Column indices corresponding with the individual values
    uint32_t* indices;
    //  Values of values
    float* values;
    uint32_t capacity;
    uint32_t n_entries;
};
#endif //JMTX_SPARSE_COLUMN_COMPRESSED_INTERNAL_H
