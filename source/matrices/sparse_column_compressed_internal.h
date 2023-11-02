//
// Created by jan on 21.7.2023.
//

#ifndef JMTX_SPARSE_COLUMN_COMPRESSED_INTERNAL_H
#define JMTX_SPARSE_COLUMN_COMPRESSED_INTERNAL_H
#include "sparse_column_compressed.h"
struct jmtx_matrix_ccs_struct
{
    jmtx_matrix base;
    //  How many elements exist in the columns left, so that column i is from index elements_before[i] ot elements_before[i + 1]
    uint32_t* elements_before;
    //  Column indices corresponding with the individual elements
    uint32_t* indices;
    //  Values of elements
    float* elements;
    uint32_t capacity;
    uint32_t n_elements;
};
#endif //JMTX_SPARSE_COLUMN_COMPRESSED_INTERNAL_H
