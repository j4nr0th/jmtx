//
// Created by jan on 21.7.2023.
//

#ifndef JMTX_SPARSE_COLUMN_COMPRESSED_INTERNAL_H
#define JMTX_SPARSE_COLUMN_COMPRESSED_INTERNAL_H
#ifndef JMTX_MATRIX_BASE_INTERNAL_H
    #include "../../matrix_base_internal.h"
#endif

#ifndef JMTX_SPARSE_COLUMN_COMPRESSED_H
    #include "../../../include/jmtx/float/matrices/sparse_column_compressed.h"
#endif
struct jmtx_matrix_ccs_struct
{
    jmtx_matrix_base base;
    //  How many values exist in the columns left, so that column i is from index end_of_column_offsets[i] ot end_of_column_offsets[i + 1]
    uint32_t* end_of_column_offsets;
    //  Column indices corresponding with the individual values
    uint32_t* indices;
    //  Values of values
    float* values;
    uint32_t n_entries;
    uint32_t capacity;
};
#endif //JMTX_SPARSE_COLUMN_COMPRESSED_INTERNAL_H
