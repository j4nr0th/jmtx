// Automatically generated from source/cfloat/matrices/sparse_column_compressed_internal.h on Fri Dec  1 18:48:06 2023
// Automatically generated from source/cdouble/matrices/sparse_column_compressed_internal.h on Fri Dec  1 17:36:03 2023
//
// Created by jan on 21.7.2023.
//

#ifndef JMTXZ_SPARSE_COLUMN_COMPRESSED_INTERNAL_H
#define JMTXZ_SPARSE_COLUMN_COMPRESSED_INTERNAL_H
#ifndef JMTXZ_MATRIX_BASE_INTERNAL_H
    #include "../../matrix_base_internal.h"
#endif

#ifndef JMTXZ_SPARSE_COLUMN_COMPRESSED_H
    #include "../../../include/jmtx/cdouble/matrices/sparse_column_compressed.h"
#endif
struct jmtxz_matrix_ccs_struct
{
    jmtx_matrix_base base;
    //  How many values exist in the columns left, so that column i is from index end_of_column_offsets[i] ot end_of_column_offsets[i + 1]
    uint32_t* end_of_column_offsets;
    //  Column indices corresponding with the individual values
    uint32_t* indices;
    //  Values of values
    _Complex double* values;
    uint32_t n_entries;
    uint32_t capacity;
};
#endif //JMTXZ_SPARSE_COLUMN_COMPRESSED_INTERNAL_H
