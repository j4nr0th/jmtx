// Automatically generated from source/cfloat/matrices/sparse_diagonal_compressed_internal.h on Fri Dec  1 18:48:06 2023
// Automatically generated from source/cdouble/matrices/sparse_diagonal_compressed_internal.h on Fri Dec  1 17:36:03 2023
//
// Created by jan on 21.7.2023.
//

#ifndef JMTXZ_SPARSE_DIAGONAL_COMPRESSED_INTERNAL_H
#define JMTXZ_SPARSE_DIAGONAL_COMPRESSED_INTERNAL_H
#ifndef JMTX_MATRIX_BASE_INTERNAL_H
    #include "../../matrix_base_internal.h"
#endif
#ifndef JMTXZ_SPARSE_DIAGONAL_COMPRESSED_H
    #include "../../../include/jmtx/cdouble/matrices/sparse_diagonal_compressed.h"
#endif

struct jmtxz_matrix_cds_diagonal_array_T
{
    uint32_t capacity;  //  Number of diagonals
    uint32_t count;     //  Capacity for diagonals
    uint32_t* indices;  //  Indices of the diagonals relative to the main
    _Complex double** diagonals;  //  Array of the diagonal array pointers. Valid for count, space for up to capacity.
};
typedef struct jmtxz_matrix_cds_diagonal_array_T jmtxz_matrix_cds_diagonal_array;

struct jmtxz_matrix_cds_struct
{
    jmtx_matrix_base base;
    jmtxz_matrix_cds_diagonal_array super_diagonals;
    jmtxz_matrix_cds_diagonal_array sub_diagonals;
    _Complex double* restrict main_diagonal;
};



#endif //JMTXZ_SPARSE_DIAGONAL_COMPRESSED_INTERNAL_H
