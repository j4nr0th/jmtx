// Automatically generated from source/float/matrices/sparse_diagonal_compressed_internal.h on Fri Dec  1 17:36:03 2023
//
// Created by jan on 21.7.2023.
//

#ifndef JMTXC_SPARSE_DIAGONAL_COMPRESSED_INTERNAL_H
#define JMTXC_SPARSE_DIAGONAL_COMPRESSED_INTERNAL_H
#ifndef JMTXC_MATRIX_BASE_INTERNAL_H
    #include "../../matrix_base_internal.h"
#endif
#ifndef JMTXC_SPARSE_DIAGONAL_COMPRESSED_H
    #include "../../../include/jmtx/cfloat/matrices/sparse_diagonal_compressed.h"
#endif

struct jmtxc_matrix_cds_diagonal_array_T
{
    uint32_t capacity;  //  Number of diagonals
    uint32_t count;     //  Capacity for diagonals
    uint32_t* indices;  //  Indices of the diagonals relative to the main
    _Complex float** diagonals;  //  Array of the diagonal array pointers. Valid for count, space for up to capacity.
};
typedef struct jmtxc_matrix_cds_diagonal_array_T jmtxc_matrix_cds_diagonal_array;

struct jmtxc_matrix_cds_struct
{
    jmtx_matrix_base base;
    jmtxc_matrix_cds_diagonal_array super_diagonals;
    jmtxc_matrix_cds_diagonal_array sub_diagonals;
    _Complex float*  main_diagonal;
};



#endif //JMTXC_SPARSE_DIAGONAL_COMPRESSED_INTERNAL_H
