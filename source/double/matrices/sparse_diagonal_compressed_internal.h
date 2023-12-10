// Automatically generated from source/float/matrices/sparse_diagonal_compressed_internal.h on Fri Dec  1 06:43:01 2023
//
// Created by jan on 21.7.2023.
//

#ifndef JMTXD_SPARSE_DIAGONAL_COMPRESSED_INTERNAL_H
#define JMTXD_SPARSE_DIAGONAL_COMPRESSED_INTERNAL_H
#ifndef JMTX_MATRIX_BASE_INTERNAL_H
    #include "../../matrix_base_internal.h"
#endif
#ifndef JMTXD_SPARSE_DIAGONAL_COMPRESSED_H
    #include "../../../include/jmtx/double/matrices/sparse_diagonal_compressed.h"
#endif

struct jmtxd_matrix_cds_diagonal_array_T
{
    uint32_t capacity;  //  Number of diagonals
    uint32_t count;     //  Capacity for diagonals
    uint32_t* indices;  //  Indices of the diagonals relative to the main
    double** diagonals;  //  Array of the diagonal array pointers. Valid for count, space for up to capacity.
};
typedef struct jmtxd_matrix_cds_diagonal_array_T jmtxd_matrix_cds_diagonal_array;

struct jmtxd_matrix_cds_struct
{
    jmtx_matrix_base base;
    jmtxd_matrix_cds_diagonal_array super_diagonals;
    jmtxd_matrix_cds_diagonal_array sub_diagonals;
    double*  main_diagonal;
};



#endif //JMTXD_SPARSE_DIAGONAL_COMPRESSED_INTERNAL_H
