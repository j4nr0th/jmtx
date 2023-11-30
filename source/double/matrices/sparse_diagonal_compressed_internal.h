// Automatically generated from source/float/matrices/sparse_diagonal_compressed_internal.h on Thu Nov 30 20:33:08 2023
//
// Created by jan on 21.7.2023.
//

#ifndef JMTX_SPARSE_DIAGONAL_COMPRESSED_INTERNAL_H
#define JMTX_SPARSE_DIAGONAL_COMPRESSED_INTERNAL_H
#ifndef JMTX_MATRIX_BASE_INTERNAL_H
    #include "../../matrix_base_internal.h"
#endif
#ifndef JMTX_JMTX_SPARSE_DIAGONAL_COMPRESSED_H
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



#endif //JMTX_SPARSE_DIAGONAL_COMPRESSED_INTERNAL_H
