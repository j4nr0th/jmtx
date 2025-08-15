//
// Created by jan on 21.7.2023.
//

#ifndef JMTX_SPARSE_DIAGONAL_COMPRESSED_INTERNAL_H
#define JMTX_SPARSE_DIAGONAL_COMPRESSED_INTERNAL_H
#ifndef JMTX_MATRIX_BASE_INTERNAL_H
#    include "../../matrix_base_internal.h"
#endif
#ifndef JMTX_SPARSE_DIAGONAL_COMPRESSED_H
#    include "../../../include/jmtx/float/matrices/sparse_diagonal_compressed.h"
#endif

struct jmtx_matrix_cds_diagonal_array_T
{
    uint32_t capacity; //  Number of diagonals
    uint32_t count;    //  Capacity for diagonals
    uint32_t *indices; //  Indices of the diagonals relative to the main
    float **diagonals; //  Array of the diagonal array pointers. Valid for count, space for up to capacity.
};
typedef struct jmtx_matrix_cds_diagonal_array_T jmtx_matrix_cds_diagonal_array;

struct jmtx_matrix_cds_struct
{
    jmtx_matrix_base base;
    jmtx_matrix_cds_diagonal_array super_diagonals;
    jmtx_matrix_cds_diagonal_array sub_diagonals;
    float *restrict main_diagonal;
};

#endif // JMTX_SPARSE_DIAGONAL_COMPRESSED_INTERNAL_H
