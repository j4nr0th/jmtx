//
// Created by jan on 21.7.2023.
//

#ifndef JMTX_DENSE_ROW_MAJOR_INTERNAL_H
#define JMTX_DENSE_ROW_MAJOR_INTERNAL_H
#ifndef JMTX_MATRIX_BASE_INTERNAL_H
    #include "../../matrix_base_internal.h"
#endif
#ifndef JMTX_DENSE_ROW_MAJOR_H
    #include "../../../include/jmtx/float/matrices/dense_row_major.h"
#endif
struct jmtx_matrix_drm_struct
{
    jmtx_matrix_base base;
    //  Length this->base.rows, contains the order in which the rows are permuted (what row is read from where)
    uint32_t* permutations;
    //  Length this->base.rows, contains the reverse order in which the rows are permuted (what row is written to where)
    uint32_t* rperm;
    //  Length this->base.rows * this->base.cols, contains every entry in row-major ordering.
    float* restrict values;
};
#endif //JMTX_DENSE_ROW_MAJOR_INTERNAL_H
