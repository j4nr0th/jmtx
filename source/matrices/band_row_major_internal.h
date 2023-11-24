//
// Created by jan on 21.7.2023.
//

#ifndef JMTX_BAND_ROW_MAJOR_INTERNAL_H
#define JMTX_BAND_ROW_MAJOR_INTERNAL_H
#include "matrix_base_internal.h"
#include "../../include/jmtx/matrices/band_row_major.h"
struct jmtx_matrix_brm_struct
{
    jmtx_matrix_base base;
    uint32_t upper_bandwidth;
    uint32_t lower_bandwidth;
    float* restrict values;
};
#endif //JMTX_BAND_ROW_MAJOR_INTERNAL_H
