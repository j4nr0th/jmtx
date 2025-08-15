// Automatically generated from source/float/matrices/band_row_major_internal.h on Fri Dec  1 17:36:03 2023
//
// Created by jan on 21.7.2023.
//

#ifndef JMTXC_BAND_ROW_MAJOR_INTERNAL_H
#define JMTXC_BAND_ROW_MAJOR_INTERNAL_H
#ifndef JMTX_MATRIX_BASE_INTERNAL_H
#include "../../matrix_base_internal.h"
#endif
#ifndef JMTXC_BAND_ROW_MAJOR_H
#include "../../../include/jmtx/cfloat/matrices/band_row_major.h"
#endif
struct jmtxc_matrix_brm_struct
{
    jmtx_matrix_base base;
    uint32_t upper_bandwidth;
    uint32_t lower_bandwidth;
    _Complex float *restrict values;
};
#endif // JMTXC_BAND_ROW_MAJOR_INTERNAL_H
