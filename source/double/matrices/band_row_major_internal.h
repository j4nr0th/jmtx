// Automatically generated from source/float/matrices/band_row_major_internal.h on Thu Nov 30 20:33:08 2023
//
// Created by jan on 21.7.2023.
//

#ifndef JMTX_BAND_ROW_MAJOR_INTERNAL_H
#define JMTX_BAND_ROW_MAJOR_INTERNAL_H
#ifndef JMTX_MATRIX_BASE_INTERNAL_H
    #include "../../matrix_base_internal.h"
#endif
#ifndef JMTX_BAND_ROW_MAJOR_H
    #include "../../../include/jmtx/double/matrices/band_row_major.h"
#endif
struct jmtxd_matrix_brm_struct
{
    jmtx_matrix_base base;
    uint32_t upper_bandwidth;
    uint32_t lower_bandwidth;
    double* restrict values;
};
#endif //JMTX_BAND_ROW_MAJOR_INTERNAL_H
