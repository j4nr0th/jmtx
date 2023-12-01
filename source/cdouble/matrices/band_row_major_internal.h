// Automatically generated from source/cfloat/matrices/band_row_major_internal.h on Fri Dec  1 18:48:06 2023
// Automatically generated from source/cdouble/matrices/band_row_major_internal.h on Fri Dec  1 17:36:03 2023
//
// Created by jan on 21.7.2023.
//

#ifndef JMTXZ_BAND_ROW_MAJOR_INTERNAL_H
#define JMTXZ_BAND_ROW_MAJOR_INTERNAL_H
#ifndef JMTXZ_MATRIX_BASE_INTERNAL_H
    #include "../../matrix_base_internal.h"
#endif
#ifndef JMTXZ_BAND_ROW_MAJOR_H
    #include "../../../include/jmtx/cdouble/matrices/band_row_major.h"
#endif
struct jmtxz_matrix_brm_struct
{
    jmtx_matrix_base base;
    uint32_t upper_bandwidth;
    uint32_t lower_bandwidth;
    _Complex double* restrict values;
};
#endif //JMTXZ_BAND_ROW_MAJOR_INTERNAL_H
