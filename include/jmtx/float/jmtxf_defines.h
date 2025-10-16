//
// Created by jan on 2025-08-26.
//

#ifndef JMTXF_JMTX_DEFINES_H
#define JMTXF_JMTX_DEFINES_H

#define JMTXF_SOURCE_ROOT_PATH "source/float"
#define JMTXF_HEADER_ROOT_PATH "source/float"
#include <math.h>

typedef float jmtxf_scalar;
typedef float jmtxf_real;

inline jmtxf_real jmtxf_scalar_abs(jmtxf_scalar x)
{
    return fabsf(x);
}

inline jmtxf_real jmtxf_scalar_mag2(jmtxf_scalar x)
{
    return x * x;
}

inline jmtxf_scalar jmtxf_scalar_sqrt(jmtxf_scalar x)
{
    return sqrtf(x);
}

#endif // JMTXF_JMTX_DEFINES_H
