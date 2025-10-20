#ifndef JMTX_JMTX_FLOAT_DEFINES_H
#define JMTX_JMTX_FLOAT_DEFINES_H

#ifdef JMTX_TYPE_DEFINED
#    error Only one type can be defined at a time.
#else
#    define JMTX_TYPE_DEFINED "float"
#    define JMTX_TYPE_DEFINED_FLOAT
#endif

#define JMTX_REAL_T float
#define JMTX_SCALAR_T float

#include "jmtx_int_defines.h"
#define JMTX_NAME_TYPED(base_name) jmtxf_##base_name
#include <math.h>
#define JMTX_REAL_ROOT(x) sqrtf(x)
#define JMTX_FULL_ROOT(x) sqrtf(x)
#define JMTX_DOT(x, y) ((x) * (y))
#define JMTX_ABS(x) fabsf(x)

#endif // JMTX_JMTX_FLOAT_DEFINES_H
