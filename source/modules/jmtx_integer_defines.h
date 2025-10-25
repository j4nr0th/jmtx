#ifndef JMTX_JMTX_INT_DEFINES_H
#define JMTX_JMTX_INT_DEFINES_H

#ifdef JMTX_TYPE_DEFINED
#    error Only one type can be defined at a time.
#else
#    define JMTX_TYPE_DEFINED "int32"
#    define JMTX_TYPE_DEFINED_INT
#endif

#include <stdint.h>
#define JMTX_REAL_T int32_t
#define JMTX_SCALAR_T int32_t

#include "jmtx_int_defines.h"
#define JMTX_NAME_TYPED(base_name) jmtxi_##base_name
#include <math.h>
#define JMTX_REAL_ROOT(x) (JMTX_REAL_T) sqrt((double)x)
#define JMTX_FULL_ROOT(x) (JMTX_SCALAR_T) sqrt((double)x)
#define JMTX_DOT(x, y) ((x) * (y))
#define JMTX_ABS(x) (JMTX_SCALAR_T) abs(x)

#endif // JMTX_JMTX_INT_DEFINES_H
