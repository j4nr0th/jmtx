#ifndef JMTX_JMTX_DOUBLE_DEFINES_H
#define JMTX_JMTX_DOUBLE_DEFINES_H

#ifdef JMTX_TYPE_DEFINED
#    error Only one type can be defined at a time.
#else
#    define JMTX_TYPE_DEFINED "double"
#    define JMTX_TYPE_DEFINED_DOUBLE
#endif

#define JMTX_REAL_T double
#define JMTX_SCALAR_T double
#include "jmtx_int_defines.h"

#define JMTX_NAME_TYPED(base_name) jmtxd_##base_name
#include <math.h>
#define JMTX_REAL_ROOT(x) sqrt(x)
#define JMTX_FULL_ROOT(x) sqrt(x)
#define JMTX_DOT(x, y) ((x) * (y))
#define JMTX_ABS(x) fabs(x)

#endif // JMTX_JMTX_DOUBLE_DEFINES_H
