#ifndef JMTX_JMTX_CDOUBLE_DEFINES_H
#define JMTX_JMTX_CDOUBLE_DEFINES_H

#ifdef JMTX_TYPE_DEFINED
#    error Only one type can be defined at a time.
#else
#    define JMTX_TYPE_DEFINED "_Complex double"
#    define JMTX_TYPE_DEFINED_CDOUBLE
#endif

#define JMTX_REAL_T double
#define JMTX_SCALAR_T _Complex double
#include "jmtx_int_defines.h"

#define JMTX_NAME_TYPED(base_name) jmtxz_##base_name

#include <complex.h>
#define JMTX_REAL_ROOT(x) csqrt(x)
#define JMTX_FULL_ROOT(x) csqrt(x)
#define JMTX_DOT(x, y) ((JMTX_REAL_T)(conj(x) * (y)))
#define JMTX_ABS(x) cabs(x)

#endif // JMTX_JMTX_CDOUBLE_DEFINES_H
