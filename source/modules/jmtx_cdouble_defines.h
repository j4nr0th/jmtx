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
#define JMTX_INDEX_T uint32_t
#define JMTX_FAST_INT_T uint_fast32_t
#define JMTX_SIZE_T size_t

#define JMTX_NAME_TYPED(base_name) jmtxz_##base_name

#include <complex.h>
#define JMTX_REAL_ROOT(x) csqrt(x)
#define JMTX_FULL_ROOT(x) csqrt(x)
#define JMTX_DOT(x, y) (conj(x) * (y))
#define JMTX_ABS(x) cabs(x)

#endif // JMTX_JMTX_CDOUBLE_DEFINES_H
