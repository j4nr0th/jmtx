#ifndef JMTX_JMTX_CFLOAT_DEFINES_H
#define JMTX_JMTX_CFLOAT_DEFINES_H

#ifdef JMTX_TYPE_DEFINED
#    error Only one type can be defined at a time.
#else
#    define JMTX_TYPE_DEFINED "_Complex float"
#    define JMTX_TYPE_DEFINED_CFLOAT
#endif

#define JMTX_REAL_T float
#define JMTX_SCALAR_T _Complex float
#define JMTX_INDEX_T uint32_t
#define JMTX_FAST_INT_T uint_fast32_t
#define JMTX_SIZE_T size_t

#define JMTX_NAME_TYPED(base_name) jmtxc_##base_name

#include <complex.h>
#define JMTX_REAL_ROOT(x) csqrtf(x)
#define JMTX_FULL_ROOT(x) csqrtf(x)
#define JMTX_DOT(x, y) (conjf(x) * (y))
#define JMTX_ABS(x) cabsf(x)

#endif // JMTX_JMTX_CFLOAT_DEFINES_H
