#ifndef JMTX_JMTX_CONVERT_H
#define JMTX_JMTX_CONVERT_H

#if defined(JMTXF_H) && defined(JMTXD_H)
// Conversion functions for float and double types

// double -> float

#    define JMTX_INPUT_TYPE double
#    define JMTX_OUTPUT_TYPE float
#    define JMTX_NAME_IN(x) jmtxd_##x
#    define JMTX_NAME_OUT(x) jmtxf_##x
#    define JMTX_NAME_CONVERSION(x) jmtx_convert_##x##_d2f

#    include "jmtx_int_defines.h"
#    include "../conversion/brm_conversion.h"
#    include "../conversion/ccs_conversion.h"
#    include "../conversion/cds_conversion.h"
#    include "../conversion/crs_conversion.h"

#    undef JMTX_INPUT_TYPE
#    undef JMTX_OUTPUT_TYPE
#    undef JMTX_NAME_IN
#    undef JMTX_NAME_OUT
#    undef JMTX_NAME_CONVERSION
#    undef JMTX_BRM_CONVERSION_H

// float -> double

#    define JMTX_INPUT_TYPE float
#    define JMTX_OUTPUT_TYPE double
#    define JMTX_NAME_IN(x) jmtxf_##x
#    define JMTX_NAME_OUT(x) jmtxd_##x
#    define JMTX_NAME_CONVERSION(x) jmtx_convert_##x##_f2d

#    include "jmtx_int_defines.h"
#    include "../conversion/brm_conversion.h"
#    include "../conversion/ccs_conversion.h"
#    include "../conversion/cds_conversion.h"
#    include "../conversion/crs_conversion.h"

#    undef JMTX_INPUT_TYPE
#    undef JMTX_OUTPUT_TYPE
#    undef JMTX_NAME_IN
#    undef JMTX_NAME_OUT
#    undef JMTX_NAME_CONVERSION
#    undef JMTX_BRM_CONVERSION_H

#endif // JMTXF_H && JMTXD_H

#if defined(JMTXD_H) && defined(JMTXZ_H)

// cdouble -> double

#    define JMTX_INPUT_TYPE _Complex double
#    define JMTX_OUTPUT_TYPE double
#    define JMTX_NAME_IN(x) jmtxz_##x
#    define JMTX_NAME_OUT(x) jmtxd_##x
#    define JMTX_NAME_CONVERSION(x) jmtx_convert_##x##_z2d

#    include "jmtx_int_defines.h"
#    include "../conversion/brm_conversion.h"
#    include "../conversion/ccs_conversion.h"
#    include "../conversion/cds_conversion.h"
#    include "../conversion/crs_conversion.h"

#    undef JMTX_INPUT_TYPE
#    undef JMTX_OUTPUT_TYPE
#    undef JMTX_NAME_IN
#    undef JMTX_NAME_OUT
#    undef JMTX_NAME_CONVERSION
#    undef JMTX_BRM_CONVERSION_H

// double -> cdouble

#    define JMTX_INPUT_TYPE double
#    define JMTX_OUTPUT_TYPE _Complex double
#    define JMTX_NAME_IN(x) jmtxd_##x
#    define JMTX_NAME_OUT(x) jmtxz_##x
#    define JMTX_NAME_CONVERSION(x) jmtx_convert_##x##_d2z

#    include "jmtx_int_defines.h"
#    include "../conversion/brm_conversion.h"
#    include "../conversion/ccs_conversion.h"
#    include "../conversion/cds_conversion.h"
#    include "../conversion/crs_conversion.h"

#    undef JMTX_INPUT_TYPE
#    undef JMTX_OUTPUT_TYPE
#    undef JMTX_NAME_IN
#    undef JMTX_NAME_OUT
#    undef JMTX_NAME_CONVERSION
#    undef JMTX_BRM_CONVERSION_H

#endif // JMTXD_H && JMTXZ_H

#if defined(JMTXZ_H) && defined(JMTXC_H)

// cdouble -> cfloat

#    define JMTX_INPUT_TYPE _Complex double
#    define JMTX_OUTPUT_TYPE _Complex float
#    define JMTX_NAME_IN(x) jmtxz_##x
#    define JMTX_NAME_OUT(x) jmtxc_##x
#    define JMTX_NAME_CONVERSION(x) jmtx_convert_##x##_z2c

#    include "jmtx_int_defines.h"
#    include "../conversion/brm_conversion.h"
#    include "../conversion/ccs_conversion.h"
#    include "../conversion/cds_conversion.h"
#    include "../conversion/crs_conversion.h"

#    undef JMTX_INPUT_TYPE
#    undef JMTX_OUTPUT_TYPE
#    undef JMTX_NAME_IN
#    undef JMTX_NAME_OUT
#    undef JMTX_NAME_CONVERSION
#    undef JMTX_BRM_CONVERSION_H

// cfloat -> cdouble

#    define JMTX_INPUT_TYPE _Complex float
#    define JMTX_OUTPUT_TYPE _Complex double
#    define JMTX_NAME_IN(x) jmtxc_##x
#    define JMTX_NAME_OUT(x) jmtxz_##x
#    define JMTX_NAME_CONVERSION(x) jmtx_convert_##x##_c2z

#    include "jmtx_int_defines.h"
#    include "../conversion/brm_conversion.h"
#    include "../conversion/ccs_conversion.h"
#    include "../conversion/cds_conversion.h"
#    include "../conversion/crs_conversion.h"

#    undef JMTX_INPUT_TYPE
#    undef JMTX_OUTPUT_TYPE
#    undef JMTX_NAME_IN
#    undef JMTX_NAME_OUT
#    undef JMTX_NAME_CONVERSION
#    undef JMTX_BRM_CONVERSION_H

#endif // JMTXZ_H && JMTXC_H

#if defined(JMTXC_H) && defined(JMTXF_H)

// cfloat -> float

#    define JMTX_INPUT_TYPE _Complex float
#    define JMTX_OUTPUT_TYPE float
#    define JMTX_NAME_IN(x) jmtxc_##x
#    define JMTX_NAME_OUT(x) jmtxf_##x
#    define JMTX_NAME_CONVERSION(x) jmtx_convert_##x##_c2f

#    include "jmtx_int_defines.h"
#    include "../conversion/brm_conversion.h"
#    include "../conversion/ccs_conversion.h"
#    include "../conversion/cds_conversion.h"
#    include "../conversion/crs_conversion.h"

#    undef JMTX_INPUT_TYPE
#    undef JMTX_OUTPUT_TYPE
#    undef JMTX_NAME_IN
#    undef JMTX_NAME_OUT
#    undef JMTX_NAME_CONVERSION
#    undef JMTX_BRM_CONVERSION_H

// float -> cfloat

#    define JMTX_INPUT_TYPE float
#    define JMTX_OUTPUT_TYPE _Complex float
#    define JMTX_NAME_IN(x) jmtxf_##x
#    define JMTX_NAME_OUT(x) jmtxc_##x
#    define JMTX_NAME_CONVERSION(x) jmtx_convert_##x##_f2c

#    include "../conversion/brm_conversion.h"
#    include "../conversion/ccs_conversion.h"
#    include "../conversion/cds_conversion.h"
#    include "../conversion/crs_conversion.h"

#    undef JMTX_INPUT_TYPE
#    undef JMTX_OUTPUT_TYPE
#    undef JMTX_NAME_IN
#    undef JMTX_NAME_OUT
#    undef JMTX_NAME_CONVERSION
#    undef JMTX_BRM_CONVERSION_H

#endif // JMTXC_H && JMTXF_H

#if defined(JMTXF_H) && defined(JMTXZ_H)

// float -> cdouble

#    define JMTX_INPUT_TYPE float
#    define JMTX_OUTPUT_TYPE _Complex double
#    define JMTX_NAME_IN(x) jmtxf_##x
#    define JMTX_NAME_OUT(x) jmtxz_##x
#    define JMTX_NAME_CONVERSION(x) jmtx_convert_##x##_f2z

#    include "jmtx_int_defines.h"
#    include "../conversion/brm_conversion.h"
#    include "../conversion/ccs_conversion.h"
#    include "../conversion/cds_conversion.h"
#    include "../conversion/crs_conversion.h"

#    undef JMTX_INPUT_TYPE
#    undef JMTX_OUTPUT_TYPE
#    undef JMTX_NAME_IN
#    undef JMTX_NAME_OUT
#    undef JMTX_NAME_CONVERSION
#    undef JMTX_BRM_CONVERSION_H

// cdouble -> float

#    define JMTX_INPUT_TYPE _Complex double
#    define JMTX_OUTPUT_TYPE float
#    define JMTX_NAME_IN(x) jmtxz_##x
#    define JMTX_NAME_OUT(x) jmtxf_##x
#    define JMTX_NAME_CONVERSION(x) jmtx_convert_##x##_z2f

#    include "jmtx_int_defines.h"
#    include "../conversion/brm_conversion.h"
#    include "../conversion/ccs_conversion.h"
#    include "../conversion/cds_conversion.h"
#    include "../conversion/crs_conversion.h"

#    undef JMTX_INPUT_TYPE
#    undef JMTX_OUTPUT_TYPE
#    undef JMTX_NAME_IN
#    undef JMTX_NAME_OUT
#    undef JMTX_NAME_CONVERSION
#    undef JMTX_BRM_CONVERSION_H

#endif // JMTXF_H && JMTXZ_H

#if defined(JMTXD_H) && defined(JMTXC_H)

// float -> cdouble

#    define JMTX_INPUT_TYPE float
#    define JMTX_OUTPUT_TYPE _Complex double
#    define JMTX_NAME_IN(x) jmtxf_##x
#    define JMTX_NAME_OUT(x) jmtxz_##x
#    define JMTX_NAME_CONVERSION(x) jmtx_convert_##x##_f2z

#    include "jmtx_int_defines.h"
#    include "../conversion/brm_conversion.h"
#    include "../conversion/ccs_conversion.h"
#    include "../conversion/cds_conversion.h"
#    include "../conversion/crs_conversion.h"

#    undef JMTX_INPUT_TYPE
#    undef JMTX_OUTPUT_TYPE
#    undef JMTX_NAME_IN
#    undef JMTX_NAME_OUT
#    undef JMTX_NAME_CONVERSION
#    undef JMTX_BRM_CONVERSION_H

// cfloat -> double

#    define JMTX_INPUT_TYPE _Complex float
#    define JMTX_OUTPUT_TYPE double
#    define JMTX_NAME_IN(x) jmtxc_##x
#    define JMTX_NAME_OUT(x) jmtxd_##x
#    define JMTX_NAME_CONVERSION(x) jmtx_convert_##x##_c2d

#    include "jmtx_int_defines.h"
#    include "../conversion/brm_conversion.h"
#    include "../conversion/ccs_conversion.h"
#    include "../conversion/cds_conversion.h"
#    include "../conversion/crs_conversion.h"

#    undef JMTX_INPUT_TYPE
#    undef JMTX_OUTPUT_TYPE
#    undef JMTX_NAME_IN
#    undef JMTX_NAME_OUT
#    undef JMTX_NAME_CONVERSION
#    undef JMTX_BRM_CONVERSION_H

#endif // JMTXD_H && JMTXC_H

#endif // JMTX_JMTX_CONVERT_H
