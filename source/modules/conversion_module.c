
// double -> float

#define JMTX_INPUT_TYPE double
#define JMTX_OUTPUT_TYPE float
#define JMTX_NAME_IN(x) jmtxd_##x
#define JMTX_NAME_OUT(x) jmtxf_##x
#define JMTX_OUTPUT_ENUM(t) JMTXF_TYPE_##t
#define JMTX_NAME_CONVERSION(x) jmtx_convert_##x##_d2f

#include "../../include/jmtx/jmtxd.h"
#include "../../include/jmtx/jmtxf.h"
#include "jmtx_int_defines.h"
#include "../conversion/brm_conversion.c"
#include "../conversion/ccs_conversion.c"
#include "../conversion/cds_conversion.c"
#include "../conversion/crs_conversion.c"

#undef JMTX_INPUT_TYPE
#undef JMTX_OUTPUT_TYPE
#undef JMTX_NAME_IN
#undef JMTX_NAME_OUT
#undef JMTX_OUTPUT_ENUM
#undef JMTX_NAME_CONVERSION
#undef JMTX_BRM_CONVERSION_H

// float -> cfloat

#define JMTX_INPUT_TYPE float
#define JMTX_OUTPUT_TYPE _Complex float
#define JMTX_NAME_IN(x) jmtxf_##x
#define JMTX_NAME_OUT(x) jmtxc_##x
#define JMTX_OUTPUT_ENUM(t) JMTXC_TYPE_##t
#define JMTX_NAME_CONVERSION(x) jmtx_convert_##x##_f2c

#include "../../include/jmtx/jmtxf.h"
#include "../../include/jmtx/jmtxc.h"
#include "jmtx_int_defines.h"
#include "../conversion/brm_conversion.c"
#include "../conversion/ccs_conversion.c"
#include "../conversion/cds_conversion.c"
#include "../conversion/crs_conversion.c"

#undef JMTX_INPUT_TYPE
#undef JMTX_OUTPUT_TYPE
#undef JMTX_NAME_IN
#undef JMTX_NAME_OUT
#undef JMTX_OUTPUT_ENUM
#undef JMTX_NAME_CONVERSION
#undef JMTX_BRM_CONVERSION_H

// cfloat -> cdouble

#define JMTX_INPUT_TYPE _Complex float
#define JMTX_OUTPUT_TYPE _Complex double
#define JMTX_NAME_IN(x) jmtxc_##x
#define JMTX_NAME_OUT(x) jmtxz_##x
#define JMTX_OUTPUT_ENUM(t) JMTXZ_TYPE_##t
#define JMTX_NAME_CONVERSION(x) jmtx_convert_##x##_c2z

#include "../../include/jmtx/jmtxc.h"
#include "../../include/jmtx/jmtxz.h"
#include "jmtx_int_defines.h"
#include "../conversion/brm_conversion.c"
#include "../conversion/ccs_conversion.c"
#include "../conversion/cds_conversion.c"
#include "../conversion/crs_conversion.c"

#undef JMTX_INPUT_TYPE
#undef JMTX_OUTPUT_TYPE
#undef JMTX_NAME_IN
#undef JMTX_NAME_OUT
#undef JMTX_OUTPUT_ENUM
#undef JMTX_NAME_CONVERSION
#undef JMTX_BRM_CONVERSION_H

// cdouble -> double

#define JMTX_INPUT_TYPE _Complex double
#define JMTX_OUTPUT_TYPE double
#define JMTX_NAME_IN(x) jmtxz_##x
#define JMTX_NAME_OUT(x) jmtxd_##x
#define JMTX_OUTPUT_ENUM(t) JMTXD_TYPE_##t
#define JMTX_NAME_CONVERSION(x) jmtx_convert_##x##_z2d

#include "../../include/jmtx/jmtxz.h"
#include "../../include/jmtx/jmtxd.h"
#include "jmtx_int_defines.h"
#include "../conversion/brm_conversion.c"
#include "../conversion/ccs_conversion.c"
#include "../conversion/cds_conversion.c"
#include "../conversion/crs_conversion.c"

#undef JMTX_INPUT_TYPE
#undef JMTX_OUTPUT_TYPE
#undef JMTX_NAME_IN
#undef JMTX_NAME_OUT
#undef JMTX_OUTPUT_ENUM
#undef JMTX_NAME_CONVERSION
#undef JMTX_BRM_CONVERSION_H

// double -> cdouble

#define JMTX_INPUT_TYPE double
#define JMTX_OUTPUT_TYPE _Complex double
#define JMTX_NAME_IN(x) jmtxd_##x
#define JMTX_NAME_OUT(x) jmtxz_##x
#define JMTX_OUTPUT_ENUM(t) JMTXZ_TYPE_##t
#define JMTX_NAME_CONVERSION(x) jmtx_convert_##x##_d2z

#include "../../include/jmtx/jmtxd.h"
#include "../../include/jmtx/jmtxz.h"
#include "jmtx_int_defines.h"
#include "../conversion/brm_conversion.c"
#include "../conversion/ccs_conversion.c"
#include "../conversion/cds_conversion.c"
#include "../conversion/crs_conversion.c"

#undef JMTX_INPUT_TYPE
#undef JMTX_OUTPUT_TYPE
#undef JMTX_NAME_IN
#undef JMTX_NAME_OUT
#undef JMTX_OUTPUT_ENUM
#undef JMTX_NAME_CONVERSION
#undef JMTX_BRM_CONVERSION_H

// cdouble -> cfloat

#define JMTX_INPUT_TYPE _Complex double
#define JMTX_OUTPUT_TYPE _Complex float
#define JMTX_NAME_IN(x) jmtxz_##x
#define JMTX_NAME_OUT(x) jmtxc_##x
#define JMTX_OUTPUT_ENUM(t) JMTXC_TYPE_##t
#define JMTX_NAME_CONVERSION(x) jmtx_convert_##x##_z2c

#include "../../include/jmtx/jmtxz.h"
#include "../../include/jmtx/jmtxc.h"
#include "jmtx_int_defines.h"
#include "../conversion/brm_conversion.c"
#include "../conversion/ccs_conversion.c"
#include "../conversion/cds_conversion.c"
#include "../conversion/crs_conversion.c"

#undef JMTX_INPUT_TYPE
#undef JMTX_OUTPUT_TYPE
#undef JMTX_NAME_IN
#undef JMTX_NAME_OUT
#undef JMTX_OUTPUT_ENUM
#undef JMTX_NAME_CONVERSION
#undef JMTX_BRM_CONVERSION_H

// cfloat -> float

#define JMTX_INPUT_TYPE _Complex float
#define JMTX_OUTPUT_TYPE float
#define JMTX_NAME_IN(x) jmtxc_##x
#define JMTX_NAME_OUT(x) jmtxf_##x
#define JMTX_OUTPUT_ENUM(t) JMTXF_TYPE_##t
#define JMTX_NAME_CONVERSION(x) jmtx_convert_##x##_c2f

#include "../../include/jmtx/jmtxc.h"
#include "../../include/jmtx/jmtxf.h"
#include "jmtx_int_defines.h"
#include "../conversion/brm_conversion.c"
#include "../conversion/ccs_conversion.c"
#include "../conversion/cds_conversion.c"
#include "../conversion/crs_conversion.c"

#undef JMTX_INPUT_TYPE
#undef JMTX_OUTPUT_TYPE
#undef JMTX_NAME_IN
#undef JMTX_NAME_OUT
#undef JMTX_OUTPUT_ENUM
#undef JMTX_NAME_CONVERSION
#undef JMTX_BRM_CONVERSION_H

// float -> double

#define JMTX_INPUT_TYPE float
#define JMTX_OUTPUT_TYPE double
#define JMTX_NAME_IN(x) jmtxf_##x
#define JMTX_NAME_OUT(x) jmtxd_##x
#define JMTX_OUTPUT_ENUM(t) JMTXD_TYPE_##t
#define JMTX_NAME_CONVERSION(x) jmtx_convert_##x##_f2d

#include "../../include/jmtx/jmtxf.h"
#include "../../include/jmtx/jmtxd.h"
#include "jmtx_int_defines.h"
#include "../conversion/brm_conversion.c"
#include "../conversion/ccs_conversion.c"
#include "../conversion/cds_conversion.c"
#include "../conversion/crs_conversion.c"

#undef JMTX_INPUT_TYPE
#undef JMTX_OUTPUT_TYPE
#undef JMTX_NAME_IN
#undef JMTX_NAME_OUT
#undef JMTX_OUTPUT_ENUM
#undef JMTX_NAME_CONVERSION
#undef JMTX_BRM_CONVERSION_H

// float -> cdouble

#define JMTX_INPUT_TYPE float
#define JMTX_OUTPUT_TYPE _Complex double
#define JMTX_NAME_IN(x) jmtxf_##x
#define JMTX_NAME_OUT(x) jmtxz_##x
#define JMTX_OUTPUT_ENUM(t) JMTXZ_TYPE_##t
#define JMTX_NAME_CONVERSION(x) jmtx_convert_##x##_f2z

#include "../../include/jmtx/jmtxf.h"
#include "../../include/jmtx/jmtxz.h"
#include "jmtx_int_defines.h"
#include "../conversion/brm_conversion.c"
#include "../conversion/ccs_conversion.c"
#include "../conversion/cds_conversion.c"
#include "../conversion/crs_conversion.c"

#undef JMTX_INPUT_TYPE
#undef JMTX_OUTPUT_TYPE
#undef JMTX_NAME_IN
#undef JMTX_NAME_OUT
#undef JMTX_OUTPUT_ENUM
#undef JMTX_NAME_CONVERSION
#undef JMTX_BRM_CONVERSION_H

// cfloat -> double

#define JMTX_INPUT_TYPE _Complex float
#define JMTX_OUTPUT_TYPE double
#define JMTX_NAME_IN(x) jmtxc_##x
#define JMTX_NAME_OUT(x) jmtxd_##x
#define JMTX_OUTPUT_ENUM(t) JMTXD_TYPE_##t
#define JMTX_NAME_CONVERSION(x) jmtx_convert_##x##_c2d

#include "../../include/jmtx/jmtxc.h"
#include "../../include/jmtx/jmtxd.h"
#include "jmtx_int_defines.h"
#include "../conversion/brm_conversion.c"
#include "../conversion/ccs_conversion.c"
#include "../conversion/cds_conversion.c"
#include "../conversion/crs_conversion.c"

#undef JMTX_INPUT_TYPE
#undef JMTX_OUTPUT_TYPE
#undef JMTX_NAME_IN
#undef JMTX_NAME_OUT
#undef JMTX_OUTPUT_ENUM
#undef JMTX_NAME_CONVERSION
#undef JMTX_BRM_CONVERSION_H

// cdouble -> float

#define JMTX_INPUT_TYPE _Complex double
#define JMTX_OUTPUT_TYPE float
#define JMTX_NAME_IN(x) jmtxz_##x
#define JMTX_NAME_OUT(x) jmtxf_##x
#define JMTX_OUTPUT_ENUM(t) JMTXF_TYPE_##t
#define JMTX_NAME_CONVERSION(x) jmtx_convert_##x##_z2f

#include "../../include/jmtx/jmtxz.h"
#include "../../include/jmtx/jmtxf.h"
#include "jmtx_int_defines.h"
#include "../conversion/brm_conversion.c"
#include "../conversion/ccs_conversion.c"
#include "../conversion/cds_conversion.c"
#include "../conversion/crs_conversion.c"

#undef JMTX_INPUT_TYPE
#undef JMTX_OUTPUT_TYPE
#undef JMTX_NAME_IN
#undef JMTX_NAME_OUT
#undef JMTX_OUTPUT_ENUM
#undef JMTX_NAME_CONVERSION
#undef JMTX_BRM_CONVERSION_H

// double -> cfloat

#define JMTX_INPUT_TYPE double
#define JMTX_OUTPUT_TYPE _Complex float
#define JMTX_NAME_IN(x) jmtxd_##x
#define JMTX_NAME_OUT(x) jmtxc_##x
#define JMTX_OUTPUT_ENUM(t) JMTXC_TYPE_##t
#define JMTX_NAME_CONVERSION(x) jmtx_convert_##x##_d2c

#include "../../include/jmtx/jmtxc.h"
#include "../../include/jmtx/jmtxd.h"
#include "jmtx_int_defines.h"
#include "../conversion/brm_conversion.c"
#include "../conversion/ccs_conversion.c"
#include "../conversion/cds_conversion.c"
#include "../conversion/crs_conversion.c"

#undef JMTX_INPUT_TYPE
#undef JMTX_OUTPUT_TYPE
#undef JMTX_NAME_IN
#undef JMTX_NAME_OUT
#undef JMTX_OUTPUT_ENUM
#undef JMTX_NAME_CONVERSION
#undef JMTX_BRM_CONVERSION_H
