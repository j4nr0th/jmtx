#ifndef JMTXI_H
#define JMTXI_H

#include "../../source/modules/jmtx_integer_defines.h"

/// Common
#include "../../source/solver_base.h"

// Includes
/// Matrix includes
#include "../../source/matrices/band_row_major.h"
#include "../../source/matrices/dense_row_major.h"
#include "../../source/matrices/sparse_column_compressed.h"
#include "../../source/matrices/sparse_conversion.h"
#include "../../source/matrices/sparse_diagonal_compressed.h"
#include "../../source/matrices/sparse_multiplication.h"
#include "../../source/matrices/sparse_row_compressed.h"

// Clear the internal header guards
/// Common
#undef JMTX_SOLVER_BASE_H
/// Matrices
#undef JMTX_BAND_ROW_MAJOR_H
#undef JMTX_DENSE_ROW_MAJOR_H
#undef JMTX_SPARSE_COLUMN_COMPRESSED_H
#undef JMTX_SPARSE_CONVERSION_H
#undef JMTX_SPARSE_DIAGONAL_COMPRESSED_H
#undef JMTX_SPARSE_MULTIPLICATION_H
#undef JMTX_SPARSE_ROW_COMPRESSED_H

#endif // JMTXI_H
