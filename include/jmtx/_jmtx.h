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
/// Decompositions
#include "../../source/decompositions/band_lu_decomposition.h"
#include "../../source/decompositions/dense_lu_decomposition.h"
#include "../../source/decompositions/incomplete_cholesky_decomposition.h"
#include "../../source/decompositions/incomplete_lu_decomposition.h"
/// Solver includes
#include "../../source/solvers/bicgstab_iteration.h"
#include "../../source/solvers/cholesky_solving.h"
#include "../../source/solvers/conjugate_gradient_iteration.h"
#include "../../source/solvers/gauss_seidel_iteration.h"
#include "../../source/solvers/generalized_minimum_residual_iteration.h"
#include "../../source/solvers/jacobi_point_iteration.h"
#include "../../source/solvers/lu_solving.h"
#include "../../source/solvers/recursive_generalized_minimum_residual_iteration.h"

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
/// Decompositions
#undef JMTX_BAND_LU_DECOMPOSITION_H
#undef JMTX_DENSE_LU_DECOMPOSITION_H
#undef JMTX_INCOMPLETE_CHOLESKY_DECOMPOSITION_H
#undef JMTX_INCOMPLETE_LU_DECOMPOSITION_H
/// Solvers
#undef JMTX_JMTX_FLOAT_DEFINES_H
#undef JMTX_BICGSTAB_ITERATION_H
#undef JMTX_CHOLESKY_SOLVING_H
#undef JMTX_CONJUGATE_GRADIENT_ITERATION_H
#undef JMTX_GAUSS_SEIDEL_ITERATION_H
#undef JMTX_GENERALIZED_MINIMUM_RESIDUAL_H
#undef JMTX_JACOBI_POINT_ITERATION_H
#undef JMTX_LU_SOLVING_H
#undef JMTX_RECURSIVE_GENERALIZED_MINIMUM_RESIDUAL_ITERATION_H

// Clear the type defines
#include "../../source/modules/jmtx_clear_defines.h"
