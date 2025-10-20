#include "jmtx_cfloat_defines.h"

// Matrix types
#include "../matrices/band_row_major.c"
// #include "../matrices/basic_io.c"
#include "../matrices/dense_row_major.c"
#include "../matrices/sparse_column_compressed.c"
#include "../matrices/sparse_conversion.c"
#include "../matrices/sparse_diagonal_compressed.c"
#include "../matrices/sparse_multiplication.c"
#include "../matrices/sparse_row_compressed.c"

// Decompositions
#include "../decompositions/band_lu_decomposition.c"
#include "../decompositions/dense_lu_decomposition.c"
#include "../decompositions/incomplete_cholesky_decomposition.c"
#include "../decompositions/incomplete_lu_decomposition.c"

// Solvers
#include "../solvers/bicgstab_iteration.c"
#include "../solvers/cholesky_solving.c"
#include "../solvers/conjugate_gradient_iteration.c"
#include "../solvers/gauss_seidel_iteration.c"
#include "../solvers/generalized_minimum_residual_iteration.c"
#include "../solvers/jacobi_point_iteration.c"
#include "../solvers/lu_solving.c"
#include "../solvers/recursive_generalized_minimum_residual_iteration.c"
