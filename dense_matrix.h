#include "matrix_base.h"

typedef struct struct_Dense_Matrix DenseMatrix;
typedef struct struct_Dense_Matrix_View *DenseMatrixView;

#define MATRIX_DENSE_COL_MAJOR 0x0400
#define MATRIX_DENSE_ROW_MAJOR 0x0401

struct struct_Dense_Matrix_View
{
    MATRIX_VIEW_STRUCT_BASE
    uint offset;
};

struct struct_Dense_Matrix
{
    MATRIX_STRUCT_BASE
    scalar_t* data;
};

DenseMatrix matrix_dense_create(uint i, uint j, int col_major);

void matrix_dense_set_element(DenseMatrix* this, uint i, uint j, scalar_t x);

void matrix_dense_row_elements(DenseMatrix* this, uint i, scalar_t* x);

void matrix_dense_col_elements(DenseMatrix* this, uint j, scalar_t* x);

MatrixIterator matrix_dense_flat(DenseMatrix* this);

MatrixIterator matrix_dense_row_iter(DenseMatrix* this);

MatrixIterator matrix_dense_col_iter(DenseMatrix* this);

MatrixIterator matrix_dense_transpose_iter(DenseMatrix* this);

int matrix_dense_row_major(DenseMatrix* this);
