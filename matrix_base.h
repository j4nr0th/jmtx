#ifndef MTXLIB_MATRIX_BASE_H
#define MTXLIB_MATRIX_BASE_H

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

typedef uint_fast32_t mtx_res_t;

typedef float scalar_t;

typedef struct interface_Matrix IMatrix;
typedef struct struct_Matrix *Matrix;
typedef struct interface_Matrix_View IMatrixView;
typedef struct struct_Matrix_View *MatrixView;
typedef struct struct_Matrix_Iterator MatrixIterator;

struct struct_Matrix_Iterator
{
    int done;
    int(*next)(MatrixIterator* this);
    scalar_t* ptr;
    uint i;
    Matrix parent;
};

typedef enum enum_Matrix_Type
{
    mtx_type_crs, //  Compressed row sparse

    mtx_type_ccs, //  Compressed column sparse
    mtx_type_drm, //  Dense row major
    mtx_type_dcm, //  Dense column major
    mtx_type_dia, //  Diagonal
    mtx_type_drs, //  Dense row major symmetric
    mtx_type_dcs, //  Dense row major symmetric
    mtx_type_utr, //  Dense upper triagonal row major
    mtx_type_ltr, //  Dense lower triagonal row major

} MatrixType;

#define MATRIX_STRUCT_BASE uint type;\
    uint rows;\
    uint columns;\
    uint n_elements;

struct struct_Matrix
{
    MATRIX_STRUCT_BASE
};

#define MATRIX_VIEW_STRUCT_BASE const IMatrixView* functions;\
uint rows;\
uint columns;\
Matrix parent;

struct struct_Matrix_View
{
    MATRIX_VIEW_STRUCT_BASE
};

struct interface_Matrix_View
{
    int view_type;

    //  Gives value of element at position
    scalar_t(*element_at)(MatrixView this, uint i, uint j);

    //  Number of elements in row i
    uint(*elements_in_row)(MatrixView this, uint i);
    
    //  Number of elements in column
    uint(*elements_in_column)(MatrixView this, uint j);

    //  Save elements in row i to buffer *elements
    void(*get_row_elements)(MatrixView this, scalar_t* elements, uint i);
    
    //  Save elements in column j to buffer *elements
    void(*get_column_elements)(MatrixView this, scalar_t* elements, uint j);

    //  Raw pointer to the elements of the view. This can then be used together with get_row_indices and get_column_indices
    scalar_t(*element_ptr)(MatrixView this);
    
    //  Save indices of elements in row i to buffer *indices
    void(*get_row_indices)(MatrixView this, uint* indices, uint i);
    
    //  Save indices of elements in column j to buffer *indices
    void(*get_column_indices)(MatrixView this, uint* indices, uint j);
    
    //  Returns non-zero if element exists at (i, j) and stores it *ptr, otherwise returns zero and does not touch ptr
    int(*element_at_if)(MatrixView this, uint i, uint j, scalar_t* ptr);

    //  Flat access into memory as a flat buffer
    scalar_t(*element_as_flat)(MatrixView this, uint i);

    //  Copies the view to a new matrix
    Matrix(*copy)(MatrixView this);
};


struct interface_Matrix
{
    int matrix_type;

    //  Gives value of element at position
    scalar_t(*element_at)(Matrix this, uint i, uint j);

    //  Save elements in row i to buffer *elements
    void(*get_row_elements)(Matrix this, scalar_t* elements, uint i);

    //  Save elements in column j to buffer *elements
    void(*get_column_elements)(Matrix this, scalar_t* elements, uint j);

    //  Save indices of elements in row i to buffer *indices
    void(*get_row_indices)(Matrix this, uint* indices, uint i);

    //  Save indices of elements in column j to buffer *indices
    void(*get_column_indices)(Matrix this, uint* indices, uint j);

    //  Get element pointer that can be indexed into using the get_row_indices and get_column_indices
    scalar_t*(*elements_ptr)(Matrix this);

    //  Returns non-zero if element exists at (i, j) and stores it *ptr, otherwise returns zero and does not touch ptr
    int(*element_at_if)(Matrix this, uint i, uint j, scalar_t* ptr);

    //  Copies the view to a new matrix
    Matrix(*copy)(Matrix this);

    //  Gives a slice of matrix as view
    MatrixView(*slice)(Matrix this, uint begin_row, uint end_row, uint begin_col, uint end_col);

    //  Get view of the entire matrix
    MatrixView(*view)(Matrix this);
};

#define Matrix_TO_BASE(mtx) ((Matrix)mtx)
#define Matrix_TYPE(mtx) (mtx->functions.matrix_type)
#endif // !MTXLIB_MATRIX_BASE_H
