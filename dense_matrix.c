#include "dense_matrix.h"

static const IMatrix dense_row_vtbl;
static const IMatrixView dense_row_view_vtbl;
static const IMatrix dense_col_vtbl;
static const IMatrixView dense_col_view_vtbl;



static scalar_t dense_row_view_element_at(MatrixView this, uint i, uint j)
{
    DenseMatrixView view = (DenseMatrixView)this;
    const DenseMatrix* const parent = (DenseMatrix*)view->parent;
    const scalar_t* const ptr = parent->data + view->offset;
    return ptr[j + i * parent->columns];
}

static uint dense_row_view_elements_in_row(MatrixView this, uint i)
{
    DenseMatrixView view = (DenseMatrixView)this;
    return view->columns;
}

static uint dense_row_view_elements_in_column(MatrixView this, uint j)
{
    DenseMatrixView view = (DenseMatrixView)this;
    return view->rows;
}

static void dense_row_view_get_row_elements(MatrixView this, scalar_t* const elements, uint i)
{
    DenseMatrixView view = (DenseMatrixView)this;
    memcpy(elements, view->offset + ((DenseMatrix*)view->parent)->data + view->parent->columns * i, view->columns * sizeof(scalar_t));
}

static void dense_row_view_get_col_elements(MatrixView this, scalar_t* const elements, uint j)
{
    DenseMatrixView view = (DenseMatrixView)this;
    const scalar_t* const ptr = ((DenseMatrix*)view->parent)->data + view->offset;
    for (uint i = 0; i < view->rows; ++i)
    {
        elements[i] = ptr[j + i * view->parent->columns];
    }
}

static scalar_t* dense_row_view_element_ptr(MatrixView this)
{
    DenseMatrixView view = (DenseMatrixView)this;
    return ((DenseMatrix*)view->parent)->data + view->offset;
}

static void dense_row_view_get_row_indices(MatrixView this, uint* const indices, uint i)
{
    DenseMatrixView view = (DenseMatrixView)this;
    for (uint j = 0; j < view->columns; ++j)
    {
        indices[j] = j + view->offset + i * view->parent->columns;
    }
}

static void dense_row_view_get_col_indices(MatrixView this, uint* const indices, uint j)
{
    DenseMatrixView view = (DenseMatrixView)this;
    for (uint i = 0; i < view->rows; ++i)
    {
        indices[i] = j + view->offset + i * view->parent->columns;
    }
}

static int dense_row_view_element_at_if(MatrixView this, uint i, uint j, scalar_t* ptr)
{
    DenseMatrixView view = (DenseMatrixView)this;
    if (ptr)
    {
        *ptr = dense_row_view_element_at(this, i, j);
    }
    return 1;
}

static scalar_t dense_row_view_element_as_flat(MatrixView this, uint j)
{
    DenseMatrixView view = (DenseMatrixView)this;
    const uint i = j / view->columns;
    j %= view->columns;
    return dense_row_view_element_at(this, i, j);
}

static Matrix dense_row_view_copy(MatrixView this)
{
    DenseMatrixView view = (DenseMatrixView)this;
    DenseMatrix* const mtx = malloc(sizeof*mtx);
    if (!mtx) return NULL;
    memset(mtx, 0, sizeof(*mtx));
    mtx->functions = &dense_row_vtbl;
    mtx->rows = view->rows;
    mtx->columns = view->columns;
    mtx->elements = view->rows * view->columns;
    mtx->data = calloc(mtx->elements, sizeof(*mtx->data));
    if (!mtx->data)
    {
        free(mtx);
        return NULL;
    }
    for (uint i = 0; i < view->rows; ++i)
    {
        dense_row_view_get_row_elements(this, mtx->data + i * mtx->columns, i);
    }

    return (Matrix)mtx;
}

static const IMatrixView dense_row_view_vtbl =
{
    .view_type = MATRIX_DENSE_ROW_MAJOR,
    .element_at = dense_row_view_element_at,
    .elements_in_row = dense_row_view_elements_in_row,
    .elements_in_column = dense_row_view_elements_in_column,
    .get_row_elements = dense_row_view_get_row_elements,
    .get_column_elements = dense_row_view_get_col_elements,
    .get_row_indices = dense_row_view_get_row_indices,
    .get_column_indices = dense_row_view_get_col_indices,
    .element_at_if = dense_row_view_element_at_if,
    .element_as_flat = dense_row_view_element_as_flat,
    .copy = dense_row_view_copy,
};




static scalar_t dense_col_view_element_at(MatrixView this, uint i, uint j)
{
    DenseMatrixView view = (DenseMatrixView)this;
    const DenseMatrix* const parent = (DenseMatrix*)view->parent;
    const scalar_t* const ptr = parent->data + view->offset;
    return ptr[i + j * parent->rows];
}

static uint dense_col_view_elements_in_row(MatrixView this, uint i)
{
    DenseMatrixView view = (DenseMatrixView)this;
    return view->columns;
}

static uint dense_col_view_elements_in_column(MatrixView this, uint j)
{
    DenseMatrixView view = (DenseMatrixView)this;
    return view->rows;
}

static void dense_col_view_get_col_elements(MatrixView this, scalar_t* const elements, uint i)
{
    DenseMatrixView view = (DenseMatrixView)this;
    memcpy(elements, view->offset + ((DenseMatrix*)view->parent)->data + view->parent->rows * i, view->rows * sizeof(scalar_t));
}

static void dense_col_view_get_row_elements(MatrixView this, scalar_t* const elements, uint j)
{
    DenseMatrixView view = (DenseMatrixView)this;
    const scalar_t* const ptr = ((DenseMatrix*)view->parent)->data + view->offset;
    for (uint i = 0; i < view->rows; ++i)
    {
        elements[i] = ptr[i + j * view->parent->rows];
    }
}

static scalar_t* dense_col_view_element_ptr(MatrixView this)
{
    DenseMatrixView view = (DenseMatrixView)this;
    return ((DenseMatrix*)view->parent)->data + view->offset;
}

static void dense_col_view_get_row_indices(MatrixView this, uint* const indices, uint i)
{
    DenseMatrixView view = (DenseMatrixView)this;
    for (uint j = 0; j < view->columns; ++j)
    {
        indices[j] = i + view->offset + j * view->parent->rows;
    }
}

static void dense_col_view_get_col_indices(MatrixView this, uint* const indices, uint j)
{
    DenseMatrixView view = (DenseMatrixView)this;
    for (uint i = 0; i < view->rows; ++i)
    {
        indices[i] = i + view->offset + j * view->parent->rows;
    }
}

static int dense_col_view_element_at_if(MatrixView this, uint i, uint j, scalar_t* ptr)
{
    DenseMatrixView view = (DenseMatrixView)this;
    if (ptr)
    {
        *ptr = dense_col_view_element_at(this, i, j);
    }
    return 1;
}

static scalar_t dense_col_view_element_as_flat(MatrixView this, uint j)
{
    DenseMatrixView view = (DenseMatrixView)this;
    const uint i = j / view->rows;
    j %= view->rows;
    return dense_col_view_element_at(this, i, j);
}

static Matrix dense_col_view_copy(MatrixView this)
{
    DenseMatrixView view = (DenseMatrixView)this;
    DenseMatrix* const mtx = malloc(sizeof*mtx);
    if (!mtx) return NULL;
    memset(mtx, 0, sizeof(*mtx));
    mtx->functions = &dense_col_vtbl;
    mtx->rows = view->rows;
    mtx->columns = view->columns;
    mtx->elements = view->rows * view->columns;
    mtx->data = calloc(mtx->elements, sizeof(*mtx->data));
    if (!mtx->data)
    {
        free(mtx);
        return NULL;
    }
    for (uint i = 0; i < view->columns; ++i)
    {
        dense_col_view_get_col_elements(this, mtx->data + i * mtx->rows, i);
    }

    return (Matrix)mtx;
}

static const IMatrixView dense_col_view_vtbl =
        {
                .view_type = MATRIX_DENSE_COL_MAJOR,
                .element_at = dense_col_view_element_at,
                .elements_in_row = dense_col_view_elements_in_row,
                .elements_in_column = dense_col_view_elements_in_column,
                .get_row_elements = dense_col_view_get_row_elements,
                .get_column_elements = dense_col_view_get_col_elements,
                .get_row_indices = dense_col_view_get_row_indices,
                .get_column_indices = dense_col_view_get_col_indices,
                .element_at_if = dense_col_view_element_at_if,
                .element_as_flat = dense_col_view_element_as_flat,
                .copy = dense_col_view_copy,
        };




static scalar_t dense_row_element_at(Matrix mtx, uint i, uint j)
{
    const DenseMatrix* const this = (DenseMatrix*)mtx;
    return this->data[j + i * this->rows];
}

static void dense_row_get_row_elements(Matrix mtx, scalar_t* const elements, uint i)
{
    const DenseMatrix* const this = (DenseMatrix*)mtx;
    memcpy(elements, this->data + i * this->columns, this->columns * sizeof(*this->data));
}

static void dense_row_get_col_elements(Matrix mtx, scalar_t* const elements, uint j)
{
    const DenseMatrix* const this = (DenseMatrix*)mtx;
    for (uint i = 0; i < this->rows; ++i)
    {
        elements[i] = this->data[j + this->columns * i];
    }
}

static void dense_row_get_row_indices(Matrix mtx, uint* const indices, uint i)
{
    const DenseMatrix* const this = (DenseMatrix*)mtx;
    for (uint j = 0; j < this->columns; ++j)
    {
        indices[j] = j + i * this->columns;
    }
}

static void dense_row_get_col_indices(Matrix mtx, uint* const indices, uint j)
{
    const DenseMatrix* const this = (DenseMatrix*)mtx;
    for (uint i = 0; i < this->rows; ++i)
    {
        indices[i] = j + i * this->columns;
    }
}

static int dense_row_element_at_if(Matrix mtx, uint i, uint j, scalar_t* ptr)
{
    const DenseMatrix* const this = (DenseMatrix*)mtx;
    if (ptr)
    {
        *ptr = this->data[j + i * this->columns];
    }
    return 1;
}

static Matrix dense_copy(Matrix mtx)
{
    const DenseMatrix* const this = (DenseMatrix*)mtx;
    DenseMatrix* const new_mtx = malloc(sizeof*new_mtx);
    if (!new_mtx) return NULL;
    scalar_t* const new_ptr = calloc(this->elements, sizeof* new_ptr);
    if (!new_ptr)
    {
        free(new_mtx);
        return NULL;
    }

//    new_mtx->columns = this->columns;
//    new_mtx->rows = this->rows;
//    new_mtx->elements = this->elements;
//    new_mtx->functions = this->functions;
    memcpy(new_mtx, this, sizeof(*this));
    new_mtx->data = new_ptr;
    memcpy(new_ptr, this->data, sizeof*new_ptr * this->elements);
    return (Matrix)new_mtx;
}

//static Matrix dense_add_scalar(Matrix mtx, scalar_t v)
//{
//    DenseMatrix* const this = (DenseMatrix*)dense_copy(mtx);
//    if (!this) return NULL;
//    for (uint i = 0; i < this->elements; ++i)
//    {
//        this->data[i] += v;
//    }
//}
//
//static Matrix dense_add_scalar_inplace(Matrix mtx, scalar_t v)
//{
//    DenseMatrix* const this = (DenseMatrix*)(mtx);
//    for (uint i = 0; i < this->elements; ++i)
//    {
//        this->data[i] += v;
//    }
//}
//
//static Matrix dense_mul_scalar(Matrix mtx, scalar_t v)
//{
//    DenseMatrix* const this = (DenseMatrix*)dense_copy(mtx);
//    if (!this) return NULL;
//    for (uint i = 0; i < this->elements; ++i)
//    {
//        this->data[i] *= v;
//    }
//}
//
//static Matrix dense_mul_scalar_inplace(Matrix mtx, scalar_t v)
//{
//    DenseMatrix* const this = (DenseMatrix*)(mtx);
//    for (uint i = 0; i < this->elements; ++i)
//    {
//        this->data[i] *= v;
//    }
//}
//
//static Matrix dense_row_add(Matrix mtx, Matrix other)
//{
//    DenseMatrix* res = (DenseMatrix*)dense_copy(mtx);
//    if (!res) return NULL;
//    if (other->functions == mtx->functions)
//    {
//        //  Assume both have to be the same type of dense matrices
//        const DenseMatrix* const b = (DenseMatrix*)other;
//        for (uint i = 0; i < res->elements; ++i)
//        {
//            res->data[i] += b->data[i];
//        }
//    }
//    else
//    {
//        for (uint i = 0; i < res->rows; ++i)
//        {
//            for (uint j = 0; j < res->columns; ++j)
//            {
//                res->data[j + i * res->rows] += other->functions->element_at(other, i, j);
//            }
//        }
//    }
//    return (Matrix)res;
//}
//
//static Matrix dense_row_add_inplace(Matrix mtx, Matrix other)
//{
//    DenseMatrix* res = (DenseMatrix*)(mtx);
//    if (other->functions == mtx->functions)
//    {
//        //  Assume both have to be the same type of dense matrices
//        const DenseMatrix* const b = (DenseMatrix*)other;
//        for (uint i = 0; i < res->elements; ++i)
//        {
//            res->data[i] += b->data[i];
//        }
//    }
//    else
//    {
//        for (uint i = 0; i < res->rows; ++i)
//        {
//            for (uint j = 0; j < res->columns; ++j)
//            {
//                res->data[j + i * res->rows] += other->functions->element_at(other, i, j);
//            }
//        }
//    }
//    return mtx;
//}
//
//static Matrix dense_row_sub(Matrix mtx, Matrix other)
//{
//    DenseMatrix* res = (DenseMatrix*)dense_copy(mtx);
//    if (!res) return NULL;
//    if (other->functions == mtx->functions)
//    {
//        //  Assume both have to be the same type of dense matrices
//        const DenseMatrix* const b = (DenseMatrix*)other;
//        for (uint i = 0; i < res->elements; ++i)
//        {
//            res->data[i] -= b->data[i];
//        }
//    }
//    else
//    {
//        for (uint i = 0; i < res->rows; ++i)
//        {
//            for (uint j = 0; j < res->columns; ++j)
//            {
//                res->data[j + i * res->rows] -= other->functions->element_at(other, i, j);
//            }
//        }
//    }
//    return (Matrix)res;
//}
//
//static Matrix dense_row_sub_inplace(Matrix mtx, Matrix other)
//{
//    DenseMatrix* res = (DenseMatrix*)(mtx);
//    if (other->functions == mtx->functions)
//    {
//        //  Assume both have to be the same type of dense matrices
//        const DenseMatrix* const b = (DenseMatrix*)other;
//        for (uint i = 0; i < res->elements; ++i)
//        {
//            res->data[i] -= b->data[i];
//        }
//    }
//    else
//    {
//        for (uint i = 0; i < res->rows; ++i)
//        {
//            for (uint j = 0; j < res->columns; ++j)
//            {
//                res->data[j + i * res->rows] -= other->functions->element_at(other, i, j);
//            }
//        }
//    }
//    return mtx;
//}
//
//static Matrix dense_row_mul(Matrix mtx, Matrix other)
//{
//    DenseMatrix res = matrix_dense_create(mtx->rows, other->columns, 0);
//    DenseMatrix* const this = malloc(sizeof res);
//    if (!this) return NULL;
//    memcpy(this, &res, sizeof res);
//    //  mtx->columns == other->rows
//    uint* const indices = calloc(2 * mtx->columns, sizeof*indices);
//    if (!indices)
//    {
//        free(this->data);
//        free(this);
//        return NULL;
//    }
//    uint* const indices2 = indices + mtx->columns;
//    const scalar_t* const ptr1 = mtx->functions->elements_ptr(mtx);
//    const scalar_t* const ptr2 = other->functions->elements_ptr(other);
//    for (uint i = 0; i < mtx->rows; ++i)
//    {
//        mtx->functions->get_row_indices(mtx, indices, i);
//        for (uint j = 0; j < other->columns; ++j)
//        {
//            other->functions->get_column_indices(other, indices2, j);
//            scalar_t v = 0.0f;
//            for (uint k = 0; k < mtx->columns; ++k)
//            {
//                v += ptr1[indices[k]] * ptr2[indices2[k]];
//            }
//            this->data[j + i * mtx->columns] = v;
//        }
//    }
//
//    free(indices);
//    return (Matrix)this;
//}
//
//static Matrix dense_row_mul_inplace(Matrix mtx, Matrix other)
//{
//    //  works only when both mtx and other are both square matrices
//    DenseMatrix* this = (DenseMatrix*)(mtx);
//
//    scalar_t* const data = calloc(mtx->elements, sizeof*data);
//    if (!data) return NULL;
//
//    //  mtx->columns == other->rows
//    uint* const indices = calloc(2 * mtx->columns, sizeof*indices);
//    if (!indices)
//    {
//        free(data);
//        return NULL;
//    }
//    uint* const indices2 = indices + mtx->columns;
//    const scalar_t* const ptr1 = mtx->functions->elements_ptr(mtx);
//    const scalar_t* const ptr2 = other->functions->elements_ptr(other);
//    for (uint i = 0; i < mtx->rows; ++i)
//    {
//        mtx->functions->get_row_indices(mtx, indices, i);
//        for (uint j = 0; j < other->columns; ++j)
//        {
//            other->functions->get_column_indices(other, indices2, j);
//            scalar_t v = 0.0f;
//            for (uint k = 0; k < mtx->columns; ++k)
//            {
//                v += ptr1[indices[k]] * ptr2[indices2[k]];
//            }
//            data[j + i * mtx->columns] = v;
//        }
//    }
//    free(this->data);
//    this->data = data;
//    free(indices);
//    return (Matrix)this;
//}

static MatrixView dense_row_slice(Matrix this, uint begin_row, uint end_row, uint begin_col, uint end_col)
{
    DenseMatrixView view = malloc(sizeof *view);
    if (!view) return NULL;
    view->parent = this;
    view->columns = end_col - begin_col;
    view->rows = end_row - begin_row;
    view->offset = begin_col + begin_row * ((DenseMatrix*)this)->columns;
    view->functions = &dense_row_view_vtbl;
    return (MatrixView)view;
}

static MatrixView dense_col_slice(Matrix this, uint begin_row, uint end_row, uint begin_col, uint end_col)
{
    DenseMatrixView view = malloc(sizeof *view);
    if (!view) return NULL;
    view->parent = this;
    view->columns = end_col - begin_col;
    view->rows = end_row - begin_row;
    view->offset = begin_row + begin_col * ((DenseMatrix*)this)->rows;
    view->functions = &dense_col_view_vtbl;
    return (MatrixView)view;
}

static MatrixView dense_row_view(Matrix this)
{
    DenseMatrixView view = malloc(sizeof *view);
    if (!view) return NULL;
    view->parent = this;
    view->columns = this->columns;
    view->rows = this->rows;
    view->offset = 0;
    view->functions = &dense_row_view_vtbl;
    return (MatrixView)view;
}

static MatrixView dense_col_view(Matrix this)
{
    DenseMatrixView view = malloc(sizeof *view);
    if (!view) return NULL;
    view->parent = this;
    view->columns = this->columns;
    view->rows = this->rows;
    view->offset = 0;
    view->functions = &dense_col_view_vtbl;
    return (MatrixView)view;
}


static const IMatrix dense_row_vtbl =
        {
        .matrix_type = MATRIX_DENSE_ROW_MAJOR,
        };


