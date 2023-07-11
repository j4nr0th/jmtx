#include "matrix_base.h"

typedef struct jmtx_dense_matrix_struct jmtx_dense_matrix;
typedef struct jmtx_dense_matrix_view_struct jmtx_dense_matrix_view;


struct jmtx_dense_matrix_view_struct
{
    jmtx_dense_matrix* base;
    jmtx_index_t b_row, e_row;
    jmtx_index_t b_col, e_col;
    int transpose;
};

struct jmtx_dense_matrix_struct
{
    jmtx_matrix base;
    jmtx_scalar_t* data;
};

jmtx_result matrix_dense_create(jmtx_index_t rows, jmtx_index_t cols, int col_major, jmtx_dense_matrix** pp_out);

jmtx_result jmtx_matrix_dense_set_element(jmtx_dense_matrix* mtx, jmtx_index_t i, jmtx_index_t j, jmtx_scalar_t x);

jmtx_result jmtx_matrix_dense_row_elements(jmtx_dense_matrix* mtx, uint32_t i, jmtx_scalar_t* x);

jmtx_result jmtx_matrix_dense_col_elements(jmtx_dense_matrix* mtx, uint32_t j, jmtx_scalar_t* x);

//MatrixIterator matrix_dense_flat(jmtx_dense_matrix* mtx);
//
//MatrixIterator matrix_dense_row_iter(jmtx_dense_matrix* mtx);
//
//MatrixIterator matrix_dense_col_iter(jmtx_dense_matrix* mtx);
//
//MatrixIterator matrix_dense_transpose_iter(jmtx_dense_matrix* mtx);

jmtx_result jmtx_matrix_dense_is_row_major(jmtx_dense_matrix* mtx);
