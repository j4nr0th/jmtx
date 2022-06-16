//
// Created by jan on 15.6.2022.
//
#ifndef MTXLIB_SPARSE_COLUMN_COMPRESSED_H
#define MTXLIB_SPARSE_COLUMN_COMPRESSED_H

//  IMPORTANT:
//  These functions are not safe in the slightest. They perform zero checking at the moment (should probably be done as
//  with macros if _DEBUG is defined and/or NDEBUG is not defined

//  Function testing 15.06.2022:
//  - matrix_crs_new : DONE
//  - matrix_crs_destroy : DONE
//  - matrix_crs_shrink : DONE
//  - matrix_crs_set_row : DONE
//  - matrix_crs_get_row : DONE
//  - matrix_crs_vector_multiply : DONE
//  - matrix_crs_set_element : DONE
//  - matrix_crs_get_element : DONE
//  - matrix_crs_beef_check : BEEF (could also replace with B16B00B5 ( ͡° ͜ʖ ͡°))
//  - matrix_crs_apply_unary_fn : DONE
//  - matrix_crs_remove_zeros : DONE
//  - matrix_crs_remove_bellow : DONE
//  - matrix_crs_elements_in_column : DONE
//  - matrix_crs_get_column : DONE
//  - matrix_crs_transpose : DONE
//  - matrix_crs_transpose : DONE
//  - matrix_crs_copy : DONE
//  The same test code was run as for crs matrices


#include "matrix_base.h"

typedef struct struct_CCS_Matrix CcsMatrix;

struct struct_CCS_Matrix
{
    MATRIX_STRUCT_BASE
    //  How many elements exist in the columns left, so that clumn i is from index elements_before[i] ot elements_before[i + 1]
    uint* elements_before;
    //  Column indices corresponding with the individual elements
    uint* indices;
    //  Values of elements
    scalar_t* elements;
    uint capacity;
};

#define matrix_ccs_memory_usage(mtx) (sizeof(mtx) + ((mtx).capacity + 1) * (sizeof(*(mtx).indices) + sizeof(*(mtx).elements))\
+ sizeof(*(mtx).elements_before) * ((mtx).columns + 1))

/**
 * Initializes a new Compressed Column Sparse matrix
 * @param mtx pointer to memory where the matrix should be initialized
 * @param columns number of columns of the sparse matrix
 * @param rows number of rows of the sparse matrix
 * @param reserved_elements how many elements should the space be reserved for in the matrix intialliy
 * @return zero if successful
 */
mtx_res_t matrix_ccs_new(CcsMatrix* mtx, uint columns, uint rows, uint reserved_elements);

/**
 * Cleans up the ccs matrix and frees all of its memory
 * @param mtx pointer to memory where the matrix is stored
 * @return zero if successful
 */
mtx_res_t matrix_ccs_destroy(CcsMatrix* mtx);

/**
 * Frees up memory which the matrix is not currently using, which is was allocated in advance
 * @param mtx pointer to the memory where the matrix is stored
 * @return zero if successful
 */
mtx_res_t matrix_ccs_shrink(CcsMatrix* mtx);

/**
 * Sets the column of the matrix. This is the most efficient way to build the matrix, as building it this way causes
 * minimum amount of memory allocation.
 * @param mtx pointer to the memory where the matrix is stored
 * @param col index of the column to set
 * @param n how many elements are in the column
 * @param indices column indices of the elements, which has to be sorted from lowest to highest
 * @param elements values of non-zero elements
 * @return zero if successful
 */
mtx_res_t matrix_ccs_set_col(CcsMatrix* mtx, uint col, uint n, const uint* indices, const scalar_t* elements);

/**
 * Returns the pointers to arrays of column indices and element values for that column
 * @param mtx pointer to the memory where the matrix is stored
 * @param col index of the column to get
 * @param n pointer which receives the number of elements in the column
 * @param p_indices pointer to column indices of elements
 * @param p_elements pointer to values of elements
 * @return zero if successful
 */
mtx_res_t matrix_ccs_get_col(const CcsMatrix* mtx, uint col, uint* n, uint** p_indices, scalar_t** p_elements);

/**
 * Multiplies a dense row vector x by the sparse matrix and stores the result at y
 * @param mtx pointer to the memory where the matrix is stored
 * @param x pointer to vector to be multiplied
 * @param y pointer to vector where the result of multiplication is to be stored
 * @return zero if successful
 */
mtx_res_t matrix_ccs_vector_multiply(const CcsMatrix* mtx, const scalar_t* restrict x, scalar_t* restrict y);

/**
 * Sets a single element in the matrix. This is about as fast as setting the entire column of the matrix at once, if the
 * element was previously zero
 * @param mtx pointer to the memory where the matrix is stored
 * @param i row index
 * @param j column index
 * @param x value to which the value is set
 * @return zero if successful
 */
mtx_res_t matrix_ccs_set_element(CcsMatrix* mtx, uint i, uint j, scalar_t x);

/**
 * Returns a single element from the matrix.
 * @param mtx pointer to the memory where the matrix is stored
 * @param i row index
 * @param j column index
 * @param x pointer which receives the value
 * @return zero if successful
 */
mtx_res_t matrix_ccs_get_element(const CcsMatrix* mtx, uint i, uint j, scalar_t* x);

/**
 * Checks if the matrix for errors and misplaced elements. Essentially checks if there is 0xDEADBEEF in the matrix
 * @param mtx pointer to the memory where the matrix is stored
 * @param p_beef_status pointer which receives the value of the beef_status (high 4 bytes is the beef count, low 4 bytes
 * are BEEF)
 * @return zero if successful
 */
mtx_res_t matrix_ccs_beef_check(const CcsMatrix* mtx, int* p_beef_status);

/**
 * Applies a unary function on the sparse matrix, only on its stored entries, which can be modified. If the user given
 * user function returns a non-zero value, that value is returned from the function and the iteration is stopped
 * @param mtx pointer to the memory where the matrix is stored
 * @param unary_fn user given function that is to be applied on the elements
 * @param param optional parameter, which is passed to the unary_fn when called
 * @return zero if successful, if the user function returns non-zero, that value is returned instead
 */
mtx_res_t matrix_ccs_apply_unary_fn(const CcsMatrix* mtx, int (*unary_fn)(uint i, uint j, scalar_t* p_element, void* param), void* param);

/**
 * Removes elements exactly equal to zero. If the element is indeed zero, it is compared to (scalar_t)0.0
 * @param mtx pointer to the memory where the matrix is stored
 * @return zero if successful
 */
mtx_res_t matrix_ccs_remove_zeros(CcsMatrix* mtx);

/**
 * Removes elements which have absolute value less than specified value
 * @param mtx pointer to the memory where the matrix is stored
 * @param v value to which to compare it to
 * @return zero if successful
 */
mtx_res_t matrix_ccs_remove_bellow(CcsMatrix* mtx, scalar_t v);




/**
 * Returns the number of elements in the column of the matrix
 * @param mtx pointer to the memory where the matrix is stored
 * @param col column index of the matrix to look at
 * @param p_n pointer to integer which receives the number of elements in the column
 * @return zero if successful
 */
mtx_res_t matrix_ccs_elements_in_row(const CcsMatrix* mtx, uint col, uint* p_n);

/**
 * Returns the values of elements in the matrix, along with what column of the matrix they were located in
 * @param mtx pointer to the memory where the matrix is stored
 * @param row column index of the matrix to look at
 * @param n number of elements in the column to be extracted
 * @param p_elements a buffer of at least n elements which receives the values of the column
 * @param p_cols a buffer of at least n elements which receives the column indices of the column
 * @return zero if successful
 */
mtx_res_t matrix_ccs_get_row(const CcsMatrix* mtx, uint row, uint n, scalar_t* p_elements, uint* p_cols);

/**
 * Creates a transpose of a matrix
 * @param mtx pointer to the memory where the input matrix is stored
 * @param out pointer to the memory where the output matrix will be stored
 * @return zero if successful
 */
mtx_res_t matrix_ccs_transpose(const CcsMatrix* restrict mtx, CcsMatrix* restrict out);

/**
 * Creates a copy of the matrix
 * @param mtx pointer to the memory where the input matrix is stored
 * @param out pointer to the memory where the output matrix will be stored
 * @return zero if successful
 */
mtx_res_t matrix_ccs_copy(const CcsMatrix* restrict mtx, CcsMatrix* restrict out);
#endif //MTXLIB_SPARSE_COLUMN_COMPRESSED_H
