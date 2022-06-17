//
// Created by jan on 13.6.2022.
//

#ifndef MTXLIB_SPARSE_ROW_COMPRESSED_H
#define MTXLIB_SPARSE_ROW_COMPRESSED_H
#include "matrix_base.h"

//  IMPORTANT:
//  These functions are not safe in the slightest. They perform zero checking at the moment (should probably be done as
//  with macros if _DEBUG is defined and/or NDEBUG is not defined

//  Function testing 14.06.2022:
//  - matrix_crs_new : DONE
//  - matrix_crs_destroy : DONE
//  - matrix_crs_shrink : DONE
//  - matrix_crs_set_row : DONE
//  - matrix_crs_get_row : DONE
//  - matrix_crs_vector_multiply : DONE
//  - matrix_crs_set_element : DONE :) 15.06.2022
//  - matrix_crs_get_element : DONE
//  - matrix_crs_beef_check : BEEF (could also replace with B16B00B5 ( ͡° ͜ʖ ͡°))
//  - matrix_crs_apply_unary_fn : DONE
//  - matrix_crs_remove_zeros : DONE
//  - matrix_crs_remove_bellow : DONE
//  All done on 14.06.2022

//  Function testing 15.06.2022:
//  - matrix_crs_elements_in_column : DONE
//  - matrix_crs_get_column : DONE
//  - matrix_crs_transpose : DONE
//  - matrix_crs_transpose : DONE
//  - matrix_crs_copy : DONE

typedef struct struct_CRS_Matrix CrsMatrix;

struct struct_CRS_Matrix
{
    MATRIX_STRUCT_BASE
    //  How many elements exist in the rows above, so that row i is from index elements_before[i] ot elements_before[i + 1]
    uint* elements_before;
    //  Column indices corresponding with the individual elements
    uint* indices;
    //  Values of elements
    scalar_t* elements;
    uint capacity;
};

#define matrix_crs_memory_usage(mtx) (sizeof(mtx) + ((mtx).capacity + 1) * (sizeof(*(mtx).indices) + sizeof(*(mtx).elements))\
+ sizeof(*(mtx).elements_before) * ((mtx).rows + 1))

/**
 * Initializes a new Compressed Row Sparse matrix
 * @param mtx pointer to memory where the matrix should be initialized
 * @param columns number of columns of the sparse matrix
 * @param rows number of rows of the sparse matrix
 * @param reserved_elements how many elements should the space be reserved for in the matrix intialliy
 * @return zero if successful
 */
mtx_res_t matrix_crs_new(CrsMatrix* mtx, uint columns, uint rows, uint reserved_elements);

/**
 * Cleans up the crs matrix and frees all of its memory
 * @param mtx pointer to memory where the matrix is stored
 * @return zero if successful
 */
mtx_res_t matrix_crs_destroy(CrsMatrix* mtx);

/**
 * Frees up memory which the matrix is not currently using, which is was allocated in advance
 * @param mtx pointer to the memory where the matrix is stored
 * @return zero if successful
 */
mtx_res_t matrix_crs_shrink(CrsMatrix* mtx);

/**
 * Sets the row of the matrix. This is the most efficient way to build the matrix, as building it this way causes
 * minimum amount of memory allocation.
 * @param mtx pointer to the memory where the matrix is stored
 * @param row index of the row to set
 * @param n how many elements are in the row
 * @param indices column indices of the elements, which has to be sorted from lowest to highest
 * @param elements values of non-zero elements
 * @return zero if successful
 */
mtx_res_t matrix_crs_set_row(CrsMatrix* mtx, uint row, uint n, const uint* indices, const scalar_t* elements);

/**
 * Version of matrix_crs_set_row which does not touch the count of elements after the current row. This is useful when
 * building a new matrix, as it avoids unnecessary setting and resetting of these values. Should be called for each row
 * in order to ensure that the matrix is properly built.
 * @param mtx pointer to the memory where the matrix is stored
 * @param row index of the row to set
 * @param n how many elements are in the row
 * @param indices column indices of the elements, which has to be sorted from lowest to highest
 * @param elements values of non-zero elements
 * @return zero if successful
 */
mtx_res_t matrix_crs_build_row(CrsMatrix* mtx, uint row, uint n, const uint* indices, const scalar_t* elements);

/**
 * Returns the pointers to arrays of column indices and element values for that row
 * @param mtx pointer to the memory where the matrix is stored
 * @param row index of the row to get
 * @param n pointer which receives the number of elements in the row
 * @param p_indices pointer to column indices of elements
 * @param p_elements pointer to values of elements
 * @return zero if successful
 */
mtx_res_t matrix_crs_get_row(const CrsMatrix* mtx, uint row, uint* n, uint** p_indices, scalar_t** p_elements);

/**
 * Multiplies a dense column vector x by the sparse matrix and stores the result at y
 * @param mtx pointer to the memory where the matrix is stored
 * @param x pointer to vector to be multiplied
 * @param y pointer to vector where the result of multiplication is to be stored
 * @return zero if successful
 */
mtx_res_t matrix_crs_vector_multiply(const CrsMatrix* mtx, const scalar_t* restrict x, scalar_t* restrict y);

/**
 * Sets a single element in the matrix. This is about as fast as setting the entire row of the matrix at once, if the
 * element was previously zero
 * @param mtx pointer to the memory where the matrix is stored
 * @param i row index
 * @param j column index
 * @param x value to which the value is set
 * @return zero if successful
 */
mtx_res_t matrix_crs_set_element(CrsMatrix* mtx, uint i, uint j, scalar_t x);

/**
 * Returns a single element from the matrix.
 * @param mtx pointer to the memory where the matrix is stored
 * @param i row index
 * @param j column index
 * @param x pointer which receives the value
 * @return zero if successful
 */
mtx_res_t matrix_crs_get_element(const CrsMatrix* mtx, uint i, uint j, scalar_t* x);

/**
 * Checks if the matrix for errors and misplaced elements. Essentially checks if there is 0xDEADBEEF in the matrix
 * @param mtx pointer to the memory where the matrix is stored
 * @param p_beef_status pointer which receives the value of the beef_status (high 4 bytes is the beef count, low 4 bytes
 * are BEEF)
 * @return zero if successful
 */
mtx_res_t matrix_crs_beef_check(const CrsMatrix* mtx, int* p_beef_status);

/**
 * Applies a unary function on the sparse matrix, only on its stored entries, which can be modified. If the user given
 * user function returns a non-zero value, that value is returned from the function and the iteration is stopped
 * @param mtx pointer to the memory where the matrix is stored
 * @param unary_fn user given function that is to be applied on the elements
 * @param param optional parameter, which is passed to the unary_fn when called
 * @return zero if successful, if the user function returns non-zero, that value is returned instead
 */
mtx_res_t matrix_crs_apply_unary_fn(const CrsMatrix* mtx, int (*unary_fn)(uint i, uint j, scalar_t* p_element, void* param), void* param);

/**
 * Removes elements exactly equal to zero. If the element is indeed zero, it is compared to (scalar_t)0.0
 * @param mtx pointer to the memory where the matrix is stored
 * @return zero if successful
 */
mtx_res_t matrix_crs_remove_zeros(CrsMatrix* mtx);

/**
 * Removes elements which have absolute value less than specified value
 * @param mtx pointer to the memory where the matrix is stored
 * @param v value to which to compare it to
 * @return zero if successful
 */
mtx_res_t matrix_crs_remove_bellow(CrsMatrix* mtx, scalar_t v);




/**
 * Returns the number of elements in the column of the matrix
 * @param mtx pointer to the memory where the matrix is stored
 * @param col column index of the matrix to look at
 * @param p_n pointer to integer which receives the number of elements in the column
 * @return zero if successful
 */
mtx_res_t matrix_crs_elements_in_column(const CrsMatrix* mtx, uint col, uint* p_n);

/**
 * Returns the values of elements in the matrix, along with what row of the matrix they were located in
 * @param mtx pointer to the memory where the matrix is stored
 * @param col column index of the matrix to look at
 * @param n number of elements in the column to be extracted
 * @param p_elements a buffer of at least n elements which receives the values of the column
 * @param p_rows a buffer of at least n elements which receives the row indices of the column
 * @return zero if successful
 */
mtx_res_t matrix_crs_get_column(const CrsMatrix* mtx, uint col, uint n, scalar_t* p_elements, uint* p_rows);

/**
 * Creates a transpose of a matrix
 * @param mtx pointer to the memory where the input matrix is stored
 * @param out pointer to the memory where the output matrix will be stored
 * @return zero if successful
 */
mtx_res_t matrix_crs_transpose(const CrsMatrix* restrict mtx, CrsMatrix* restrict out);

/**
 * Creates a copy of the matrix
 * @param mtx pointer to the memory where the input matrix is stored
 * @param out pointer to the memory where the output matrix will be stored
 * @return zero if successful
 */
mtx_res_t matrix_crs_copy(const CrsMatrix* restrict mtx, CrsMatrix* restrict out);

#endif //MTXLIB_SPARSE_ROW_COMPRESSED_H
