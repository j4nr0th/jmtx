//
// Created by jan on 13.6.2022.
//

#ifndef MTXLIB_SPARSE_ROW_COMPRESSED_H
#define MTXLIB_SPARSE_ROW_COMPRESSED_H
#include "matrix_base.h"

//  TODO: check if refactoring and getting rid of calloc broke anything, due to memory no longer being zeroed before use

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

typedef struct jmtx_matrix_crs_struct jmtx_matrix_crs;

struct jmtx_matrix_crs_struct
{
    jmtx_matrix base;
    //  How many elements exist in the rows above, so that row i is from index elements_before[i] ot elements_before[i + 1]
    uint32_t* elements_before;
    //  Column indices corresponding with the individual elements
    uint32_t* indices;
    //  Values of elements
    jmtx_scalar_t* elements;
    uint32_t n_elements;
    uint32_t capacity;
};

#define jmtx_matrix_crs_memory_usage(mtx) (sizeof(mtx) + ((mtx).capacity + 1) * (sizeof(*(mtx).indices) + sizeof(*(mtx).elements))\
+ sizeof(*(mtx).elements_before) * ((mtx).rows + 1))

/**
 * Initializes a new Compressed Row Sparse matrix
 * @param mtx pointer to memory where the matrix should be initialized
 * @param columns number of columns of the sparse matrix
 * @param rows number of rows of the sparse matrix
 * @param reserved_elements how many elements should the space be reserved for in the matrix intialliy
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result matrix_crs_new(
        jmtx_matrix_crs* mtx, uint32_t columns, uint32_t rows, uint32_t reserved_elements,
        const jmtx_allocator_callbacks* allocator_callbacks);

/**
 * Cleans up the crs matrix and frees all of its memory
 * @param mtx pointer to memory where the matrix is stored
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result matrix_crs_destroy(jmtx_matrix_crs* mtx);

/**
 * Frees up memory which the matrix is not currently using, which is was allocated in advance
 * @param mtx pointer to the memory where the matrix is stored
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result matrix_crs_shrink(jmtx_matrix_crs* mtx);

/**
 * Sets the row of the matrix. This is the most efficient way to build the matrix, as building it this way causes
 * minimum amount of memory allocation.
 * @param mtx pointer to the memory where the matrix is stored
 * @param row index of the row to set
 * @param n how many elements are in the row
 * @param indices column indices of the elements, which has to be sorted from lowest to highest
 * @param elements values of non-zero elements
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result matrix_crs_set_row(jmtx_matrix_crs* mtx, uint32_t row, uint32_t n, const uint32_t* indices, const jmtx_scalar_t* elements);

/**
 * Version of matrix_crs_set_row which does not touch the count of elements after the current row. This is useful when
 * building a new matrix, as it avoids unnecessary setting and resetting of these values. Should be called for each row
 * in order to ensure that the matrix is properly built.
 * @param mtx pointer to the memory where the matrix is stored
 * @param row index of the row to set
 * @param n how many elements are in the row
 * @param indices column indices of the elements, which has to be sorted from lowest to highest
 * @param elements values of non-zero elements
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result matrix_crs_build_row(jmtx_matrix_crs* mtx, uint32_t row, uint32_t n, const uint32_t* indices, const jmtx_scalar_t* elements);

/**
 * Returns the pointers to arrays of column indices and element values for that row
 * @param mtx pointer to the memory where the matrix is stored
 * @param row index of the row to get
 * @param n pointer which receives the number of elements in the row
 * @param p_indices pointer to column indices of elements
 * @param p_elements pointer to values of elements
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result matrix_crs_get_row(const jmtx_matrix_crs* mtx, uint32_t row, uint32_t* n, uint32_t** p_indices, jmtx_scalar_t** p_elements);

/**
 * Multiplies a dense column vector x by the sparse matrix and stores the result at y
 * @param mtx pointer to the memory where the matrix is stored
 * @param x pointer to vector to be multiplied
 * @param y pointer to vector where the result of multiplication is to be stored
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result matrix_crs_vector_multiply(const jmtx_matrix_crs* mtx, const jmtx_scalar_t* restrict x, jmtx_scalar_t* restrict y);

/**
 * Sets a single element in the matrix. This is about as fast as setting the entire row of the matrix at once, if the
 * element was previously zero
 * @param mtx pointer to the memory where the matrix is stored
 * @param i row index
 * @param j column index
 * @param x value to which the value is set
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result matrix_crs_set_element(jmtx_matrix_crs* mtx, uint32_t i, uint32_t j, jmtx_scalar_t x);

/**
 * Returns a single element from the matrix.
 * @param mtx pointer to the memory where the matrix is stored
 * @param i row index
 * @param j column index
 * @param x pointer which receives the value
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result matrix_crs_get_element(const jmtx_matrix_crs* mtx, uint32_t i, uint32_t j, jmtx_scalar_t* x);

/**
 * Checks if the matrix for errors and misplaced elements. Essentially checks if there is 0xDEADBEEF in the matrix
 * @param mtx pointer to the memory where the matrix is stored
 * @param p_beef_status pointer which receives the value of the beef_status (high 4 bytes is the beef count, low 4 bytes
 * are BEEF)
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result matrix_crs_beef_check(const jmtx_matrix_crs* mtx, int* p_beef_status);

/**
 * Applies a unary function on the sparse matrix, only on its stored entries, which can be modified. If the user given
 * user function returns a non-zero value, that value is returned from the function and the iteration is stopped
 * @param mtx pointer to the memory where the matrix is stored
 * @param unary_fn user given function that is to be applied on the elements
 * @param param optional parameter, which is passed to the unary_fn when called
 * @return JMTX_RESULT_SUCCESS if successful, if the user function returns non-zero, that value is returned instead
 */
jmtx_result matrix_crs_apply_unary_fn(const jmtx_matrix_crs* mtx, int (*unary_fn)(uint32_t i, uint32_t j, jmtx_scalar_t* p_element, void* param), void* param);

/**
 * Removes elements exactly equal to zero. If the element is indeed zero, it is compared to (scalar_t)0.0
 * @param mtx pointer to the memory where the matrix is stored
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result matrix_crs_remove_zeros(jmtx_matrix_crs* mtx);

/**
 * Removes elements which have absolute value less than specified value
 * @param mtx pointer to the memory where the matrix is stored
 * @param v value to which to compare it to
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result matrix_crs_remove_bellow(jmtx_matrix_crs* mtx, jmtx_scalar_t v);




/**
 * Returns the number of elements in the column of the matrix
 * @param mtx pointer to the memory where the matrix is stored
 * @param col column index of the matrix to look at
 * @param p_n pointer to integer which receives the number of elements in the column
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result matrix_crs_elements_in_column(const jmtx_matrix_crs* mtx, uint32_t col, uint32_t* p_n);

/**
 * Returns the values of elements in the matrix, along with what row of the matrix they were located in
 * @param mtx pointer to the memory where the matrix is stored
 * @param col column index of the matrix to look at
 * @param n number of elements in the column to be extracted
 * @param p_elements a buffer of at least n elements which receives the values of the column
 * @param p_rows a buffer of at least n elements which receives the row indices of the column
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result matrix_crs_get_column(const jmtx_matrix_crs* mtx, uint32_t col, uint32_t n, jmtx_scalar_t* p_elements, uint32_t* p_rows);

/**
 * Creates a transpose of a matrix
 * @param mtx pointer to the memory where the input matrix is stored
 * @param out pointer to the memory where the output matrix will be stored
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result matrix_crs_transpose(const jmtx_matrix_crs* restrict mtx, jmtx_matrix_crs* restrict out);

/**
 * Creates a copy of the matrix
 * @param mtx pointer to the memory where the input matrix is stored
 * @param out pointer to the memory where the output matrix will be stored
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result matrix_crs_copy(const jmtx_matrix_crs* restrict mtx, jmtx_matrix_crs* restrict out);

/**
 * Computes one entry of Ax. This function only computes the i-th entry to make it possible to compute it in parallel
 * @param mtx pointer to the memory where the matrix A is stored
 * @param x pointer to the memory where the vector x is stored
 * @param i what entry of the residual to compute
 * @param p_r pointer to the memory which receives the computed residual
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result matrix_crs_vector_multiply_row(const jmtx_matrix_crs* mtx, const jmtx_scalar_t* x, uint32_t i, jmtx_scalar_t* p_r);

#endif //MTXLIB_SPARSE_ROW_COMPRESSED_H