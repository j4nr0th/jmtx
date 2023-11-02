//
// Created by jan on 15.6.2022.
//

#ifndef JMTX_SPARSE_COLUMN_COMPRESSED_H
#define JMTX_SPARSE_COLUMN_COMPRESSED_H
#include "matrix_base.h"

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



typedef struct jmtx_matrix_ccs_struct jmtx_matrix_ccs;

/**
 * Initializes a new Compressed Column Sparse matrix
 * @param mtx pointer to memory where the matrix should be initialized
 * @param columns number of columns of the sparse matrix
 * @param rows number of rows of the sparse matrix
 * @param reserved_entries how many values should the space be reserved for in the matrix intialliy
 * @return zero if successful
 */
jmtx_result jmtx_matrix_ccs_new(
        jmtx_matrix_ccs** mtx, uint32_t columns, uint32_t rows, uint32_t reserved_entries,
        const jmtx_allocator_callbacks* allocator_callbacks);

/**
 * Cleans up the ccs matrix and frees all of its memory
 * @param mtx pointer to memory where the matrix is stored
 * @return zero if successful
 */
jmtx_result jmtx_matrix_ccs_destroy(jmtx_matrix_ccs* mtx);

/**
 * Frees up memory which the matrix is not currently using, which is was allocated in advance
 * @param mtx pointer to the memory where the matrix is stored
 * @return zero if successful
 */
jmtx_result jmtx_matrix_ccs_shrink(jmtx_matrix_ccs* mtx);

/**
 * Sets the column of the matrix. This is the most efficient way to build the matrix, as building it this way causes
 * minimum amount of memory allocation.
 * @param mtx pointer to the memory where the matrix is stored
 * @param col index of the column to set
 * @param n how many values are in the column
 * @param indices column indices of the values, which has to be sorted from lowest to highest
 * @param values values of non-zero values
 * @return zero if successful
 */
jmtx_result jmtx_matrix_ccs_set_col(jmtx_matrix_ccs* mtx, uint32_t col, uint32_t n, const uint32_t* indices, const float* values);

/**
 * Version of jmtx_matrix_ccs_set_col which does not touch the count of entries after the current column. This is useful
 * when building a new matrix, as it avoids unnecessary setting and resetting of these entries. Must be called for each
 * column in order to ensure that the matrix is properly built. Makes no checks on the input parameters
 * @param mtx pointer to the memory where the matrix is stored
 * @param col index of the column to set
 * @param n how many entries are in the column
 * @param indices row indices of the values, which has to be sorted from lowest to highest
 * @param values values of non-zero values
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtx_matrix_ccs_build_col(jmtx_matrix_ccs* mtx, uint32_t col, uint32_t n, const uint32_t* indices, const float* values);

/**
 * Returns the pointers to arrays of column indices and element values for that column
 * @param mtx pointer to the memory where the matrix is stored
 * @param col index of the column to get
 * @param n pointer which receives the number of values in the column
 * @param p_indices pointer to column indices of values
 * @param p_elements pointer to values of values
 * @return zero if successful
 */
jmtx_result jmtx_matrix_ccs_get_col(const jmtx_matrix_ccs* mtx, uint32_t col, uint32_t* n, uint32_t** p_indices, float** p_elements);

/**
 * Multiplies a dense <b>row</b> vector x by the sparse matrix and stores the result at y
 * @param mtx pointer to the memory where the matrix is stored
 * @param x pointer to vector to be multiplied
 * @param y pointer to vector where the result of multiplication is to be stored
 * @return zero if successful
 */
jmtx_result jmtx_matrix_ccs_vector_multiply(const jmtx_matrix_ccs* mtx, const float* restrict x, float* restrict y);

/**
 * Sets a single element in the matrix. This is about as fast as setting the entire column of the matrix at once, if the
 * element was previously zero
 * @param mtx pointer to the memory where the matrix is stored
 * @param i row index
 * @param j column index
 * @param value value to which the value is set
 * @return zero if successful
 */
jmtx_result jmtx_matrix_ccs_set_entry(jmtx_matrix_ccs* mtx, uint32_t i, uint32_t j, float value);

/**
 * Returns a single element from the matrix.
 * @param mtx pointer to the memory where the matrix is stored
 * @param i row index
 * @param j column index
 * @param p_value pointer which receives the value
 * @return zero if successful
 */
jmtx_result jmtx_matrix_ccs_get_entry(const jmtx_matrix_ccs* mtx, uint32_t i, uint32_t j, float* p_value);

/**
 * Checks if the matrix for errors and misplaced values. Essentially checks if there is 0xDEADBEEF in the matrix
 * @param mtx pointer to the memory where the matrix is stored
 * @param p_beef_status pointer which receives the value of the beef_status (high 4 bytes is the beef count, low 4 bytes
 * are BEEF)
 * @return zero if successful
 */
jmtx_result jmtx_matrix_ccs_beef_check(const jmtx_matrix_ccs* mtx, int* p_beef_status);

/**
 * Applies a unary function on the sparse matrix, only on its stored entries, which can be modified. If the user given
 * user function returns a non-zero value, that value is returned from the function and the iteration is stopped
 * @param mtx pointer to the memory where the matrix is stored
 * @param unary_fn user given function that is to be applied on the values
 * @param param optional parameter, which is passed to the unary_fn when called
 * @return zero if successful, if the user function returns non-zero, that value is returned instead
 */
jmtx_result jmtx_matrix_ccs_apply_unary_fn(const jmtx_matrix_ccs* mtx, int (*unary_fn)(uint32_t i, uint32_t j, float* p_element, void* param), void* param);

/**
 * Removes values exactly equal to zero. If the element is indeed zero, it is compared to (scalar_t)0.0
 * @param mtx pointer to the memory where the matrix is stored
 * @return zero if successful
 */
jmtx_result jmtx_matrix_ccs_remove_zeros(jmtx_matrix_ccs* mtx);

/**
 * Removes values which have absolute value less than specified value
 * @param mtx pointer to the memory where the matrix is stored
 * @param v value to which to compare it to
 * @return zero if successful
 */
jmtx_result jmtx_matrix_ccs_remove_bellow(jmtx_matrix_ccs* mtx, float v);

/**
 * Zeros all values within a matrix, but does not remove them in case they need to be reused
 * @param mtx matrix to zero
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result jmtx_matrix_ccs_zero_all_elements(jmtx_matrix_ccs* mtx);

/**
 * Similar to jmtx_matrix_ccs_zero_all_elements, but slower, since it can not use memset. On the other hand, it allows for
 * the value to be other than 0
 * @param mtx matrix to set
 * @param x value to which to set to
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result jmtx_matrix_ccs_set_all_elements(jmtx_matrix_ccs* mtx, float x);

/**
 * Returns the number of values in the column of the matrix
 * @param mtx pointer to the memory where the matrix is stored
 * @param row column index of the matrix to look at
 * @param p_n pointer to integer which receives the number of values in the column
 * @return zero if successful
 */
jmtx_result jmtx_matrix_ccs_elements_in_row(const jmtx_matrix_ccs* mtx, uint32_t row, uint32_t* p_n);

/**
 * Returns the values of values in the matrix, along with what column of the matrix they were located in
 * @param mtx pointer to the memory where the matrix is stored
 * @param row column index of the matrix to look at
 * @param n number of values in the column to be extracted
 * @param p_values a buffer of at least n values which receives the values of the column
 * @param p_count a buffer of at least n values which receives the column indices of the column
 * @return zero if successful
 */
jmtx_result jmtx_matrix_ccs_get_row(
        const jmtx_matrix_ccs* mtx, uint32_t row, uint32_t n, float* p_values, uint32_t* p_count, uint32_t* p_columns);

/**
 * Creates a transpose of a matrix
 * @param mtx pointer to the memory where the input matrix is stored
 * @param out pointer to the memory where the output matrix will be stored
 * @return zero if successful
 */
jmtx_result jmtx_matrix_ccs_transpose(const jmtx_matrix_ccs* mtx, jmtx_matrix_ccs** p_out,
                                      const jmtx_allocator_callbacks* allocator_callbacks);

/**
 * Creates a copy of the matrix
 * @param mtx pointer to the memory where the input matrix is stored
 * @param out pointer to the memory where the output matrix will be stored
 * @return zero if successful
 */
jmtx_result jmtx_matrix_ccs_copy(const jmtx_matrix_ccs* restrict mtx, jmtx_matrix_ccs** p_out, const jmtx_allocator_callbacks* allocator_callbacks);
#endif //JMTX_SPARSE_COLUMN_COMPRESSED_H
