//
// Created by jan on 13.6.2022.
//

#ifndef JMTX_SPARSE_ROW_COMPRESSED_H
#define JMTX_SPARSE_ROW_COMPRESSED_H
#include "matrix_base.h"

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

//  Function refactoring on 21.7.2023

typedef struct jmtx_matrix_crs_struct jmtx_matrix_crs;


/**
 * Initializes a new Compressed Row Sparse matrix
 * @param p_mtx address that receives the pointer to the matrix
 * @param cols number of columns of the sparse matrix
 * @param rows number of rows of the sparse matrix
 * @param reserved_entries how many entries should the space be reserved for in the matrix initially
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result jmtx_matrix_crs_new(
        jmtx_matrix_crs** p_mtx, uint32_t cols, uint32_t rows, uint32_t reserved_entries,
        const jmtx_allocator_callbacks* allocator_callbacks);

/**
 * Cleans up the crs matrix and frees all of its memory
 * @param mtx pointer to memory where the matrix is stored
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result jmtx_matrix_crs_destroy(jmtx_matrix_crs* mtx);

/**
 * Frees up memory which the matrix is not currently using, which is was allocated in advance
 * @param mtx pointer to the memory where the matrix is stored
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result jmtx_matrix_crs_shrink(jmtx_matrix_crs* mtx);

/**
 * Sets the row of the matrix. More efficient than setting it element by element
 * @param mtx pointer to the memory where the matrix is stored
 * @param row index of the row to set
 * @param n how many entries are in the row
 * @param indices column indices of the entries, which has to be sorted from lowest to highest
 * @param values values of non-zero entries
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result jmtx_matrix_crs_set_row(jmtx_matrix_crs* mtx, uint32_t row, uint32_t n, const uint32_t* indices, const float* values);

/**
 * Version of jmtx_matrix_crs_set_row which does not touch the count of entries after the current row. This is useful when
 * building a new matrix, as it avoids unnecessary setting and resetting of these entries. Must be called for each row
 * in order to ensure that the matrix is properly built. Makes no checks on the input parameters
 * @param mtx pointer to the memory where the matrix is stored
 * @param row index of the row to set
 * @param n how many entries are in the row
 * @param indices column indices of the values, which has to be sorted from lowest to highest
 * @param values values of non-zero values
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtx_matrix_crs_build_row(jmtx_matrix_crs* mtx, uint32_t row, uint32_t n, const uint32_t* indices, const float* values);

/**
 * Returns the pointers to arrays of column indices and element values for that row
 * @param mtx pointer to the memory where the matrix is stored
 * @param row index of the row to get
 * @param n pointer which receives the number of elements in the row
 * @param p_indices pointer to column indices of elements
 * @param p_elements pointer to values of elements
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result jmtx_matrix_crs_get_row(const jmtx_matrix_crs* mtx, uint32_t row, uint32_t* n, uint32_t** p_indices, float** p_elements);

/**
 * Multiplies a dense column vector x by the sparse matrix and stores the result at y
 * @param mtx pointer to the memory where the matrix is stored
 * @param x pointer to vector to be multiplied
 * @param y pointer to vector where the result of multiplication is to be stored
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result jmtx_matrix_crs_vector_multiply(const jmtx_matrix_crs* mtx, const float* x, float* y);

/**
 * Sets a single entry in the matrix. This is about as fast as setting the entire row of the matrix at once, if the
 * element was previously zero
 * @param mtx pointer to the memory where the matrix is stored
 * @param i row index
 * @param j column index
 * @param value value to which the value is set
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result jmtx_matrix_crs_set_entry(jmtx_matrix_crs* mtx, uint32_t i, uint32_t j, float value);

/**
 * Returns a single entry from the matrix.
 * @param mtx pointer to the memory where the matrix is stored
 * @param i row index
 * @param j column index
 * @param p_value pointer which receives the value
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result jmtx_matrix_crs_get_entry(const jmtx_matrix_crs* mtx, uint32_t i, uint32_t j, float* p_value);

/**
 * Adds a value to an entry in the matrix when it exists or sets it to that value if it does not. This is about as
 * fast as setting the entire row of the matrix at once, if the entry was non-existent.
 * @param mtx pointer to the memory where the matrix is stored
 * @param i row index
 * @param j column index
 * @param value value to which the value is to be added
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result jmtx_matrix_crs_add_to_entry(jmtx_matrix_crs* mtx, uint32_t i, uint32_t j, float value);

/**
 * Checks if the matrix for errors and misplaced entries. Essentially checks if there is 0xDEADBEEF in the matrix
 * @param mtx pointer to the memory where the matrix is stored
 * @param p_beef_status pointer which receives the value of the beef_status (high 4 bytes is the beef count, low 4 bytes
 * are BEEF)
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result jmtx_matrix_crs_beef_check(const jmtx_matrix_crs* mtx, int* p_beef_status);

/**
 * Applies a unary function on the sparse matrix, only on its stored entries, which can be modified. If the user given
 * user function returns a non-zero value, that value is returned from the function and the iteration is stopped
 * @param mtx pointer to the memory where the matrix is stored
 * @param unary_fn user given function that is to be applied on the elements
 * @param param optional parameter, which is passed to the unary_fn when called
 * @return JMTX_RESULT_SUCCESS if successful, if the user function returns non-zero, that value is returned instead
 */
jmtx_result jmtx_matrix_crs_apply_unary_fn(const jmtx_matrix_crs* mtx, int (*unary_fn)(uint32_t i, uint32_t j, float* p_value, void* param), void* param);

/**
 * Removes entries exactly equal to zero. If the element is indeed zero, it is compared to (scalar_t)0.0
 * @param mtx pointer to the memory where the matrix is stored
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result jmtx_matrix_crs_remove_zeros(jmtx_matrix_crs* mtx);

/**
 * Removes entries which have absolute value less than specified value
 * @param mtx pointer to the memory where the matrix is stored
 * @param v value to which to compare it to
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result jmtx_matrix_crs_remove_bellow_magnitude(jmtx_matrix_crs* mtx, float v);

/**
 * Zeros all entries within a matrix, but does not remove them in case they need to be reused
 * @param mtx matrix to zero
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result jmtx_matrix_crs_zero_all_entries(jmtx_matrix_crs* mtx);

/**
 * Similar to jmtx_matrix_crs_zero_all_entries, but slower, since it can not use memset. On the other hand, it allows for
 * the value to be other than 0
 * @param mtx matrix to set
 * @param x value to which to set all entries to
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result jmtx_matrix_crs_set_all_entries(jmtx_matrix_crs* mtx, float x);

/**
 * Keeps the memory for the matrix, but sets the entry count to 0, so that matrix can be rebuilt.
 * @param mtx matrix to clear
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result jmtx_matrix_crs_clear(jmtx_matrix_crs* mtx);

/**
 * Returns the number of entries in the column of the matrix
 * @param mtx pointer to the memory where the matrix is stored
 * @param col column index of the matrix to look at
 * @param p_n pointer to integer which receives the number of entries in the column
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result jmtx_matrix_crs_entries_in_col(const jmtx_matrix_crs* mtx, uint32_t col, uint32_t* p_n);

/**
 * Returns the values of entries in the matrix, along with what row of the matrix they were located in
 * @param mtx pointer to the memory where the matrix is stored
 * @param col column index of the matrix to look at
 * @param n number of elements in the column to be extracted
 * @param p_values a buffer of at least n elements which receives the values of the column
 * @param p_rows a buffer of at least n elements which receives the row indices of the column
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result jmtx_matrix_crs_get_col(
        const jmtx_matrix_crs* mtx, uint32_t col, uint32_t n, uint32_t* p_count, float* p_values, uint32_t* p_rows);

/**
 * Creates a transpose of a matrix
 * @param mtx pointer to the memory where the input matrix is stored
 * @param p_out address where the pointer to the output matrix will be returned
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result jmtx_matrix_crs_transpose(
        const jmtx_matrix_crs* mtx, jmtx_matrix_crs** p_out, const jmtx_allocator_callbacks* allocator_callbacks);

/**
 * Creates a copy of the matrix
 * @param mtx pointer to the memory where the input matrix is stored
 * @param p_out address where the pointer to the output matrix will be returned
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result jmtx_matrix_crs_copy(const jmtx_matrix_crs* mtx, jmtx_matrix_crs** p_out, const jmtx_allocator_callbacks* allocator_callbacks);

/**
 * Computes one entry of Ax. This function only computes the i-th entry to make it possible to compute it in parallel.
 * Unless JMTX_NO_VERIFY_PARAMS is defined, the function performs checks on parameters.
 * @param mtx pointer to the memory where the matrix A is stored.
 * @param x pointer to the memory where the vector x is stored
 * @param i what entry of the residual to compute
 * @param p_r pointer to the memory which receives the result
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result jmtx_matrix_crs_vector_multiply_row(const jmtx_matrix_crs* mtx, const float* x, uint32_t i, float* p_r);


/**
 * Computes one entry of Ax. This function only computes the i-th entry to make it possible to compute it in parallel.
 * Performs no parameter checking.
 * @param mtx pointer to the memory where the matrix A is stored.
 * @param x pointer to the memory where the vector x is stored
 * @param i what entry of the residual to compute
 * @return result of the multiplication
 */
float jmtx_matrix_crs_vector_multiply_row_raw(const jmtx_matrix_crs* mtx, const float* x, uint32_t i);

/**
 * Removes a row from a CRS matrix, reducing it's row count by 1
 * @param mtx matrix to remove the row from
 * @param row row to remove
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result jmtx_matrix_crs_remove_row(jmtx_matrix_crs* mtx, uint32_t row);

/**
 * Removes a column from a CRS matrix, reducing it's column count by 1. Slower than removing a row, so remove rows before
 * columns if both have to be done.
 * @param mtx matrix to remove the column from
 * @param col column to remove
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result jmtx_matrix_crs_remove_column(jmtx_matrix_crs* mtx, uint32_t col);

/**
 * Combines k matrices of size N_i x M into a single Sum(N_i) x M matrix by vertically stacking them
 * @param output receives pointer to the resulting matrix
 * @param allocators memory allocators to use for the output matrix
 * @param k number of the matrices
 * @param matrix_list array of matrices to be joined together. Must have the same number of columns.
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_DIMS_MISMATCH if the number of columns is not the same for all
 * input matrices, return value of jmtx_matrix_crs_new if that fails
 */
jmtx_result jmtx_matrix_crs_join_vertically(jmtx_matrix_crs** output, const jmtx_allocator_callbacks* allocators, unsigned k, const jmtx_matrix_crs** matrix_list);

#endif //JMTX_SPARSE_ROW_COMPRESSED_H
