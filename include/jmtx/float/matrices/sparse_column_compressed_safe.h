//
// Created by jan on 15.6.2022.
//
/**
 * Functions declared here perform more checking of the parameters. Faster "unsafe" versions of these functions,
 * which do perform parameter validation are in the "sparse_row_compressed.h" header.
 *
 */

#ifndef JMTX_SPARSE_COLUMN_COMPRESSED_SAFE_H
#define JMTX_SPARSE_COLUMN_COMPRESSED_SAFE_H
#ifndef JMTX_SPARSE_COLUMN_COMPRESSED_H
    #include "sparse_column_compressed.h"
#endif


/**
 * Initializes a new Compressed Column Sparse matrix
 * @param p_mtx address that receives the pointer to the matrix
 * @param rows number of rows of the sparse matrix
 * @param cols number of columns of the sparse matrix
 * @param reserved_entries how many entries should the space be reserved for in the matrix initially
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result jmtxs_matrix_ccs_new(
 jmtx_matrix_ccs** p_mtx, uint32_t rows, uint32_t cols, uint32_t reserved_entries,
 const jmtx_allocator_callbacks* allocator_callbacks);

/**
 * Cleans up the ccs matrix and frees all of its memory
 * @param mtx pointer to memory where the matrix is stored
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result jmtxs_matrix_ccs_destroy(jmtx_matrix_ccs* mtx);

/**
 * Frees up memory which the matrix is not currently using, which is was allocated in advance
 * @param mtx pointer to the memory where the matrix is stored
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result jmtxs_matrix_ccs_shrink(jmtx_matrix_ccs* mtx);

/**
 * Sets the column of the matrix. This is the most efficient way to build the matrix, as building it this way causes
 * minimum amount of memory allocation.
 * @param mtx pointer to the memory where the matrix is stored
 * @param col index of the column to set
 * @param n how many entries are in the column
 * @param indices row indices of the entries, which has to be sorted from lowest to highest
 * @param values values of non-zero entries
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result jmtxs_matrix_ccs_set_col(jmtx_matrix_ccs* mtx, uint32_t col, uint32_t n, const uint32_t* indices, const float* values);

/**
 * Returns the pointers to arrays of column indices and element values for that column
 * @param mtx pointer to the memory where the matrix is stored
 * @param col index of the row to get
* @param n pointer that receives the number of elements in the column, which is the number of valid elements in arrays given to p_indices and
 * p_elements
 * @param p_indices pointer to row indices of values
 * @param p_elements pointer to values of values
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result jmtxs_matrix_ccs_get_col(const jmtx_matrix_ccs* mtx, uint32_t col, uint32_t* n, uint32_t** p_indices, float** p_elements);

/**
 * Multiplies a dense row vector x by the sparse matrix and stores the result at y
 * @param mtx pointer to the memory where the matrix is stored
 * @param x pointer to vector to be multiplied
 * @param y pointer to vector where the result of multiplication is to be stored
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result jmtxs_matrix_ccs_vector_multiply(const jmtx_matrix_ccs* mtx, const float* restrict x, float* restrict y);

/**
 * Sets a single entry in the matrix. This is about as fast as setting the entire column of the matrix at once, if the
 * element was previously zero
 * @param mtx pointer to the memory where the matrix is stored
 * @param i row index
 * @param j column index
 * @param value value to which the value is set
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result jmtxs_matrix_ccs_set_entry(jmtx_matrix_ccs* mtx, uint32_t i, uint32_t j, float value);

/**
 * Returns a single element from the matrix.
 * @param mtx pointer to the memory where the matrix is stored
 * @param i row index
 * @param j column index
 * @param p_value pointer which receives the value
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result jmtxs_matrix_ccs_get_entry(const jmtx_matrix_ccs* mtx, uint32_t i, uint32_t j, float* p_value);

/**
 * Counts the number of times a specific value occurs in the matrix
 * @param mtx pointer to the memory where the matrix is stored
 * @param v value which to search for
 * @param p_count pointer which receives the number of times the value appeared in the matrix
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result jmtxs_matrix_ccs_count_values(const jmtx_matrix_ccs* mtx, float v, uint32_t* p_count);

/**
 * Counts the number of times a specific row index occurs in the matrix
 * @param mtx pointer to the memory where the matrix is stored
 * @param v row index which to search for
 * @param p_count pointer which receives the number of times the value appeared in the matrix
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result jmtxs_matrix_ccs_count_indices(const jmtx_matrix_ccs* mtx, uint32_t v, uint32_t* p_count);

/**
 * Applies a unary function on the sparse matrix, only on its stored entries, which can be modified. If the user given
 * user function returns a non-zero value, that value is returned from the function and the iteration is stopped
 * @param mtx pointer to the memory where the matrix is stored
 * @param unary_fn user given function that is to be applied on the values
 * @param param optional parameter, which is passed to the unary_fn when called
 * @return JMTX_RESULT_SUCCESS if user function never returned non-zero, JMTX_RESULT_UNARY_RETURN as soon as the user
 * function returns non-zero
 */
jmtx_result jmtxs_matrix_ccs_apply_unary_fn(const jmtx_matrix_ccs* mtx, int (*unary_fn)(uint32_t i, uint32_t j, float* p_element, void* param), void* param);

/**
 * Removes entries exactly equal to zero. If the element is indeed zero, it is compared to (float)0.0
 * @param mtx pointer to the memory where the matrix is stored
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result jmtxs_matrix_ccs_remove_zeros(jmtx_matrix_ccs* mtx);

/**
 * Removes values which have absolute value less than specified value
 * @param mtx pointer to the memory where the matrix is stored
 * @param v value to which to compare it to
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result jmtxs_matrix_ccs_remove_bellow(jmtx_matrix_ccs* mtx, float v);

/**
 * Zeros all values within a matrix, but does not remove them in case they need to be reused
 * @param mtx matrix to zero
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result jmtxs_matrix_zero_all_entries(const jmtx_matrix_ccs* mtx);

/**
 * Similar to jmtx_matrix_zero_all_entries, but slower, since it can not use memset. On the other hand, it allows for
 * the value to be other than 0
 * @param mtx matrix to set
 * @param x value to which to set to
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result jmtxs_matrix_ccs_set_all_entries(jmtx_matrix_ccs* mtx, float x);

/**
 * Returns the number of values in the column of the matrix
 * @param mtx pointer to the memory where the matrix is stored
 * @param row column index of the matrix to look at
 * @param p_n pointer to integer which receives the number of values in the column
 * @return zero if successful
 */
jmtx_result jmtxs_matrix_ccs_elements_in_row(const jmtx_matrix_ccs* mtx, uint32_t row, uint32_t* p_n);

/**
 * Returns the values of values in the matrix, along with what column of the matrix they were located in
 * @param mtx pointer to the memory where the matrix is stored
 * @param row column index of the matrix to look at
 * @param n number of values in the column to be extracted
 * @param p_values a buffer of at least n values which receives the values of the column
 * @param p_columns a buffer of at least n values which receives the row indices of the row
 * @param p_count pointer which receives the number of entries that were extracted from the row (may be less than are
 * really in the row if n was too small)
 * @return zero if successful
 */
jmtx_result jmtxs_matrix_ccs_get_row(
const jmtx_matrix_ccs* mtx, uint32_t row, uint32_t n, float* p_values, uint32_t* p_count, uint32_t* p_columns);

/**
 * Creates a transpose of a matrix
 * @param mtx pointer to the memory where the input matrix is stored
 * @param p_out address where the pointer to the output matrix will be returned
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result jmtxs_matrix_ccs_transpose(const jmtx_matrix_ccs* mtx, jmtx_matrix_ccs** p_out,
                                      const jmtx_allocator_callbacks* allocator_callbacks);

/**
 * Creates a copy of the matrix
 * @param mtx pointer to the memory where the input matrix is stored
 * @param p_out address where the pointer to the output matrix will be returned
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result jmtxs_matrix_ccs_copy(const jmtx_matrix_ccs* restrict mtx, jmtx_matrix_ccs** p_out, const jmtx_allocator_callbacks* allocator_callbacks);
#endif //JMTX_SPARSE_COLUMN_COMPRESSED_SAFE_H
