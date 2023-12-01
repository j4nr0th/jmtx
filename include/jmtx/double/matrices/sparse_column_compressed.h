// Automatically generated from include/jmtx/float/matrices/sparse_column_compressed.h on Fri Dec  1 06:43:05 2023
//
// Created by jan on 15.6.2022.
//
/**
 * Functions declared here perform minimum checking of the parameters and assume that values that were passed to them
 * were valid (matrix have proper dimensions, indices are within required bounds, etc.). "Safe" versions of these functions,
 * which do perform parameter validation are in the "sparse_row_compressed_safe.h" header.
 */

#ifndef JMTXD_SPARSE_COLUMN_COMPRESSED_H
#define JMTXD_SPARSE_COLUMN_COMPRESSED_H
#ifndef JMTXD_MATRIX_BASE_H
    #include "../../matrix_base.h"
#endif



typedef struct jmtxd_matrix_ccs_struct jmtxd_matrix_ccs;

/**
 * Initializes a new Compressed Column Sparse matrix
 * @param p_mtx address that receives the pointer to the matrix
 * @param cols number of columns of the sparse matrix
 * @param rows number of rows of the sparse matrix
 * @param reserved_entries how many entries should the space be reserved for in the matrix initially
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxd_matrix_ccs_new(
        jmtxd_matrix_ccs** p_mtx, uint32_t cols, uint32_t rows, uint32_t reserved_entries,
        const jmtx_allocator_callbacks* allocator_callbacks);

/**
 * Cleans up the ccs matrix and frees all of its memory
 * @param mtx pointer to memory where the matrix is stored
 */
void jmtxd_matrix_ccs_destroy(jmtxd_matrix_ccs* mtx);

/**
 * Frees up memory which the matrix is not currently using, which is was allocated in advance
 * @param mtx pointer to the memory where the matrix is stored
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxd_matrix_ccs_shrink(jmtxd_matrix_ccs* mtx);

/**
 * Sets the column of the matrix. This is the most efficient way to build the matrix, as building it this way causes
 * minimum amount of memory allocation.
 * @param mtx pointer to the memory where the matrix is stored
 * @param col index of the column to set
 * @param n how many entries are in the column
 * @param indices row indices of the entries, which has to be sorted from lowest to highest
 * @param values values of non-zero entries
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxd_matrix_ccs_set_col(jmtxd_matrix_ccs* mtx, uint32_t col, uint32_t n, const uint32_t* indices, const double* values);

/**
 * Version of jmtxd_matrix_ccs_set_col which does not touch the count of entries after the current column. This is useful
 * when building a new matrix, as it avoids unnecessary setting and resetting of these entries. Must be called for each
 * column in order to ensure that the matrix is properly built. Makes no checks on the input parameters
 * @param mtx pointer to the memory where the matrix is stored
 * @param col index of the column to set
 * @param n how many entries are in the column
 * @param indices column indices of the values, which has to be sorted from lowest to highest
 * @param values values of non-zero values
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxd_matrix_ccs_build_col(jmtxd_matrix_ccs* mtx, uint32_t col, uint32_t n, const uint32_t* indices, const double* values);

/**
 * Returns the pointers to arrays of column indices and element values for that column
 * @param mtx pointer to the memory where the matrix is stored
 * @param col index of the row to get
 * @param p_indices pointer to row indices of values
 * @param p_elements pointer to values of values
 * @return number of elements in the column, which is the number of valid elements in arrays given to p_indices and
 * p_elements
 */
uint32_t jmtxd_matrix_ccs_get_col(const jmtxd_matrix_ccs* mtx, uint32_t col, uint32_t** p_indices, double** p_elements);

/**
 * Multiplies a dense row vector x by the sparse matrix and stores the result at y
 * @param mtx pointer to the memory where the matrix is stored
 * @param x pointer to vector to be multiplied
 * @param y pointer to vector where the result of multiplication is to be stored
 */
void jmtxd_matrix_ccs_vector_multiply(const jmtxd_matrix_ccs* mtx, const double* restrict x, double* restrict y);

/**
 * Sets a single entry in the matrix. This is about as fast as setting the entire column of the matrix at once, if the
 * element was previously zero
 * @param mtx pointer to the memory where the matrix is stored
 * @param i row index
 * @param j column index
 * @param value value to which the value is set
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxd_matrix_ccs_set_entry(jmtxd_matrix_ccs* mtx, uint32_t i, uint32_t j, double value);

/**
 * Returns a single entry from the matrix.
 * @param mtx pointer to the memory where the matrix is stored
 * @param i row index
 * @param j column index
 * @return value of the entry (0 if the entry was not manually set to anything else)
 */
double jmtxd_matrix_ccs_get_entry(const jmtxd_matrix_ccs* mtx, uint32_t i, uint32_t j);

/**
 * Counts the number of times a specific value occurs in the matrix
 * @param mtx pointer to the memory where the matrix is stored
 * @param v value which to search for
 * @return number of times the value appeared in the matrix
 */
uint32_t jmtxd_matrix_ccs_count_values(const jmtxd_matrix_ccs* mtx, double v);

/**
 * Counts the number of times a specific row index occurs in the matrix
 * @param mtx pointer to the memory where the matrix is stored
 * @param v row index which to search for
 * @return number of times the row index appeared in the matrix
 */
uint32_t jmtxd_matrix_ccs_count_indices(const jmtxd_matrix_ccs* mtx, uint32_t v);

/**
 * Applies a unary function on the sparse matrix, only on its stored entries, which can be modified. If the user given
 * user function returns a non-zero value, that value is returned from the function and the iteration is stopped
 * @param mtx pointer to the memory where the matrix is stored
 * @param unary_fn user given function that is to be applied on the values
 * @param param optional parameter, which is passed to the unary_fn when called
 * @return JMTX_RESULT_SUCCESS if user function never returned non-zero, JMTX_RESULT_UNARY_RETURN as soon as the user
 * function returns non-zero
 */
jmtx_result jmtxd_matrix_ccs_apply_unary_fn(const jmtxd_matrix_ccs* mtx, int (*unary_fn)(uint32_t i, uint32_t j, double* p_element, void* param), void* param);

/**
 * Removes entries exactly equal to zero. If the element is indeed zero, it is compared to (double)0.0
 * @param mtx pointer to the memory where the matrix is stored
 */
void jmtxd_matrix_ccs_remove_zeros(jmtxd_matrix_ccs* mtx);

/**
 * Removes entries which have absolute value less than specified value
 * @param mtx pointer to the memory where the matrix is stored
 * @param v value to which to compare it to
 */
void jmtxd_matrix_ccs_remove_bellow(jmtxd_matrix_ccs* mtx, double v);

/**
 * Zeros all entries within a matrix, but does not remove them in case they need to be reused
 * @param mtx matrix to zero
 */
void jmtxd_matrix_ccs_zero_all_elements(const jmtxd_matrix_ccs* mtx);

/**
 * Similar to jmtxd_matrix_crs_zero_all_entries, but slower, since it can not use memset. On the other hand, it allows for
 * the value to be other than 0
 * @param mtx matrix to set
 * @param x value to which to set all entries to
 */
void jmtxd_matrix_ccs_set_all_elements(const jmtxd_matrix_ccs* mtx, double x);

/**
 * Returns the number of entries in the row of the matrix
 * @param mtx pointer to the memory where the matrix is stored
 * @param row column index of the matrix to look at
 * @return number of entries in the row
 */
uint32_t jmtxd_matrix_ccs_elements_in_row(const jmtxd_matrix_ccs* mtx, uint32_t row);

/**
 * Returns the values of entries in the matrix, along with what column of the matrix they were located in
 * @param mtx pointer to the memory where the matrix is stored
 * @param row row index of the matrix to look at
 * @param n number of values in the row to be extracted at most
 * @param p_values a buffer of at least n values which receives the values of the row
 * @param p_columns a buffer of at least n values which receives the row indices of the row
 * @return number of entries that were extracted from the row (may be less than are really in the row if n was too
 * small)
 */
uint32_t jmtxd_matrix_ccs_get_row(
 const jmtxd_matrix_ccs* mtx, uint32_t row, uint32_t n, double* p_values, uint32_t* p_columns);

/**
 * Creates a transpose of a matrix
 * @param mtx pointer to the memory where the input matrix is stored
 * @param p_out address where the pointer to the output matrix will be returned
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxd_matrix_ccs_transpose(const jmtxd_matrix_ccs* mtx, jmtxd_matrix_ccs** p_out,
                                      const jmtx_allocator_callbacks* allocator_callbacks);

/**
 * Creates a copy of the matrix
 * @param mtx pointer to the memory where the input matrix is stored
 * @param p_out address where the pointer to the output matrix will be returned
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxd_matrix_ccs_copy(const jmtxd_matrix_ccs* restrict mtx, jmtxd_matrix_ccs** p_out, const jmtx_allocator_callbacks* allocator_callbacks);
#endif //JMTXD_SPARSE_COLUMN_COMPRESSED_H
