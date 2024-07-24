// Automatically generated from include/jmtx/float/matrices/sparse_row_compressed_safe.h on Fri Dec  1 06:43:05 2023
//
// Created by jan on 13.6.2022.
//
/**
 * Functions declared here perform more checking of the parameters. Faster "unsafe" versions of these functions,
 * which do perform parameter validation are in the "sparse_row_compressed.h" header.
 *
 */

#ifndef JMTXD_SPARSE_ROW_COMPRESSED_SAFE_H
#define JMTXD_SPARSE_ROW_COMPRESSED_SAFE_H
#ifndef JMTXD_SPARSE_ROW_COMPRESSED_H
    #include "sparse_row_compressed.h"
#endif

/**
 * Initializes a new Compressed Row Sparse matrix
 * @param p_mtx address that receives the pointer to the matrix
 * @param rows number of rows of the sparse matrix
 * @param cols number of columns of the sparse matrix
 * @param reserved_entries how many entries should the space be reserved for in the matrix initially
 * @return JMTX_RESULT_SUCCESS if successful
 */
JMTX_NODISCARD_FUNCTION
jmtx_result jmtxds_matrix_crs_new(
 jmtxd_matrix_crs** p_mtx, uint32_t rows, uint32_t cols, uint32_t reserved_entries,
 const jmtx_allocator_callbacks* allocator_callbacks);

/**
 * Cleans up the crs matrix and frees all of its memory
 * @param mtx pointer to memory where the matrix is stored
 * @return JMTX_RESULT_SUCCESS if successful
 */
JMTX_NODISCARD_FUNCTION
jmtx_result jmtxds_matrix_crs_destroy(jmtxd_matrix_crs* mtx);

/**
 * Frees up memory which the matrix is not currently using, which is was allocated in advance
 * @param mtx pointer to the memory where the matrix is stored
 * @return JMTX_RESULT_SUCCESS if successful
 */
JMTX_NODISCARD_FUNCTION
jmtx_result jmtxds_matrix_crs_shrink(jmtxd_matrix_crs* mtx);

/**
 * Sets the row of the matrix. More efficient than setting it element by element
 * @param mtx pointer to the memory where the matrix is stored
 * @param row index of the row to set
 * @param n how many entries are in the row
 * @param indices column indices of the entries, which has to be sorted from lowest to highest
 * @param values values of non-zero entries
 * @return JMTX_RESULT_SUCCESS if successful
 */
JMTX_NODISCARD_FUNCTION
jmtx_result jmtxds_matrix_crs_set_row(jmtxd_matrix_crs* mtx, uint32_t row, uint32_t n, const uint32_t indices[JMTX_ARRAY_ATTRIB(static n)], const double values[JMTX_ARRAY_ATTRIB(static n)]);

/**
 * Returns the pointers to arrays of column indices and element values for that row
 * @param mtx pointer to the memory where the matrix is stored
 * @param row index of the row to get
 * @param n pointer which receives the number of values in the row
 * @param p_indices pointer to column indices of values
 * @param p_elements pointer to values of values
 * @return JMTX_RESULT_SUCCESS if successful
 */
JMTX_NODISCARD_FUNCTION
jmtx_result jmtxds_matrix_crs_get_row(const jmtxd_matrix_crs* mtx, uint32_t row, uint32_t* n, uint32_t* p_indices[1], double* p_elements[1]);

/**
 * Multiplies a dense column vector x by the sparse matrix and stores the result at y
 * @param mtx pointer to the memory where the matrix is stored
 * @param x pointer to vector to be multiplied
 * @param y pointer to vector where the result of multiplication is to be stored
 * @return JMTX_RESULT_SUCCESS if successful
 */
JMTX_NODISCARD_FUNCTION
jmtx_result jmtxds_matrix_crs_vector_multiply(const jmtxd_matrix_crs* mtx, const double* restrict x, double* restrict y);

/**
 * Sets a single entry in the matrix. This is about as fast as setting the entire row of the matrix at once, if the
 * element was previously zero
 * @param mtx pointer to the memory where the matrix is stored
 * @param i row index
 * @param j column index
 * @param value value to which the value is set
 * @return JMTX_RESULT_SUCCESS if successful
 */
JMTX_NODISCARD_FUNCTION
jmtx_result jmtxds_matrix_crs_set_entry(jmtxd_matrix_crs* mtx, uint32_t i, uint32_t j, double value);

/**
 * Returns a single entry from the matrix.
 * @param mtx pointer to the memory where the matrix is stored
 * @param i row index
 * @param j column index
 * @param p_value pointer which receives the value
 * @return JMTX_RESULT_SUCCESS if successful
 */
JMTX_NODISCARD_FUNCTION
jmtx_result jmtxds_matrix_crs_get_entry(const jmtxd_matrix_crs* mtx, uint32_t i, uint32_t j, double* p_value);

/**
 * Adds a value to an entry in the matrix when it exists or sets it to that value if it does not. This is about as
 * fast as setting the entire row of the matrix at once, if the entry was non-existent.
 * @param mtx pointer to the memory where the matrix is stored
 * @param i row index
 * @param j column index
 * @param value value to which the value is to be added
 * @return JMTX_RESULT_SUCCESS if successful
 */
JMTX_NODISCARD_FUNCTION
jmtx_result jmtxds_matrix_crs_add_to_entry(jmtxd_matrix_crs* mtx, uint32_t i, uint32_t j, double value);

/**
 * Counts the number of times a specific value occurs in the matrix
 * @param mtx pointer to the memory where the matrix is stored
 * @param v value which to search for
 * @param p_count pointer that receives the number of times the value appeared in the matrix
 * @return JMTX_RESULT_SUCCESS if successful
 */
JMTX_NODISCARD_FUNCTION
jmtx_result jmtxds_matrix_crs_count_values(const jmtxd_matrix_crs* mtx, double v, uint32_t* p_count);

/**
 * Counts the number of times a specific column index occurs in the matrix
 * @param mtx pointer to the memory where the matrix is stored
 * @param v column index which to search for
 * @param p_count pointer that receives the number of times the column index appeared in the matrix
 * @return JMTX_RESULT_SUCCESS if successful
 */
JMTX_NODISCARD_FUNCTION
jmtx_result jmtxds_matrix_crs_count_indices(const jmtxd_matrix_crs* mtx, uint32_t v, uint32_t* p_count);

/**
 * Applies a unary function on the sparse matrix, only on its stored entries, which can be modified. If the user given
 * user function returns a non-zero value, that value is returned from the function and the iteration is stopped
 * @param mtx pointer to the memory where the matrix is stored
 * @param unary_fn user given function that is to be applied on the values
 * @param param optional parameter, which is passed to the unary_fn when called
 * @return JMTX_RESULT_SUCCESS if successful, if the user function returns non-zero, that value is returned instead
 */
JMTX_NODISCARD_FUNCTION
jmtx_result jmtxds_matrix_crs_apply_unary_fn(const jmtxd_matrix_crs* mtx, int (*unary_fn)(uint32_t i, uint32_t j, double* p_value, void* param), void* param);

/**
 * Removes entries exactly equal to zero. If the element is indeed zero, it is compared to (scalar_t)0.0
 * @param mtx pointer to the memory where the matrix is stored
 * @return JMTX_RESULT_SUCCESS if successful
 */
JMTX_NODISCARD_FUNCTION
jmtx_result jmtxds_matrix_crs_remove_zeros(jmtxd_matrix_crs* mtx);

/**
 * Removes entries which have absolute value less than specified value
 * @param mtx pointer to the memory where the matrix is stored
 * @param v value to which to compare it to
 * @return JMTX_RESULT_SUCCESS if successful
 */
JMTX_NODISCARD_FUNCTION
jmtx_result jmtxds_matrix_crs_remove_bellow_magnitude(jmtxd_matrix_crs* mtx, double v);

/**
 * Zeros all entries within a matrix, but does not remove them in case they need to be reused
 * @param mtx matrix to zero
 * @return JMTX_RESULT_SUCCESS if successful
 */
JMTX_NODISCARD_FUNCTION
jmtx_result jmtxds_matrix_crs_zero_all_entries(const jmtxd_matrix_crs* mtx);

/**
 * Similar to jmtxd_matrix_crs_zero_all_entries, but slower, since it can not use memset. On the other hand, it allows for
 * the value to be other than 0
 * @param mtx matrix to set
 * @param x value to which to set all entries to
 * @return JMTX_RESULT_SUCCESS if successful
 */
JMTX_NODISCARD_FUNCTION
jmtx_result jmtxds_matrix_crs_set_all_entries(jmtxd_matrix_crs* mtx, double x);

/**
 * Keeps the memory for the matrix, but sets the entry count to 0, so that matrix can be rebuilt.
 * @param mtx matrix to clear
 * @return JMTX_RESULT_SUCCESS if successful
 */
JMTX_NODISCARD_FUNCTION
jmtx_result jmtxds_matrix_crs_clear(jmtxd_matrix_crs* mtx);

/**
 * Returns the number of entries in the column of the matrix
 * @param mtx pointer to the memory where the matrix is stored
 * @param col column index of the matrix to look at
 * @param p_n pointer to integer which receives the number of entries in the column
 * @return JMTX_RESULT_SUCCESS if successful
 */
JMTX_NODISCARD_FUNCTION
jmtx_result jmtxds_matrix_crs_entries_in_col(const jmtxd_matrix_crs* mtx, uint32_t col, uint32_t* p_n);

/**
 * Returns the values of entries in the matrix, along with what row of the matrix they were located in
 * @param mtx pointer to the memory where the matrix is stored
 * @param col column index of the matrix to look at
 * @param n number of values in the column to be extracted
 * @param p_values a buffer of at least n values which receives the values of the column
 * @param p_rows a buffer of at least n values which receives the row indices of the column
 * @return JMTX_RESULT_SUCCESS if successful
 */
JMTX_NODISCARD_FUNCTION
jmtx_result jmtxds_matrix_crs_get_col(
        const jmtxd_matrix_crs* mtx, uint32_t col, uint32_t n, uint32_t* p_count, double p_values[JMTX_ARRAY_ATTRIB(n)], uint32_t p_rows[JMTX_ARRAY_ATTRIB(n)]);

/**
 * Creates a transpose of a matrix
 * @param mtx pointer to the memory where the input matrix is stored
 * @param p_out address where the pointer to the output matrix will be returned
 * @return JMTX_RESULT_SUCCESS if successful
 */
JMTX_NODISCARD_FUNCTION
jmtx_result jmtxds_matrix_crs_transpose(
        const jmtxd_matrix_crs* mtx, jmtxd_matrix_crs** p_out, const jmtx_allocator_callbacks* allocator_callbacks);

/**
 * Creates a copy of the matrix
 * @param mtx pointer to the memory where the input matrix is stored
 * @param p_out address where the pointer to the output matrix will be returned
 * @return JMTX_RESULT_SUCCESS if successful
 */
JMTX_NODISCARD_FUNCTION
jmtx_result jmtxds_matrix_crs_copy(const jmtxd_matrix_crs* mtx, jmtxd_matrix_crs** p_out, const jmtx_allocator_callbacks* allocator_callbacks);

/**
 * Computes one entry of Ax. This function only computes the i-th entry to make it possible to compute it in parallel.
 * @param mtx pointer to the memory where the matrix A is stored.
 * @param x pointer to the memory where the vector x is stored
 * @param i what entry of the residual to compute
 * @param p_r pointer to the memory which receives the result
 * @return JMTX_RESULT_SUCCESS if successful
 */
JMTX_NODISCARD_FUNCTION
jmtx_result jmtxds_matrix_crs_vector_multiply_row(const jmtxd_matrix_crs* mtx, const double* restrict x, uint32_t i, double* restrict p_r);

/**
 * Removes a row from a CRS matrix, reducing it's row count by 1
 * @param mtx matrix to remove the row from
 * @param row row to remove
 * @return JMTX_RESULT_SUCCESS if successful
 */
JMTX_NODISCARD_FUNCTION
jmtx_result jmtxds_matrix_crs_remove_row(jmtxd_matrix_crs* mtx, uint32_t row);

/**
 * Removes a column from a CRS matrix, reducing it's column count by 1. Slower than removing a row, so remove rows before
 * columns if both have to be done.
 * @param mtx matrix to remove the column from
 * @param col column to remove
 * @return JMTX_RESULT_SUCCESS if successful
 */
JMTX_NODISCARD_FUNCTION
jmtx_result jmtxds_matrix_crs_remove_column(jmtxd_matrix_crs* mtx, uint32_t col);

/**
 * Combines k matrices of size N_i x M into a single Sum(N_i) x M matrix by vertically stacking them
 * @param output receives pointer to the resulting matrix
 * @param allocators memory allocators to use for the output matrix
 * @param k number of the matrices
 * @param matrix_list array of matrices to be joined together. Must have the same number of columns.
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_DIMS_MISMATCH if the number of columns is not the same for all
 * input matrices, return value of jmtxd_matrix_crs_new if that fails
 */
JMTX_NODISCARD_FUNCTION
jmtx_result jmtxds_matrix_crs_join_vertically(jmtxd_matrix_crs** output, const jmtx_allocator_callbacks* allocators, unsigned k, const jmtxd_matrix_crs* matrix_list[JMTX_ARRAY_ATTRIB(static k)]);

#endif //JMTXD_SPARSE_ROW_COMPRESSED_SAFE_H
