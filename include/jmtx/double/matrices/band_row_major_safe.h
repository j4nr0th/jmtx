// Automatically generated from include/jmtx/float/matrices/band_row_major_safe.h on Fri Dec  1 06:43:05 2023
//
// Created by jan on 13.6.2022.
//
/**
 * Functions declared here perform more checking of the parameters. Faster "unsafe" versions of these functions,
 * which do perform parameter validation are in the "band_row_major.h" header.
 */

#ifndef JMTXD_BAND_ROW_MAJOR_SAFE_H
#define JMTXD_BAND_ROW_MAJOR_SAFE_H
#ifndef JMTXD_BAND_ROW_MAJOR_H
    #include "band_row_major.h"
#endif
/**
 * @paragraph
 * Band Row-Major matrix (BRM) is a matrix which has constant upper bandwidth (ubw) and lower bandwidths (lbw): constant
 * number of entries above and bellow the diagonal. The rows are stored contiguously. The memory required is bounded by
 * lbw + 1 + ubw. Memory access is performed in constant time for both rows and columns, with rows being slightly
 * faster and doesn't need a separate function for setting it, since a pointer to a row allows modification.
 *
 * @paragraph
 * Main advantage of BRM matrices is the fact that for a BRM with bandwidths ubw and lbw, its exact LU decomposition
 * results in the matrix L having bandwidths 0 and lbw, and matrix U having bandwidths ubw and 0. This means that for
 * a chosen matrix number of elements that need to be computed for full LU decomposition becomes (lbw + 1 + ubw) * N,
 * with N being the size of the matrix. This means that using full LU decomposition can be viable, given that the
 * bandwidth is low (most PDEs on a 1D domain).
 */


/**
 * Initializes a new Band Row Major matrix
 * @param p_mtx address that receives the pointer to the matrix
 * @param rows number of rows of the matrix
 * @param cols number of columns of the matrix
 * @param ubw upper bandwidth
 * @param lbw lower bandwidth
 * @param set_value pointer to value with which to initialize all entries. If NULL, then matrix is left uninitialized
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxds_matrix_brm_new(
        jmtxd_matrix_brm** p_mtx, uint32_t cols, uint32_t rows, uint32_t ubw, uint32_t lbw, const double* set_value,
        const jmtx_allocator_callbacks* allocator_callbacks);

/**
 * Cleans up the BRM matrix and frees all of its memory
 * @param mtx pointer to memory where the matrix is stored
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result jmtxds_matrix_brm_destroy(jmtxd_matrix_brm* mtx);

/**
 * Sets the row of the matrix. More efficient than setting it element by element
 * @param mtx pointer to the memory where the matrix is stored
 * @param row index of the row to set
 * @param values values of non-zero entries
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result jmtxds_matrix_brm_set_row(const jmtxd_matrix_brm* mtx, uint32_t row, double values[]);

/**
 * Sets the column of the matrix. More efficient than setting it element by element
 * @param mtx pointer to the memory where the matrix is stored
 * @param row index of the column to set
 * @param values values of entries
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result jmtxds_matrix_brm_set_col(const jmtxd_matrix_brm* mtx, uint32_t col, const double* values);

/**
 * Returns the pointers to arrays of column indices and element values for that row
 * @param mtx pointer to the memory where the matrix is stored
 * @param row index of the row to get
 * @param p_elements pointer to array of values
 * @param n pointer that receives the number of elements in the row, which is the number of valid elements in arrays
 * given to p_indices and p_elements
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result jmtxds_matrix_brm_get_row(const jmtxd_matrix_brm* mtx, uint_fast32_t* n, uint32_t row, double* p_elements[1]);

/**
 * Multiplies a dense column vector x by the sparse matrix and stores the result at y
 * @param mtx pointer to the memory where the matrix is stored
 * @param x pointer to vector to be multiplied
 * @param y pointer to vector where the result of multiplication is to be stored
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result jmtxds_matrix_brm_vector_multiply(const jmtxd_matrix_brm* mtx, const double* restrict x, double* restrict y);

/**
 * Sets a single entry in the matrix. This is about as fast as setting the entire row of the matrix at once, if the
 * element was previously zero
 * @param mtx pointer to the memory where the matrix is stored
 * @param i row index
 * @param j column index
 * @param value value to which the value is set
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result jmtxds_matrix_brm_set_entry(const jmtxd_matrix_brm* mtx, uint32_t i, uint32_t j, double value);

/**
 * Returns a single entry from the matrix.
 * @param mtx pointer to the memory where the matrix is stored
 * @param i row index
 * @param j column index
 * @param p_value pointer that receives the value of the entry (0 if the entry was not manually set to anything else)
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result jmtxds_matrix_brm_get_entry(const jmtxd_matrix_brm* mtx, uint32_t i, uint32_t j, double* p_value);

/**
 * Adds a value to an entry in the matrix when it exists or sets it to that value if it does not. This is about as
 * fast as setting the entire row of the matrix at once, if the entry was non-existent.
 * @param mtx pointer to the memory where the matrix is stored
 * @param i row index
 * @param j column index
 * @param value value to which the value is to be added
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result jmtxds_matrix_brm_add_to_entry(const jmtxd_matrix_brm* mtx, uint32_t i, uint32_t j, double value);

/**
 * Counts the number of times a specific value occurs in the matrix
 * @param mtx pointer to the memory where the matrix is stored
 * @param v value which to search for
 * @param p_count pointer which receives the number of times the value appeared in the matrix
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result jmtxds_matrix_brm_count_values(const jmtxd_matrix_brm* mtx, double v, uint32_t* p_count);

/**
 * Zeros all entries within a matrix, but does not remove them in case they need to be reused
 * @param mtx matrix to zero
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result jmtxds_matrix_brm_zero_all_entries(const jmtxd_matrix_brm* mtx);

/**
 * Similar to jmtxd_matrix_brm_zero_all_entries, but slower, since it can not use memset. On the other hand, it allows for
 * the value to be other than 0
 * @param mtx matrix to set
 * @param x value to which to set all entries to
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result jmtxds_matrix_brm_set_all_entries(const jmtxd_matrix_brm* mtx, double x);

/**
 * Returns the number of entries in the column of the matrix
 * @param mtx pointer to the memory where the matrix is stored
 * @param col column index of the matrix to look at
 * @param p_n pointer which receives the number of entries in the column
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result jmtxds_matrix_brm_entries_in_col(const jmtxd_matrix_brm* mtx, uint32_t col, uint32_t* p_n);

/**
 * Returns the values of entries in the matrix, along with what row of the matrix they were located in
 * @param mtx pointer to the memory where the matrix is stored
 * @param col column index of the matrix to look at
 * @param p_values a buffer of at least n values which receives the values of the column
 * @param p_count the pointer which receives number of entries that were extracted from the column
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result jmtxds_matrix_brm_get_col(
        const jmtxd_matrix_brm* mtx, uint32_t col, uint32_t* p_count, double p_values[]);

/**
 * Creates a transpose of a matrix
 * @param mtx pointer to the memory where the input matrix is stored
 * @param p_out address where the pointer to the output matrix will be returned
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxds_matrix_brm_transpose(
        const jmtxd_matrix_brm* mtx, jmtxd_matrix_brm** p_out, const jmtx_allocator_callbacks* allocator_callbacks);

/**
 * Creates a copy of the matrix
 * @param mtx pointer to the memory where the input matrix is stored
 * @param p_out address where the pointer to the output matrix will be returned
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxds_matrix_brm_copy(const jmtxd_matrix_brm* mtx, jmtxd_matrix_brm** p_out, const jmtx_allocator_callbacks* allocator_callbacks);


/**
 * Computes one entry of Ax. This function only computes the i-th entry to make it possible to compute it in parallel.
 * @param mtx pointer to the memory where the matrix A is stored.
 * @param x pointer to the memory where the vector x is stored
 * @param i what entry of the residual to compute
 * @param p_r pointer which receives the result of the multiplication
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result jmtxds_matrix_brm_vector_multiply_row(const jmtxd_matrix_brm* mtx, const double* restrict x, uint32_t i, double* restrict p_r);


#endif //JMTXD_BAND_ROW_MAJOR_SAFE_H
