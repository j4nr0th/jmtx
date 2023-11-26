//
// Created by jan on 13.6.2022.
//

#ifndef JMTX_BAND_ROW_MAJOR_H
#define JMTX_BAND_ROW_MAJOR_H
#ifndef JMTX_MATRIX_BASE_H
    #include "matrix_base.h"
#endif
/**
 * Functions declared here perform minimum checking of the parameters and assume that values that were passed to them
 * were valid (matrix have proper dimensions, indices are within required bounds, etc.). "Safe" versions of these functions,
 * which do perform parameter validation are in the "band_row_major_safe.h" header.
 */




typedef struct jmtx_matrix_brm_struct jmtx_matrix_brm;


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
jmtx_result jmtx_matrix_brm_new(
        jmtx_matrix_brm** p_mtx, uint32_t cols, uint32_t rows, uint32_t ubw, uint32_t lbw, const float* set_value,
        const jmtx_allocator_callbacks* allocator_callbacks);

/**
 * Cleans up the BRM matrix and frees all of its memory
 * @param mtx pointer to memory where the matrix is stored
 */
void jmtx_matrix_brm_destroy(jmtx_matrix_brm* mtx);

/**
 * Sets the row of the matrix. More efficient than setting it element by element
 * @param mtx pointer to the memory where the matrix is stored
 * @param row index of the row to set
 * @param values values of entries
 */
void jmtx_matrix_brm_set_row(const jmtx_matrix_brm* mtx, uint32_t row, float values[]);

/**
 * Sets the column of the matrix. More efficient than setting it element by element
 * @param mtx pointer to the memory where the matrix is stored
 * @param row index of the column to set
 * @param values values of entries
 */
void jmtx_matrix_brm_set_col(const jmtx_matrix_brm* mtx, uint32_t col, const float* values);

/**
 * Returns the pointers to arrays of column indices and element values for that row
 * @param mtx pointer to the memory where the matrix is stored
 * @param row index of the row to get
 * @param p_elements pointer to array of values
 * @return number of elements in the row, which is the number of valid elements in arrays given to p_indices and
 * p_elements
 */
uint_fast32_t jmtx_matrix_brm_get_row(const jmtx_matrix_brm* mtx, uint32_t row, float* p_elements[1]);

uint_fast32_t jmtx_matrix_brm_first_pos_in_row(const jmtx_matrix_brm* mtx, uint32_t row);

uint_fast32_t jmtx_matrix_brm_last_pos_in_row(const jmtx_matrix_brm* mtx, uint32_t row);

uint_fast32_t jmtx_matrix_brm_length_of_row(const jmtx_matrix_brm* mtx, uint32_t row);

uint_fast32_t jmtx_matrix_brm_first_pos_in_col(const jmtx_matrix_brm* mtx, uint32_t col);

uint_fast32_t jmtx_matrix_brm_last_pos_in_col(const jmtx_matrix_brm* mtx, uint32_t col);

uint_fast32_t jmtx_matrix_brm_length_of_col(const jmtx_matrix_brm* mtx, uint32_t col);


/**
 * Multiplies a dense column vector x by the sparse matrix and stores the result at y
 * @param mtx pointer to the memory where the matrix is stored
 * @param x pointer to vector to be multiplied
 * @param y pointer to vector where the result of multiplication is to be stored
 */
void jmtx_matrix_brm_vector_multiply(const jmtx_matrix_brm* mtx, const float* restrict x, float* restrict y);

/**
 * Sets a single entry in the matrix. This is about as fast as setting the entire row of the matrix at once, if the
 * element was previously zero
 * @param mtx pointer to the memory where the matrix is stored
 * @param i row index
 * @param j column index
 * @param value value to which the value is set
 */
void jmtx_matrix_brm_set_entry(const jmtx_matrix_brm* mtx, uint32_t i, uint32_t j, float value);

/**
 * Returns a single entry from the matrix.
 * @param mtx pointer to the memory where the matrix is stored
 * @param i row index
 * @param j column index
 * @return value of the entry (0 if the entry was not manually set to anything else)
 */
float jmtx_matrix_brm_get_entry(const jmtx_matrix_brm* mtx, uint32_t i, uint32_t j);

/**
 * Adds a value to an entry in the matrix when it exists or sets it to that value if it does not. This is about as
 * fast as setting the entire row of the matrix at once, if the entry was non-existent.
 * @param mtx pointer to the memory where the matrix is stored
 * @param i row index
 * @param j column index
 * @param value value to which the value is to be added
 */
void jmtx_matrix_brm_add_to_entry(const jmtx_matrix_brm* mtx, uint32_t i, uint32_t j, float value);

/**
 * Counts the number of times a specific value occurs in the matrix
 * @param mtx pointer to the memory where the matrix is stored
 * @param v value which to search for
 * @return number of times the value appeared in the matrix
 */
uint32_t jmtx_matrix_brm_count_values(const jmtx_matrix_brm* mtx, float v);

/**
 * Zeros all entries within a matrix, but does not remove them in case they need to be reused
 * @param mtx matrix to zero
 */
void jmtx_matrix_brm_zero_all_entries(const jmtx_matrix_brm* mtx);

/**
 * Similar to jmtx_matrix_brm_zero_all_entries, but slower, since it can not use memset. On the other hand, it allows for
 * the value to be other than 0
 * @param mtx matrix to set
 * @param x value to which to set all entries to
 */
void jmtx_matrix_brm_set_all_entries(const jmtx_matrix_brm* mtx, float x);

/**
 * Returns the number of entries in the column of the matrix
 * @param mtx pointer to the memory where the matrix is stored
 * @param col column index of the matrix to look at
 * @return number of entries in the column
 */
uint32_t jmtx_matrix_brm_entries_in_col(const jmtx_matrix_brm* mtx, uint32_t col);

/**
 * Returns the values of entries in the matrix, along with what row of the matrix they were located in
 * @param mtx pointer to the memory where the matrix is stored
 * @param col column index of the matrix to look at
 * @param p_values a buffer of at least n values which receives the values of the column
 * @return number of entries that were extracted from the column (may be less than are really in the column if n was too
 * small)
 */
uint32_t
jmtx_matrix_brm_get_col(const jmtx_matrix_brm* mtx, uint32_t col, float values[]);

/**
 * Creates a transpose of a matrix
 * @param mtx pointer to the memory where the input matrix is stored
 * @param p_out address where the pointer to the output matrix will be returned
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtx_matrix_brm_transpose(
        const jmtx_matrix_brm* mtx, jmtx_matrix_brm** p_out, const jmtx_allocator_callbacks* allocator_callbacks);

/**
 * Creates a copy of the matrix
 * @param mtx pointer to the memory where the input matrix is stored
 * @param p_out address where the pointer to the output matrix will be returned
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtx_matrix_brm_copy(const jmtx_matrix_brm* mtx, jmtx_matrix_brm** p_out, const jmtx_allocator_callbacks* allocator_callbacks);


/**
 * Computes one entry of Ax. This function only computes the i-th entry to make it possible to compute it in parallel.
 * @param mtx pointer to the memory where the matrix A is stored.
 * @param x pointer to the memory where the vector x is stored
 * @param i what entry of the residual to compute
 * @return result of the multiplication
 */
float jmtx_matrix_brm_vector_multiply_row(const jmtx_matrix_brm* mtx, const float* x, uint32_t i);


#endif //JMTX_BAND_ROW_MAJOR_H
