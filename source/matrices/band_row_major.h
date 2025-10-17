#ifndef JMTX_BAND_ROW_MAJOR_H
#define JMTX_BAND_ROW_MAJOR_H
#include "../matrix_base.h"
#include "../matrix_base_internal.h"

struct JMTX_NAME_TYPED(matrix_brm_struct)
{
    jmtx_matrix_base base;
    JMTX_INDEX_T upper_bandwidth;
    JMTX_INDEX_T lower_bandwidth;
    JMTX_SCALAR_T *restrict values;
};

/**
 * @paragraph
 * Band Row-Major matrix (BRM) is a matrix, which has constant upper bandwidth (ubw) and lower bandwidths (lbw):
 * constant number of entries above and below the diagonal. The rows are stored contiguously. The memory required is
 * bounded by lbw + 1 + ubw. Memory access is performed in constant time for both rows and columns, with rows being
 * slightly faster and doesn't need a separate function for setting it, since a pointer to a row allows modification.
 *
 * @paragraph
 * Main advantage of BRM matrices is the fact that for a BRM with bandwidths ubw and lbw, its exact LU decomposition
 * results in the matrix L having bandwidths 0 and lbw, and matrix U having bandwidths ubw and 0. This means that for
 * a chosen matrix number of elements that need to be computed for full LU decomposition becomes (lbw + 1 + ubw) * N,
 * with N being the size of the matrix. This means that using full LU decomposition can be viable, given that the
 * bandwidth is low (most PDEs on a 1D domain).
 */
typedef struct JMTX_NAME_TYPED(matrix_brm_struct) JMTX_NAME_TYPED(matrix_brm);

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
jmtx_result JMTX_NAME_TYPED(matrix_brm_new)(JMTX_NAME_TYPED(matrix_brm) * *p_mtx, JMTX_INDEX_T rows, JMTX_INDEX_T cols,
                                            JMTX_INDEX_T ubw, JMTX_INDEX_T lbw, const JMTX_SCALAR_T *set_value,
                                            const jmtx_allocator_callbacks *allocator_callbacks);

/**
 * Cleans up the BRM matrix and frees all of its memory
 * @param mtx pointer to memory where the matrix is stored
 */
void JMTX_NAME_TYPED(matrix_brm_destroy)(JMTX_NAME_TYPED(matrix_brm) * mtx);

/**
 * Sets the row of the matrix. More efficient than setting it element by element
 * @param mtx pointer to the memory where the matrix is stored
 * @param row index of the row to set
 * @param values values of entries
 */
void JMTX_NAME_TYPED(matrix_brm_set_row)(const JMTX_NAME_TYPED(matrix_brm) * mtx, JMTX_INDEX_T row,
                                         const JMTX_SCALAR_T values[]);

/**
 * Sets the column of the matrix. More efficient than setting it element by element
 * @param mtx pointer to the memory where the matrix is stored
 * @param col index of the column to set
 * @param values values of entries
 */
void JMTX_NAME_TYPED(matrix_brm_set_col)(const JMTX_NAME_TYPED(matrix_brm) * mtx, JMTX_INDEX_T col,
                                         const JMTX_SCALAR_T values[]);

/**
 * Returns the pointers to arrays of column indices and element values for that row
 * @param mtx pointer to the memory where the matrix is stored
 * @param row index of the row to get
 * @param p_elements pointer to array of values
 * @return number of elements in the row, which is the number of valid elements in arrays given to p_indices and
 * p_elements
 */
JMTX_FAST_INT_T JMTX_NAME_TYPED(matrix_brm_get_row)(const JMTX_NAME_TYPED(matrix_brm) * mtx, JMTX_INDEX_T row,
                                                    JMTX_SCALAR_T *p_elements[1]);

/**
 * Returns the index of the first non-zero column in a specific row
 * @param mtx matrix for which this is to be determined
 * @param row row of the matrix for which this is to be determined
 * @return index of the first non-zero column in the row
 */
JMTX_FAST_INT_T JMTX_NAME_TYPED(matrix_brm_first_pos_in_row)(const JMTX_NAME_TYPED(matrix_brm) * mtx, JMTX_INDEX_T row);

/**
 * Returns the index of the last non-zero column in a specific row
 * @param mtx matrix for which this is to be determined
 * @param row row of the matrix for which this is to be determined
 * @return index of the last non-zero column in the row
 */
JMTX_FAST_INT_T JMTX_NAME_TYPED(matrix_brm_last_pos_in_row)(const JMTX_NAME_TYPED(matrix_brm) * mtx, JMTX_INDEX_T row);

/**
 * Returns the number of non-zero elements in a row of a matrix
 * @param mtx matrix for which this is to be determined
 * @param row row of the matrix for which this is to be determined
 * @return index of the last non-zero column in the row
 */
JMTX_FAST_INT_T JMTX_NAME_TYPED(matrix_brm_length_of_row)(const JMTX_NAME_TYPED(matrix_brm) * mtx, JMTX_INDEX_T row);

/**
 * Returns the index of the first non-zero row in a specific column
 * @param mtx matrix for which this is to be determined
 * @param col column of the matrix for which this is to be determined
 * @return index of the first non-zero row in the column
 */
JMTX_FAST_INT_T JMTX_NAME_TYPED(matrix_brm_first_pos_in_col)(const JMTX_NAME_TYPED(matrix_brm) * mtx, JMTX_INDEX_T col);

/**
 * Returns the index of the last non-zero row in a specific column
 * @param mtx matrix for which this is to be determined
 * @param col column of the matrix for which this is to be determined
 * @return index of the last non-zero row in the column
 */
JMTX_FAST_INT_T JMTX_NAME_TYPED(matrix_brm_last_pos_in_col)(const JMTX_NAME_TYPED(matrix_brm) * mtx, JMTX_INDEX_T col);

/**
 * Multiplies a dense column vector x by the sparse matrix and stores the result at y
 * @param mtx pointer to the memory where the matrix is stored
 * @param x pointer to vector to be multiplied
 * @param y pointer to vector where the result of multiplication is to be stored
 */
void JMTX_NAME_TYPED(matrix_brm_vector_multiply)(const JMTX_NAME_TYPED(matrix_brm) * mtx,
                                                 const JMTX_SCALAR_T *restrict x, JMTX_SCALAR_T *restrict y);

/**
 * Sets a single entry in the matrix. This is about as fast as setting the entire row of the matrix at once, if the
 * element was previously zero
 * @param mtx pointer to the memory where the matrix is stored
 * @param i row index
 * @param j column index
 * @param value value to which the value is set
 */
void JMTX_NAME_TYPED(matrix_brm_set_entry)(const JMTX_NAME_TYPED(matrix_brm) * mtx, JMTX_INDEX_T i, JMTX_INDEX_T j,
                                           JMTX_SCALAR_T value);

/**
 * Returns a single entry from the matrix.
 * @param mtx pointer to the memory where the matrix is stored
 * @param i row index
 * @param j column index
 * @return value of the entry (0 if the entry was not manually set to anything else)
 */
JMTX_SCALAR_T JMTX_NAME_TYPED(matrix_brm_get_entry)(const JMTX_NAME_TYPED(matrix_brm) * mtx, JMTX_INDEX_T i,
                                                    JMTX_INDEX_T j);

/**
 * Adds a value to an entry in the matrix when it exists or sets it to that value if it does not. This is about as
 * fast as setting the entire row of the matrix at once, if the entry was non-existent.
 * @param mtx pointer to the memory where the matrix is stored
 * @param i row index
 * @param j column index
 * @param value value to which the value is to be added
 */
void JMTX_NAME_TYPED(matrix_brm_add_to_entry)(const JMTX_NAME_TYPED(matrix_brm) * mtx, JMTX_INDEX_T i, JMTX_INDEX_T j,
                                              JMTX_SCALAR_T value);

/**
 * Counts the number of times a specific value occurs in the matrix
 * @param mtx pointer to the memory where the matrix is stored
 * @param v value which to search for
 * @return number of times the value appeared in the matrix
 */
JMTX_INDEX_T JMTX_NAME_TYPED(matrix_brm_count_values)(const JMTX_NAME_TYPED(matrix_brm) * mtx, JMTX_SCALAR_T v);

/**
 * Zeros all entries within a matrix, but does not remove them in case they need to be reused
 * @param mtx matrix to zero
 */
void JMTX_NAME_TYPED(matrix_brm_zero_all_entries)(const JMTX_NAME_TYPED(matrix_brm) * mtx);

/**
 * Similar to JMTX_NAME_TYPED(matrix_brm_zero_all_entries, but slower, since it can not use memset. On the other hand,
 * it allows for the value to be other than 0
 * @param mtx matrix to set
 * @param x value to which to set all entries to
 */
void JMTX_NAME_TYPED(matrix_brm_set_all_entries)(const JMTX_NAME_TYPED(matrix_brm) * mtx, JMTX_SCALAR_T x);

/**
 * Returns the number of entries in the column of the matrix
 * @param mtx pointer to the memory where the matrix is stored
 * @param col column index of the matrix to look at
 * @return number of entries in the column
 */
JMTX_INDEX_T JMTX_NAME_TYPED(matrix_brm_length_of_col)(const JMTX_NAME_TYPED(matrix_brm) * mtx, JMTX_INDEX_T col);

/**
 * Returns the values of entries in the matrix, along with what row of the matrix they were located in
 * @param mtx pointer to the memory where the matrix is stored
 * @param col column index of the matrix to look at
 * @param values a buffer of at least n values which receives the values of the column
 * @return number of entries that were extracted from the column (may be less than are really in the column if n was too
 * small)
 */
JMTX_INDEX_T JMTX_NAME_TYPED(matrix_brm_get_col)(const JMTX_NAME_TYPED(matrix_brm) * mtx, JMTX_INDEX_T col,
                                                 JMTX_SCALAR_T values[]);

/**
 * Creates a transpose of a matrix
 * @param mtx pointer to the memory where the input matrix is stored
 * @param p_out address where the pointer to the output matrix will be returned
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result JMTX_NAME_TYPED(matrix_brm_transpose)(const JMTX_NAME_TYPED(matrix_brm) * mtx,
                                                  JMTX_NAME_TYPED(matrix_brm) * *p_out,
                                                  const jmtx_allocator_callbacks *allocator_callbacks);

/**
 * Creates a copy of the matrix
 * @param mtx pointer to the memory where the input matrix is stored
 * @param p_out address where the pointer to the output matrix will be returned
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result JMTX_NAME_TYPED(matrix_brm_copy)(const JMTX_NAME_TYPED(matrix_brm) * mtx,
                                             JMTX_NAME_TYPED(matrix_brm) * *p_out,
                                             const jmtx_allocator_callbacks *allocator_callbacks);

/**
 * Computes one entry of Ax. This function only computes the i-th entry to make it possible to compute it in parallel.
 * @param mtx pointer to the memory where the matrix A is stored.
 * @param x pointer to the memory where the vector x is stored
 * @param i what entry of the residual to compute
 * @return result of the multiplication
 */
JMTX_SCALAR_T JMTX_NAME_TYPED(matrix_brm_vector_multiply_row)(const JMTX_NAME_TYPED(matrix_brm) * mtx,
                                                              const JMTX_SCALAR_T *x, JMTX_INDEX_T i);

/**
 * Returns the upper bandwidth and the lower bandwidth of the BRM matrix
 * @param mtx pointer to the memory where the matrix is stored
 * @param ubw pointer which receives the upper bandwidth of the matrix
 * @param lbw pointer which receives the lower bandwidth of the matrix
 */
void JMTX_NAME_TYPED(matrix_brm_get_bandwidths)(const JMTX_NAME_TYPED(matrix_brm) * mtx, JMTX_INDEX_T *ubw,
                                                JMTX_INDEX_T *lbw);

#endif //  JMTX_BAND_ROW_MAJOR_H
