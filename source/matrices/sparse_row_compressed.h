#ifndef JMTX_SPARSE_ROW_COMPRESSED_H
#define JMTX_SPARSE_ROW_COMPRESSED_H
#include "../matrix_base.h"

struct JMTX_NAME_TYPED(matrix_crs_struct)
{
    jmtx_matrix_base base;
    //  end_of_row_offsets[i] has the number of values are there before the end of the row i
    JMTX_INDEX_T *end_of_row_offsets;
    //  Column indices corresponding with the individual values
    JMTX_INDEX_T *indices;
    //  Values of values
    JMTX_SCALAR_T *values;
    JMTX_SIZE_T n_entries;
    JMTX_SIZE_T capacity;
};
typedef struct JMTX_NAME_TYPED(matrix_crs_struct) JMTX_NAME_TYPED(matrix_crs);

/**
 * Initializes a new Compressed Row Sparse matrix
 * @param p_mtx address that receives the pointer to the matrix
 * @param rows number of rows of the sparse matrix
 * @param cols number of columns of the sparse matrix
 * @param reserved_entries how many entries should the space be reserved for in the matrix initially
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result JMTX_NAME_TYPED(matrix_crs_new)(JMTX_NAME_TYPED(matrix_crs) * *p_mtx, JMTX_INDEX_T rows, JMTX_INDEX_T cols,
                                            JMTX_INDEX_T reserved_entries,
                                            const jmtx_allocator_callbacks *allocator_callbacks);

/**
 * Cleans up the crs matrix and frees all of its memory
 * @param mtx pointer to memory where the matrix is stored
 */
void JMTX_NAME_TYPED(matrix_crs_destroy)(JMTX_NAME_TYPED(matrix_crs) * mtx);

/**
 * Frees up memory which the matrix is not currently using, which is was allocated in advance
 * @param mtx pointer to the memory where the matrix is stored
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result JMTX_NAME_TYPED(matrix_crs_shrink)(JMTX_NAME_TYPED(matrix_crs) * mtx);

/**
 * Sets the row of the matrix. More efficient than setting it element by element
 * @param mtx pointer to the memory where the matrix is stored
 * @param row index of the row to set
 * @param n how many entries are in the row
 * @param indices column indices of the entries, which has to be sorted from lowest to highest
 * @param values values of non-zero entries
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result JMTX_NAME_TYPED(matrix_crs_set_row)(JMTX_NAME_TYPED(matrix_crs) * mtx, JMTX_INDEX_T row, JMTX_INDEX_T n,
                                                const JMTX_INDEX_T indices[JMTX_ARRAY_ATTRIB(static n)],
                                                const JMTX_SCALAR_T values[JMTX_ARRAY_ATTRIB(static n)]);

/**
 * Version of JMTX_NAME_TYPED(matrix_crs_set_row which does not touch the count of entries after the current row. This
 * is useful when building a new matrix, as it avoids unnecessary setting and resetting of these entries. Must be called
 * for each row in order to ensure that the matrix is properly built. Makes no checks on the input parameters
 * @param mtx pointer to the memory where the matrix is stored
 * @param row index of the row to set
 * @param n how many entries are in the row
 * @param indices column indices of the values, which has to be sorted from lowest to highest
 * @param values values of non-zero values
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result JMTX_NAME_TYPED(matrix_crs_build_row)(JMTX_NAME_TYPED(matrix_crs) * mtx, JMTX_INDEX_T row, JMTX_INDEX_T n,
                                                  const JMTX_INDEX_T indices[JMTX_ARRAY_ATTRIB(static n)],
                                                  const JMTX_SCALAR_T values[JMTX_ARRAY_ATTRIB(static n)]);

/**
 * Returns the pointers to arrays of column indices and element values for that row
 * @param mtx pointer to the memory where the matrix is stored
 * @param row index of the row to get
 * @param p_indices pointer to column indices of values
 * @param p_elements pointer to values of values
 * @return number of elements in the row, which is the number of valid elements in arrays given to p_indices and
 * p_elements
 */
JMTX_INDEX_T JMTX_NAME_TYPED(matrix_crs_get_row)(const JMTX_NAME_TYPED(matrix_crs) * mtx, JMTX_INDEX_T row,
                                                 JMTX_INDEX_T *p_indices[1], JMTX_SCALAR_T *p_elements[1]);

/**
 * Multiplies a dense column vector x by the sparse matrix and stores the result at y
 * @param mtx pointer to the memory where the matrix is stored
 * @param x pointer to vector to be multiplied
 * @param y pointer to vector where the result of multiplication is to be stored
 */
void JMTX_NAME_TYPED(matrix_crs_vector_multiply)(const JMTX_NAME_TYPED(matrix_crs) * mtx,
                                                 const JMTX_SCALAR_T *restrict x, JMTX_SCALAR_T *restrict y);

/**
 * Sets a single entry in the matrix. This is about as fast as setting the entire row of the matrix at once, if the
 * element was previously zero
 * @param mtx pointer to the memory where the matrix is stored
 * @param i row index
 * @param j column index
 * @param value value to which the value is set
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result JMTX_NAME_TYPED(matrix_crs_set_entry)(JMTX_NAME_TYPED(matrix_crs) * mtx, JMTX_INDEX_T i, JMTX_INDEX_T j,
                                                  JMTX_SCALAR_T value);

/**
 * Returns a single entry from the matrix.
 * @param mtx pointer to the memory where the matrix is stored
 * @param i row index
 * @param j column index
 * @return value of the entry (0 if the entry was not manually set to anything else)
 */
JMTX_SCALAR_T JMTX_NAME_TYPED(matrix_crs_get_entry)(const JMTX_NAME_TYPED(matrix_crs) * mtx, JMTX_INDEX_T i,
                                                    JMTX_INDEX_T j);

/**
 * Adds a value to an entry in the matrix when it exists or sets it to that value if it does not. This is about as
 * fast as setting the entire row of the matrix at once, if the entry was non-existent.
 * @param mtx pointer to the memory where the matrix is stored
 * @param i row index
 * @param j column index
 * @param value value to which the value is to be added
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure when new entry was
 * inserted
 */
jmtx_result JMTX_NAME_TYPED(matrix_crs_add_to_entry)(JMTX_NAME_TYPED(matrix_crs) * mtx, JMTX_INDEX_T i, JMTX_INDEX_T j,
                                                     JMTX_SCALAR_T value);

/**
 * Counts the number of times a specific value occurs in the matrix
 * @param mtx pointer to the memory where the matrix is stored
 * @param v value which to search for
 * @return number of times the value appeared in the matrix
 */
JMTX_INDEX_T JMTX_NAME_TYPED(matrix_crs_count_values)(const JMTX_NAME_TYPED(matrix_crs) * mtx, JMTX_SCALAR_T v);

/**
 * Counts the number of times a specific column index occurs in the matrix
 * @param mtx pointer to the memory where the matrix is stored
 * @param v column index which to search for
 * @return number of times the column index appeared in the matrix
 */
JMTX_INDEX_T JMTX_NAME_TYPED(matrix_crs_count_indices)(const JMTX_NAME_TYPED(matrix_crs) * mtx, JMTX_INDEX_T v);

/**
 * Applies a unary function on the sparse matrix, only on its stored entries, which can be modified. If the user given
 * user function returns a non-zero value, that value is returned from the function and the iteration is stopped
 * @param mtx pointer to the memory where the matrix is stored
 * @param unary_fn user given function that is to be applied on the values
 * @param param optional parameter, which is passed to the unary_fn when called
 * @return JMTX_RESULT_SUCCESS if user function never returned non-zero, JMTX_RESULT_UNARY_RETURN as soon as the user
 * function returns non-zero
 */
jmtx_result JMTX_NAME_TYPED(matrix_crs_apply_unary_fn)(const JMTX_NAME_TYPED(matrix_crs) * mtx,
                                                       int (*unary_fn)(JMTX_INDEX_T i, JMTX_INDEX_T j,
                                                                       JMTX_SCALAR_T *p_value, void *param),
                                                       void *param);

/**
 * Removes entries exactly equal to zero. If the element is indeed zero, it is compared to (scalar_t)0.0
 * @param mtx pointer to the memory where the matrix is stored
 */
void JMTX_NAME_TYPED(matrix_crs_remove_zeros)(JMTX_NAME_TYPED(matrix_crs) * mtx);

/**
 * Removes entries which have absolute value less than specified value
 * @param mtx pointer to the memory where the matrix is stored
 * @param v value to which to compare it to
 */
void JMTX_NAME_TYPED(matrix_crs_remove_bellow_magnitude)(JMTX_NAME_TYPED(matrix_crs) * mtx, JMTX_REAL_T v);

/**
 * Zeros all entries within a matrix, but does not remove them in case they need to be reused
 * @param mtx matrix to zero
 */
void JMTX_NAME_TYPED(matrix_crs_zero_all_entries)(const JMTX_NAME_TYPED(matrix_crs) * mtx);

/**
 * Similar to JMTX_NAME_TYPED(matrix_crs_zero_all_entries, but slower, since it can not use memset. On the other hand,
 * it allows for the value to be other than 0
 * @param mtx matrix to set
 * @param x value to which to set all entries to
 */
void JMTX_NAME_TYPED(matrix_crs_set_all_entries)(const JMTX_NAME_TYPED(matrix_crs) * mtx, JMTX_SCALAR_T x);

/**
 * Keeps the memory for the matrix, but sets the entry count to 0, so that matrix can be rebuilt.
 * @param mtx matrix to clear
 */
void JMTX_NAME_TYPED(matrix_crs_clear)(JMTX_NAME_TYPED(matrix_crs) * mtx);

/**
 * Returns the number of entries in the column of the matrix
 * @param mtx pointer to the memory where the matrix is stored
 * @param col column index of the matrix to look at
 * @return number of entries in the column
 */
JMTX_INDEX_T JMTX_NAME_TYPED(matrix_crs_entries_in_col)(const JMTX_NAME_TYPED(matrix_crs) * mtx, JMTX_INDEX_T col);

/**
 * Returns the values of entries in the matrix, along with what row of the matrix they were located in
 * @param mtx pointer to the memory where the matrix is stored
 * @param col column index of the matrix to look at
 * @param n number of values in the column to be extracted at most
 * @param p_values a buffer of at least n values which receives the values of the column
 * @param p_rows a buffer of at least n values which receives the row indices of the column
 * @return number of entries that were extracted from the column (may be less than are really in the column if n was too
 * small)
 */
JMTX_INDEX_T JMTX_NAME_TYPED(matrix_crs_get_col)(const JMTX_NAME_TYPED(matrix_crs) * mtx, JMTX_INDEX_T col,
                                                 JMTX_INDEX_T n, JMTX_SCALAR_T p_values[JMTX_ARRAY_ATTRIB(n)],
                                                 JMTX_INDEX_T p_rows[JMTX_ARRAY_ATTRIB(n)]);

/**
 * Creates a transpose of a matrix
 * @param mtx pointer to the memory where the input matrix is stored
 * @param p_out address where the pointer to the output matrix will be returned
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result JMTX_NAME_TYPED(matrix_crs_transpose)(const JMTX_NAME_TYPED(matrix_crs) * mtx,
                                                  JMTX_NAME_TYPED(matrix_crs) * *p_out,
                                                  const jmtx_allocator_callbacks *allocator_callbacks);

/**
 * Creates a copy of the matrix
 * @param mtx pointer to the memory where the input matrix is stored
 * @param p_out address where the pointer to the output matrix will be returned
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result JMTX_NAME_TYPED(matrix_crs_copy)(const JMTX_NAME_TYPED(matrix_crs) * mtx,
                                             JMTX_NAME_TYPED(matrix_crs) * *p_out,
                                             const jmtx_allocator_callbacks *allocator_callbacks);

/**
 * Computes one entry of Ax. This function only computes the i-th entry to make it possible to compute it in parallel.
 * @param mtx pointer to the memory where the matrix A is stored.
 * @param x pointer to the memory where the vector x is stored
 * @param i what entry of the residual to compute
 * @return result of the multiplication
 */
JMTX_SCALAR_T JMTX_NAME_TYPED(matrix_crs_vector_multiply_row)(const JMTX_NAME_TYPED(matrix_crs) * mtx,
                                                              const JMTX_SCALAR_T *x, JMTX_INDEX_T i);

/**
 * Removes a row from a CRS matrix, reducing it's row count by 1
 * @param mtx matrix to remove the row from
 * @param row row to remove
 */
void JMTX_NAME_TYPED(matrix_crs_remove_row)(JMTX_NAME_TYPED(matrix_crs) * mtx, JMTX_INDEX_T row);

/**
 * Removes a column from a CRS matrix, reducing it's column count by 1. Slower than removing a row, so remove rows
 * before columns if both have to be done.
 * @param mtx matrix to remove the column from
 * @param col column to remove
 */
void JMTX_NAME_TYPED(matrix_crs_remove_column)(JMTX_NAME_TYPED(matrix_crs) * mtx, JMTX_INDEX_T col);

/**
 * Combines k matrices of size N_i x M into a single Sum(N_i) x M matrix by vertically stacking them
 * @param output receives pointer to the resulting matrix
 * @param allocator_callbacks memory allocators to use for the output matrix
 * @param k number of the matrices
 * @param matrix_list array of matrices to be joined together. Must have the same number of columns.
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_DIMS_MISMATCH if the number of columns is not the same for all
 * input matrices, return value of JMTX_NAME_TYPED(matrix_crs_new if that fails
 */
jmtx_result JMTX_NAME_TYPED(matrix_crs_join_vertically)(JMTX_NAME_TYPED(matrix_crs) * *output,
                                                        const jmtx_allocator_callbacks *allocator_callbacks, unsigned k,
                                                        const JMTX_NAME_TYPED(matrix_crs) *
                                                            matrix_list[JMTX_ARRAY_ATTRIB(static k)]);

jmtx_result JMTX_NAME_TYPED(matrix_crs_new_like)(const JMTX_NAME_TYPED(matrix_crs) * mtx,
                                                 JMTX_NAME_TYPED(matrix_crs) * *p_out,
                                                 const jmtx_allocator_callbacks *allocator_callbacks,
                                                 const JMTX_SCALAR_T *p_val);

/**
 * Finds the upper bandwidth of the matrix; what is the furthest distance of and entry above the main diagonal
 * @param mtx matrx to find the upper bandwidth of
 * @return upper bandwidth of the matrix
 */
JMTX_INDEX_T JMTX_NAME_TYPED(matrix_crs_find_upper_bandwidth)(const JMTX_NAME_TYPED(matrix_crs) * mtx);

/**
 * Finds the lower bandwidth of the matrix; what is the furthest distance of and entry below the main diagonal
 * @param mtx matrx to find the lower bandwidth of
 * @return lower bandwidth of the matrix
 */
JMTX_INDEX_T JMTX_NAME_TYPED(matrix_crs_find_lower_bandwidth)(const JMTX_NAME_TYPED(matrix_crs) * mtx);

#endif // JMTX_SPARSE_ROW_COMPRESSED_H
