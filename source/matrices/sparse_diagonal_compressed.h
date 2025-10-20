// Automatically generated from include/jmtx/float/matrices/sparse_diagonal_compressed.h on Fri Dec  1 06:43:05 2023
//
// Created by jan on 27.11.2023.
//

#ifndef JMTX_SPARSE_DIAGONAL_COMPRESSED_H
#define JMTX_SPARSE_DIAGONAL_COMPRESSED_H
#include "../matrix_base.h"
#include "../matrix_base_internal.h"

typedef struct JMTX_NAME_TYPED(matrix_cds_struct) JMTX_NAME_TYPED(matrix_cds);

struct JMTX_NAME_TYPED(matrix_cds_diagonal_array_t)
{
    JMTX_INDEX_T capacity;     //  Number of diagonals
    JMTX_INDEX_T count;        //  Capacity for diagonals
    JMTX_INDEX_T *indices;     //  Indices of the diagonals relative to the main
    JMTX_SCALAR_T **diagonals; //  Array of the diagonal array pointers. Valid for count, space for up to capacity.
};
typedef struct JMTX_NAME_TYPED(matrix_cds_diagonal_array_t) JMTX_NAME_TYPED(matrix_cds_diagonal_array);

struct JMTX_NAME_TYPED(matrix_cds_struct)
{
    jmtx_matrix_base base;
    JMTX_NAME_TYPED(matrix_cds_diagonal_array) super_diagonals;
    JMTX_NAME_TYPED(matrix_cds_diagonal_array) sub_diagonals;
    JMTX_SCALAR_T *restrict main_diagonal;
};

/**
 * Initializes a new Compressed Diagonal Sparse matrix
 * @param p_mtx address that receives the pointer to the matrix
 * @param rows number of rows of the sparse matrix
 * @param cols number of columns of the sparse matrix
 * @param n_diagonals how already filled diagonals are given to the matrix to initialize with
 * @param p_dia_idx SORTED indices of locations of the diagonals to reserve, indices being relative to the main
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result JMTX_NAME_TYPED(matrix_cds_new)(JMTX_NAME_TYPED(matrix_cds) * *p_mtx, JMTX_INDEX_T rows, JMTX_INDEX_T cols,
                                            JMTX_INDEX_T n_diagonals,
                                            const int32_t p_dia_idx[JMTX_ARRAY_ATTRIB(static n_diagonals)],
                                            const jmtx_allocator_callbacks *allocator_callbacks);

/**
 * Cleans up the cds matrix and frees all of its memory
 * @param mtx pointer to memory where the matrix is stored
 */
void JMTX_NAME_TYPED(matrix_cds_destroy)(JMTX_NAME_TYPED(matrix_cds) * mtx);

/**
 * Returns the total number of diagonals in the matrix
 * @param mtx matrix for which to return this value for
 * @return the number of non-zero diagonals
 */
JMTX_FAST_INT_T JMTX_NAME_TYPED(matrix_cds_diagonal_count)(const JMTX_NAME_TYPED(matrix_cds) * mtx);

/**
 * Returns the number of entries in the diagonal of the matrix
 * @param mtx pointer to the memory where the matrix is stored
 * @param dia index of the diagonal of the matrix to look at
 * @return number of entries in the diagonal
 */
JMTX_FAST_INT_T JMTX_NAME_TYPED(matrix_cds_entries_in_dia)(JMTX_NAME_TYPED(matrix_cds) * mtx, int32_t dia);

/**
 * Sets the entire diagonal to the values provided in the array. Will allocate a new diagonal if it was previously not
 * @param mtx pointer to the memory where the matrix is stored
 * @param dia_idx index of the diagonal to set
 * @param values values of non-zero values
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result JMTX_NAME_TYPED(matrix_cds_set_diagonal_full)(JMTX_NAME_TYPED(matrix_cds) * mtx, int32_t dia_idx,
                                                          const JMTX_SCALAR_T values[]);

/**
 * Sets the specified number of entries in the diagonal, starting at the given offset for up to n entries at most
 * @param mtx pointer to the memory where the matrix is stored
 * @param dia_idx index of the diagonal to set
 * @param offset offset from the start of the diagonal to start at
 * @param n maximum number of entries to set
 * @param p_count (optional) pointer which receives the number of entries that were actually set (may be less than n if
 * it would go out of bounds)
 * @param values values of non-zero values
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result JMTX_NAME_TYPED(matrix_cds_set_diagonal_part)(JMTX_NAME_TYPED(matrix_cds) * mtx, int32_t dia_idx,
                                                          JMTX_INDEX_T offset, JMTX_INDEX_T n, JMTX_INDEX_T *p_count,
                                                          const JMTX_SCALAR_T values[JMTX_ARRAY_ATTRIB(static n)]);

/**
 * Allocates a new diagonal if one does not already exits, otherwise it returns the pointer to the existing one.
 * Creating a new one does not initialize it.
 * @param mtx matrix which to allocate the diagonal for
 * @param dia offset of the diagonal from the main diagonal
 * @param p_size pointer which receives the number of elements in the allocated diagonal. May be NULL
 * @return pointer to the newly allocated diagonal, or NULL in case it failed to do so
 */
JMTX_SCALAR_T *JMTX_NAME_TYPED(matrix_cds_allocate_diagonal)(JMTX_NAME_TYPED(matrix_cds) * mtx, int32_t dia,
                                                             JMTX_INDEX_T *p_size);

/**
 * Allocates a new diagonal if one does not already exits, otherwise it returns the pointer to the existing one.
 * Sets the diagonal to zero if it was not allocated before.
 * @param mtx matrix which to allocate the diagonal for
 * @param dia offset of the diagonal from the main diagonal
 * @param p_size pointer which receives the number of elements in the allocated diagonal. May be NULL
 * @return pointer to the newly allocated diagonal, or NULL in case it failed to do so
 */
JMTX_SCALAR_T *JMTX_NAME_TYPED(matrix_cds_allocate_zero_diagonal)(JMTX_NAME_TYPED(matrix_cds) * mtx, int32_t dia,
                                                                  JMTX_INDEX_T *p_size);

/**
 * Sets a row of the matrix to zero. Does not free any memory and does it as efficiently as possible
 * @param mtx matrix the row of which to clear
 * @param row the index of the row to clear
 */
void JMTX_NAME_TYPED(matrix_cds_zero_row)(const JMTX_NAME_TYPED(matrix_cds) * mtx, JMTX_INDEX_T row);

/**
 * Sets a column of the matrix to zero. Does not free any memory and does it as efficiently as possible
 * @param mtx matrix the column of which to clear
 * @param col the index of the column to clear
 */
void JMTX_NAME_TYPED(matrix_cds_zero_col)(const JMTX_NAME_TYPED(matrix_cds) * mtx, JMTX_INDEX_T col);

/**
 * Returns the pointer to the diagonal if one already exits, otherwise it returns NULL
 * @param mtx matrix which to allocate the diagonal for
 * @param dia offset of the diagonal from the main diagonal
 * @return pointer to the diagonal, or NULL in case it does not exist
 */
JMTX_SCALAR_T *JMTX_NAME_TYPED(matrix_cds_get_diagonal)(const JMTX_NAME_TYPED(matrix_cds) * mtx, int32_t dia);

/**
 * Returns the number of entries in the row of the matrix
 * @param mtx pointer to the memory where the matrix is stored
 * @param row row index of the matrix to look at
 * @return number of entries in the row
 */
JMTX_INDEX_T JMTX_NAME_TYPED(matrix_cds_entries_in_row)(const JMTX_NAME_TYPED(matrix_cds) * mtx, JMTX_INDEX_T row);

/**
 * Returns the values of entries in the matrix, along with what column of the matrix they were located in
 * @param mtx pointer to the memory where the matrix is stored
 * @param row row index of the matrix to look at
 * @param n number of values in the row to be extracted at most
 * @param p_values a buffer of at least n values which receives the values of the row
 * @param p_cols a buffer of at least n values which receives the column indices of the row
 * @return number of entries that were extracted from the column (may be less than are really in the column if n was too
 * small)
 */
JMTX_INDEX_T JMTX_NAME_TYPED(matrix_cds_get_row)(const JMTX_NAME_TYPED(matrix_cds) * mtx, JMTX_INDEX_T row,
                                                 JMTX_INDEX_T n, JMTX_SCALAR_T p_values[JMTX_ARRAY_ATTRIB(restrict n)],
                                                 JMTX_INDEX_T p_cols[JMTX_ARRAY_ATTRIB(restrict n)]);

/**
 * Returns the values of entries in the matrix, along with what column of the matrix they were located in
 * @param mtx pointer to the memory where the matrix is stored
 * @param row row index of the matrix to look at
 * @param n number of values to be set in the row
 * @param p_values a buffer of at least n values which contain the values to be set to
 * @param p_cols a buffer of at least n values which receives the column indices of the row
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result JMTX_NAME_TYPED(matrix_cds_set_row)(JMTX_NAME_TYPED(matrix_cds) * mtx, JMTX_INDEX_T row, JMTX_INDEX_T n,
                                                const JMTX_SCALAR_T p_values[JMTX_ARRAY_ATTRIB(restrict static n)],
                                                const JMTX_INDEX_T p_cols[JMTX_ARRAY_ATTRIB(restrict static n)]);

/**
 * Returns the number of entries in the column of the matrix
 * @param mtx pointer to the memory where the matrix is stored
 * @param col column index of the matrix to look at
 * @return number of entries in the column
 */
JMTX_INDEX_T JMTX_NAME_TYPED(matrix_cds_entries_in_col)(const JMTX_NAME_TYPED(matrix_cds) * mtx, JMTX_INDEX_T col);

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
JMTX_INDEX_T JMTX_NAME_TYPED(matrix_cds_get_col)(const JMTX_NAME_TYPED(matrix_cds) * mtx, JMTX_INDEX_T col,
                                                 JMTX_INDEX_T n, JMTX_SCALAR_T p_values[JMTX_ARRAY_ATTRIB(restrict n)],
                                                 JMTX_INDEX_T p_rows[JMTX_ARRAY_ATTRIB(restrict n)]);

/**
 * Sets all the entries in the column of the matrix, zeroing non-specified entries and allocating new diagonals as
 * needed
 * @param mtx pointer to the memory where the matrix is stored
 * @param col column index of the matrix to look at
 * @param n number of values in the column to be set
 * @param p_values a buffer of n values of column entries
 * @param p_rows a buffer of n indices of column indices of columns to set
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on allocation failure
 */
jmtx_result JMTX_NAME_TYPED(matrix_cds_set_col)(JMTX_NAME_TYPED(matrix_cds) * mtx, JMTX_INDEX_T col, JMTX_INDEX_T n,
                                                const JMTX_SCALAR_T p_values[JMTX_ARRAY_ATTRIB(restrict static n)],
                                                const JMTX_INDEX_T p_rows[JMTX_ARRAY_ATTRIB(restrict static n)]);

/**
 * Multiplies a dense column vector x by the sparse matrix and stores the result at y
 * @param mtx pointer to the memory where the matrix is stored
 * @param x pointer to vector to be multiplied
 * @param y pointer to vector where the result of multiplication is to be stored
 */
void JMTX_NAME_TYPED(matrix_cds_vector_multiply)(const JMTX_NAME_TYPED(matrix_cds) * mtx,
                                                 const JMTX_SCALAR_T *restrict x, JMTX_SCALAR_T *restrict y);

/**
 * Sets a single entry in the matrix. The diagonal that it is to be inserted on MUST exist before
 * @param mtx pointer to the memory where the matrix is stored
 * @param i row index
 * @param j column index
 * @param value value to which the value is set
 */
void JMTX_NAME_TYPED(matrix_cds_set_entry)(const JMTX_NAME_TYPED(matrix_cds) * mtx, JMTX_INDEX_T i, JMTX_INDEX_T j,
                                           JMTX_SCALAR_T value);

/**
 * Returns a single entry from the matrix.
 * @param mtx pointer to the memory where the matrix is stored
 * @param i row index
 * @param j column index
 * @return value of the entry (0 if the entry was not manually set to anything else)
 */
JMTX_SCALAR_T JMTX_NAME_TYPED(matrix_cds_get_entry)(const JMTX_NAME_TYPED(matrix_cds) * mtx, JMTX_INDEX_T i,
                                                    JMTX_INDEX_T j);

/**
 * Inserts and entry into the matrix, even if diagonal was not present before in the matrix
 * @param mtx pointer to the memory where the matrix is stored
 * @param i row index
 * @param j column index
 * @param value value to which the value is to be added
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result JMTX_NAME_TYPED(matrix_cds_insert_entry)(JMTX_NAME_TYPED(matrix_cds) * mtx, JMTX_INDEX_T i, JMTX_INDEX_T j,
                                                     JMTX_SCALAR_T value);

/**
 * Adds a value to an entry in the matrix when it exists or sets it to that value if it does not.
 * @param mtx pointer to the memory where the matrix is stored
 * @param i row index
 * @param j column index
 * @param value value to which the value is to be added
 */
jmtx_result JMTX_NAME_TYPED(matrix_cds_add_to_entry)(JMTX_NAME_TYPED(matrix_cds) * mtx, JMTX_INDEX_T i, JMTX_INDEX_T j,
                                                     JMTX_SCALAR_T value);

/**
 * Counts the number of times a specific value occurs in the matrix
 * @param mtx pointer to the memory where the matrix is stored
 * @param v value which to search for
 * @return number of times the value appeared in the matrix
 */
JMTX_INDEX_T JMTX_NAME_TYPED(matrix_cds_count_values)(const JMTX_NAME_TYPED(matrix_cds) * mtx, JMTX_SCALAR_T v);

/**
 * Zeros all entries within a matrix, but does not remove them in case they need to be reused
 * @param mtx matrix to zero
 */
void JMTX_NAME_TYPED(matrix_cds_zero_all_entries)(const JMTX_NAME_TYPED(matrix_cds) * mtx);

/**
 * Similar to JMTX_NAME_TYPED(matrix_cds_zero_all_entries, but slower, since it can not use memset. On the other hand,
 * it allows for the value to be other than 0
 * @param mtx matrix to set
 * @param x value to which to set all entries to
 */
void JMTX_NAME_TYPED(matrix_cds_set_all_entries)(const JMTX_NAME_TYPED(matrix_cds) * mtx, JMTX_SCALAR_T x);

/**
 * Keeps the memory for the matrix, but sets the entry count to 0, so that matrix can be rebuilt.
 * @param mtx matrix to clear
 */
void JMTX_NAME_TYPED(matrix_cds_clear)(JMTX_NAME_TYPED(matrix_cds) * mtx);

/**
 * Creates a transpose of a matrix
 * @param mtx pointer to the memory where the input matrix is stored
 * @param p_out address where the pointer to the output matrix will be returned
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result JMTX_NAME_TYPED(matrix_cds_transpose)(const JMTX_NAME_TYPED(matrix_cds) * mtx,
                                                  JMTX_NAME_TYPED(matrix_cds) * *p_out,
                                                  const jmtx_allocator_callbacks *allocator_callbacks);

/**
 * Creates a copy of the matrix
 * @param mtx pointer to the memory where the input matrix is stored
 * @param p_out address where the pointer to the output matrix will be returned
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result JMTX_NAME_TYPED(matrix_cds_copy)(const JMTX_NAME_TYPED(matrix_cds) * mtx,
                                             JMTX_NAME_TYPED(matrix_cds) * *p_out,
                                             const jmtx_allocator_callbacks *allocator_callbacks);

/**
 * Sets all the values on the diagonal to same value
 * @param mtx matrix to remove the diagonal from from
 * @param dia_idx diagonal to remove
 * @param v value to set the diagonal to
 */
void JMTX_NAME_TYPED(matrix_cds_set_diagonal)(JMTX_NAME_TYPED(matrix_cds) * mtx, int32_t dia_idx, JMTX_SCALAR_T v);

/**
 * Sets all the values on the diagonal to zero
 * @param mtx matrix to remove the diagonal from from
 * @param dia_idx diagonal to remove
 */
void JMTX_NAME_TYPED(matrix_cds_zero_diagonal)(JMTX_NAME_TYPED(matrix_cds) * mtx, int32_t dia_idx);

#endif // JMTX_SPARSE_DIAGONAL_COMPRESSED_H
