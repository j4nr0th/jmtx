//
// Created by jan on 3.1.2024.
//

/**
 * Functions declared here perform minimum checking of the parameters and assume that values that were passed to them
 * were valid (matrix have proper dimensions, indices are within required bounds, etc.). "Safe" versions of these
 * functions, which do perform parameter validation are in the <TBD> header.
 */
#ifndef JMTXZ_DENSE_ROW_MAJOR_H
#define JMTXZ_DENSE_ROW_MAJOR_H
#ifndef JMTX_MATRIX_BASE_H
    #include "../../matrix_base.h"
#endif
/**
 * @paragraph
 * Dense Row-Major matrix (DRM) is a matrix which has only a few or no zeros. As such, all entries are stored in memory
 * in a row-contiguous way. This allows for fast row access, but slower column access.
 *
 * @paragraph
 * Main advantage of DRM matrices is the constant time random access to any and all elements. With all elements possible
 * to utilize, it is also possible to compute decompositions exactly, making the useful for solving smaller system in an
 * exact way.
 */
typedef struct jmtxz_matrix_drm_struct jmtxz_matrix_drm;


/**
 * Initializes a new Dense Row-Major matrix
 * @param p_mtx address that receives the pointer to the matrix
 * @param rows number of rows of the matrix
 * @param cols number of columns of the matrix
 * @param set_value pointer to value with which to initialize all entries. If NULL, then matrix is left uninitialized
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxz_matrix_drm_new(
 jmtxz_matrix_drm** p_mtx, uint32_t rows, uint32_t cols, const _Complex double* set_value,
 const jmtx_allocator_callbacks* allocator_callbacks);


/**
 * Cleans up the DRM matrix and frees all of its memory
 * @param mtx pointer to memory where the matrix is stored
 */
void jmtxz_matrix_drm_destroy(jmtxz_matrix_drm* mtx);

/**
 * Sets the row of the matrix. More efficient than setting it element by element
 * @param mtx pointer to the memory where the matrix is stored
 * @param row index of the row to set
 * @param values values of entries
 */
void jmtxz_matrix_drm_set_row(const jmtxz_matrix_drm* mtx, uint32_t row, const _Complex double values[]);

/**
 * Sets the column of the matrix. More efficient than setting it element by element
 * @param mtx pointer to the memory where the matrix is stored
 * @param col index of the column to set
 * @param values values of entries
 */
void jmtxz_matrix_drm_set_col(const jmtxz_matrix_drm* mtx, uint32_t col, const _Complex double values[]);

/**
 * Returns the pointers to arrays of column indices and element values for that row
 * @param mtx pointer to the memory where the matrix is stored
 * @param row index of the row to get
 * @param p_elements pointer to array of values
 * @return number of elements in the row, which is the number of valid elements in arrays given to p_indices and
 * p_elements
 */
uint_fast32_t jmtxz_matrix_drm_get_row(const jmtxz_matrix_drm* mtx, uint32_t row, _Complex double* p_elements[1]);

/**
 * Multiplies a dense column vector x by the sparse matrix and stores the result at y
 * @param mtx pointer to the memory where the matrix is stored
 * @param x pointer to vector to be multiplied
 * @param y pointer to vector where the result of multiplication is to be stored
 */
void jmtxz_matrix_drm_vector_multiply(const jmtxz_matrix_drm* mtx, const _Complex double* restrict x, _Complex double* restrict y);

/**
 * Sets a single entry in the matrix.
 * @param mtx pointer to the memory where the matrix is stored
 * @param i row index
 * @param j column index
 * @param value value to which the value is set
 */
void jmtxz_matrix_drm_set_entry(const jmtxz_matrix_drm* mtx, uint32_t i, uint32_t j, _Complex double value);

/**
 * Returns a single entry from the matrix.
 * @param mtx pointer to the memory where the matrix is stored
 * @param i row index
 * @param j column index
 * @return value of the entry
 */
_Complex double jmtxz_matrix_drm_get_entry(const jmtxz_matrix_drm* mtx, uint32_t i, uint32_t j);

/**
 * Returns a pointer to an entry in the matrix.
 * @param mtx pointer to the memory where the matrix is stored
 * @param i row index
 * @param j column index
 * @return pointer to an entry
 */
_Complex double* jmtxz_matrix_drm_entry_ptr(const jmtxz_matrix_drm* mtx, uint32_t i, uint32_t j);


/**
 * Zeros all entries within a matrix, but does not remove them in case they need to be reused
 * @param mtx matrix to zero
 */
void jmtxz_matrix_drm_zero_all_entries(const jmtxz_matrix_drm* mtx);

/**
 * Similar to jmtxz_matrix_drm_set_all_entries, but slower, since it can not use memset. On the other hand, it allows for
 * the value to be other than 0
 * @param mtx matrix to set
 * @param x value to which to set all entries to
 */
void jmtxz_matrix_drm_set_all_entries(const jmtxz_matrix_drm* mtx, _Complex double x);

/**
 * Returns the values of entries in the matrix, along with what row of the matrix they were located in
 * @param mtx pointer to the memory where the matrix is stored
 * @param col column index of the matrix to look at
 * @param p_values a buffer of at least n values which receives the values of the column
 * @return number of entries that were extracted from the column (may be less than are really in the column if n was too
 * small)
 */
uint32_t
jmtxz_matrix_drm_get_col(const jmtxz_matrix_drm* mtx, uint32_t col, _Complex double values[]);

/**
 * Creates a transpose of a matrix
 * @param mtx pointer to the memory where the input matrix is stored
 * @param p_out address where the pointer to the output matrix will be returned
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxz_matrix_drm_transpose(
        const jmtxz_matrix_drm* mtx, jmtxz_matrix_drm** p_out, const jmtx_allocator_callbacks* allocator_callbacks);

/**
 * Creates a transpose of a matrix
 * @param mtx pointer to the memory where the input matrix is stored
 * @param p_out address where the pointer to the output matrix will be returned
 * @param aux_row auxiliary memory to use for intermediate storage
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result jmtxz_matrix_drm_transpose_inplace(jmtxz_matrix_drm* mtx, _Complex double* aux_row);

/**
 * Creates a copy of the matrix
 * @param mtx pointer to the memory where the input matrix is stored
 * @param p_out address where the pointer to the output matrix will be returned
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result jmtxz_matrix_drm_copy(const jmtxz_matrix_drm* mtx, jmtxz_matrix_drm** p_out, const jmtx_allocator_callbacks* allocator_callbacks);


/**
 * Computes one entry of Ax. This function only computes the i-th entry to make it possible to compute it in parallel.
 * @param mtx pointer to the memory where the matrix A is stored.
 * @param x pointer to the memory where the vector x is stored
 * @param i what entry of the residual to compute
 * @return result of the multiplication
 */
_Complex double jmtxz_matrix_drm_vector_multiply_row(const jmtxz_matrix_drm* mtx, const _Complex double* x, uint32_t i);

/**
 * Returns the upper bandwidth and the lower bandwidth of the BRM matrix
 * @param mtx pointer to the memory where the matrix is stored
 * @param ubw pointer which receives the upper bandwidth of the matrix
 * @param lbw pointer which receives the lower bandwidth of the matrix
 */
void jmtxz_matrix_drm_get_bandwidths(const jmtxz_matrix_drm* mtx, uint32_t* ubw, uint32_t* lbw);

/**
 * Exchanges row of matrix through the use of a permutation matrix. This does not cause memory copying and is thus fast
 * @param mtx pointer to the memory where the matrix is stored
 * @param row1 index of a row to exchange with row2
 * @param row2 index of a row to exchange with row1
 */
jmtx_result jmtxz_matrix_drm_swap_rows(jmtxz_matrix_drm* mtx, uint32_t row1, uint32_t row2);

/**
 * Permutes all rows at once. The provided list must contain all entries in the set [0, n) exactly once, where n is the
 * number or rows of the matrix.
 * @param mtx pointer to the memory where the matrix is stored
 * @param perm list of permutation indices
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on allocation failure
 */
jmtx_result jmtxz_matrix_drm_set_permutation(jmtxz_matrix_drm* mtx, const uint32_t* perm);

/**
 * Executes permutations of matrix rows and reorders rows in memory. This may allow for better memory access on a
 * permuted matrix. Should not be done on decompositions. This function requires additional memory, but that allows
 * it to swap whole rows at once, which makes it more cache friendly and potentially faster for large matrices.
 * @param mtx pointer to the memory where the matrix is stored
 * @param aux_row memory that can be used by the function to store an intermediate matrix row
 */
void jmtxz_matrix_drm_commit_permutations(jmtxz_matrix_drm* mtx);

/**
 * Executes permutations of matrix rows and reorders rows in memory. This may allow for better memory access on a
 * permuted matrix. Should not be done on decompositions. This function requires additional memory, but that allows
 * it to swap whole rows at once, which makes it more cache friendly and potentially faster for large matrices.
 * @param mtx pointer to the memory where the matrix is stored
 * @param aux_row memory that can be used by the function to store an intermediate matrix row
 */
void jmtxz_matrix_drm_commit_permutations2(jmtxz_matrix_drm* mtx, _Complex double* aux_row);

#endif //JMTXZ_DENSE_ROW_MAJOR_H
