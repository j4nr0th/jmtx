/**
 * Functions declared here perform minimum checking of the parameters and assume that values that were passed to them
 * were valid (matrix have proper dimensions, indices are within required bounds, etc.). "Safe" versions of these
 * functions, which do perform parameter validation are in the <TBD> header.
 */
#ifndef JMTX_DENSE_ROW_MAJOR_H
#define JMTX_DENSE_ROW_MAJOR_H
#include "../matrix_base.h"
#include "../matrix_base_internal.h"

struct JMTX_NAME_TYPED(matrix_drm_struct)
{
    jmtx_matrix_base base;
    //  Length this->base.rows, contains the order in which the rows are permuted (what row is read from where)
    JMTX_INDEX_T *permutations;
    //  Length this->base.rows, contains the reverse order in which the rows are permuted (what row is written to where)
    JMTX_INDEX_T *rperm;
    //  Length this->base.rows * this->base.cols, contains every entry in row-major ordering.
    JMTX_SCALAR_T *restrict values;
};
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
typedef struct JMTX_NAME_TYPED(matrix_drm_struct) JMTX_NAME_TYPED(matrix_drm);

/**
 * @brief Constructs a dense row-major matrix from the provided data.
 *
 * This function initializes a dense row-major matrix with specified dimensions and values.
 * It sets up the base structure with the appropriate type and dimensions. This matrix
 * cannot be destroyed, permuted, or resized, since allocator callbacks are not valid.
 *
 * The purpose of this function is to create inputs for functions by wrapping a raw pointer.
 *
 * @param rows The number of rows in the matrix.
 * @param cols The number of columns in the matrix.
 * @param values A pointer to an array of size (rows * cols) that contains the matrix elements in row-major order.
 *               The array must remain valid for the lifetime of the returned matrix.
 *
 * @return An initialized dense row-major matrix.
 */
JMTX_NAME_TYPED(matrix_drm)
JMTX_NAME_TYPED(matrix_drm_from_data)(JMTX_INDEX_T rows, JMTX_INDEX_T cols,
                                      JMTX_SCALAR_T values[JMTX_ARRAY_ATTRIB(static rows * cols)]);

/**
 * Initializes a new Dense Row-Major matrix
 * @param p_mtx address that receives the pointer to the matrix
 * @param rows number of rows of the matrix
 * @param cols number of columns of the matrix
 * @param set_value pointer to value with which to initialize all entries. If NULL, then matrix is left
 * uninitialized
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL
 * to use malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result JMTX_NAME_TYPED(matrix_drm_new)(JMTX_NAME_TYPED(matrix_drm) * *p_mtx, JMTX_INDEX_T rows, JMTX_INDEX_T cols,
                                            const JMTX_SCALAR_T *set_value,
                                            const jmtx_allocator_callbacks *allocator_callbacks);

/**
 * Cleans up the DRM matrix and frees all of its memory
 * @param mtx pointer to memory where the matrix is stored
 */
void JMTX_NAME_TYPED(matrix_drm_destroy)(JMTX_NAME_TYPED(matrix_drm) * mtx);

/**
 * Sets the row of the matrix. More efficient than setting it element by element
 * @param mtx pointer to the memory where the matrix is stored
 * @param row index of the row to set
 * @param values values of entries
 */
void JMTX_NAME_TYPED(matrix_drm_set_row)(const JMTX_NAME_TYPED(matrix_drm) * mtx, JMTX_INDEX_T row,
                                         const JMTX_SCALAR_T values[]);

/**
 * Sets the column of the matrix. More efficient than setting it element by element
 * @param mtx pointer to the memory where the matrix is stored
 * @param col index of the column to set
 * @param values values of entries
 */
void JMTX_NAME_TYPED(matrix_drm_set_col)(const JMTX_NAME_TYPED(matrix_drm) * mtx, JMTX_INDEX_T col,
                                         const JMTX_SCALAR_T values[]);

/**
 * Returns the pointers to arrays of column indices and element values for that row
 * @param mtx pointer to the memory where the matrix is stored
 * @param row index of the row to get
 * @param p_elements pointer to array of values
 * @return number of elements in the row, which is the number of valid elements in arrays given to p_indices and
 * p_elements
 */
JMTX_FAST_INT_T JMTX_NAME_TYPED(matrix_drm_get_row)(const JMTX_NAME_TYPED(matrix_drm) * mtx, JMTX_INDEX_T row,
                                                    JMTX_SCALAR_T *p_elements[1]);

/**
 * Multiplies a dense column vector x by the sparse matrix and stores the result at y
 * @param mtx pointer to the memory where the matrix is stored
 * @param x pointer to vector to be multiplied
 * @param y pointer to vector where the result of multiplication is to be stored
 */
void JMTX_NAME_TYPED(matrix_drm_vector_multiply)(const JMTX_NAME_TYPED(matrix_drm) * mtx,
                                                 const JMTX_SCALAR_T *restrict x, JMTX_SCALAR_T *restrict y);

/**
 * Sets a single entry in the matrix.
 * @param mtx pointer to the memory where the matrix is stored
 * @param i row index
 * @param j column index
 * @param value value to which the value is set
 */
void JMTX_NAME_TYPED(matrix_drm_set_entry)(const JMTX_NAME_TYPED(matrix_drm) * mtx, JMTX_INDEX_T i, JMTX_INDEX_T j,
                                           JMTX_SCALAR_T value);

/**
 * Returns a single entry from the matrix.
 * @param mtx pointer to the memory where the matrix is stored
 * @param i row index
 * @param j column index
 * @return value of the entry
 */
JMTX_SCALAR_T JMTX_NAME_TYPED(matrix_drm_get_entry)(const JMTX_NAME_TYPED(matrix_drm) * mtx, JMTX_INDEX_T i,
                                                    JMTX_INDEX_T j);

/**
 * Returns a pointer to an entry in the matrix.
 * @param mtx pointer to the memory where the matrix is stored
 * @param i row index
 * @param j column index
 * @return pointer to an entry
 */
JMTX_SCALAR_T *JMTX_NAME_TYPED(matrix_drm_entry_ptr)(const JMTX_NAME_TYPED(matrix_drm) * mtx, JMTX_INDEX_T i,
                                                     JMTX_INDEX_T j);

/**
 * Zeros all entries within a matrix, but does not remove them in case they need to be reused
 * @param mtx matrix to zero
 */
void JMTX_NAME_TYPED(matrix_drm_zero_all_entries)(const JMTX_NAME_TYPED(matrix_drm) * mtx);

/**
 * Similar to JMTX_NAME_TYPED(matrix_drm_set_all_entries, but slower, since it can not use memset. On the other hand, it
 * allows for the value to be other than 0
 * @param mtx matrix to set
 * @param x value to which to set all entries to
 */
void JMTX_NAME_TYPED(matrix_drm_set_all_entries)(const JMTX_NAME_TYPED(matrix_drm) * mtx, JMTX_SCALAR_T x);

/**
 * Returns the values of entries in the matrix, along with what row of the matrix they were located in
 * @param mtx pointer to the memory where the matrix is stored
 * @param col column index of the matrix to look at
 * @param values a buffer of at least n values which receives the values of the column
 * @return number of entries that were extracted from the column (may be less than are really in the column if n was too
 * small)
 */
JMTX_INDEX_T JMTX_NAME_TYPED(matrix_drm_get_col)(const JMTX_NAME_TYPED(matrix_drm) * mtx, JMTX_INDEX_T col,
                                                 JMTX_SCALAR_T values[]);

/**
 * Creates a transpose of a matrix
 * @param mtx pointer to the memory where the input matrix is stored
 * @param p_out address where the pointer to the output matrix will be returned
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result JMTX_NAME_TYPED(matrix_drm_transpose)(const JMTX_NAME_TYPED(matrix_drm) * mtx,
                                                  JMTX_NAME_TYPED(matrix_drm) * *p_out,
                                                  const jmtx_allocator_callbacks *allocator_callbacks);

/**
 * Creates a transpose of a matrix
 * @param mtx pointer to the memory where the input matrix is stored
 * @param aux_row auxiliary memory to use for intermediate storage
 * @return JMTX_RESULT_SUCCESS if successful
 */
jmtx_result JMTX_NAME_TYPED(matrix_drm_transpose_inplace)(JMTX_NAME_TYPED(matrix_drm) * mtx, JMTX_SCALAR_T *aux_row);

/**
 * Creates a copy of the matrix
 * @param mtx pointer to the memory where the input matrix is stored
 * @param p_out address where the pointer to the output matrix will be returned
 * @param allocator_callbacks pointer to a struct with callbacks and state to use for memory allocation or NULL to use
 * malloc, free, and realloc
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on memory allocation failure
 */
jmtx_result JMTX_NAME_TYPED(matrix_drm_copy)(const JMTX_NAME_TYPED(matrix_drm) * mtx,
                                             JMTX_NAME_TYPED(matrix_drm) * *p_out,
                                             const jmtx_allocator_callbacks *allocator_callbacks);

/**
 * Computes one entry of Ax. This function only computes the i-th entry to make it possible to compute it in parallel.
 * @param mtx pointer to the memory where the matrix A is stored.
 * @param x pointer to the memory where the vector x is stored
 * @param i what entry of the residual to compute
 * @return result of the multiplication
 */
JMTX_SCALAR_T JMTX_NAME_TYPED(matrix_drm_vector_multiply_row)(const JMTX_NAME_TYPED(matrix_drm) * mtx,
                                                              const JMTX_SCALAR_T *x, JMTX_INDEX_T i);

/**
 * Returns the upper bandwidth and the lower bandwidth of the BRM matrix
 * @param mtx pointer to the memory where the matrix is stored
 * @param ubw pointer which receives the upper bandwidth of the matrix
 * @param lbw pointer which receives the lower bandwidth of the matrix
 */
void JMTX_NAME_TYPED(matrix_drm_get_bandwidths)(const JMTX_NAME_TYPED(matrix_drm) * mtx, JMTX_INDEX_T *ubw,
                                                JMTX_INDEX_T *lbw);

/**
 * Exchanges row of matrix through the use of a permutation matrix. This does not cause memory copying and is thus fast
 * @param mtx pointer to the memory where the matrix is stored
 * @param row1 index of a row to exchange with row2
 * @param row2 index of a row to exchange with row1
 */
jmtx_result JMTX_NAME_TYPED(matrix_drm_swap_rows)(JMTX_NAME_TYPED(matrix_drm) * mtx, JMTX_INDEX_T row1,
                                                  JMTX_INDEX_T row2);

/**
 * Permutes all rows at once. The provided list must contain all entries in the set [0, n) exactly once, where n is the
 * number or rows of the matrix.
 * @param mtx pointer to the memory where the matrix is stored
 * @param perm list of permutation indices
 * @return JMTX_RESULT_SUCCESS if successful, JMTX_RESULT_BAD_ALLOC on allocation failure
 */
jmtx_result JMTX_NAME_TYPED(matrix_drm_set_permutation)(JMTX_NAME_TYPED(matrix_drm) * mtx, const JMTX_INDEX_T *perm);

/**
 * Executes permutations of matrix rows and reorders rows in memory. This may allow for better memory access on a
 * permuted matrix. Should not be done on decompositions. This function requires additional memory, but that allows
 * it to swap whole rows at once, which makes it more cache-friendly and potentially faster for large matrices.
 * @param mtx pointer to the memory where the matrix is stored
 */
void JMTX_NAME_TYPED(matrix_drm_commit_permutations)(JMTX_NAME_TYPED(matrix_drm) * mtx);

/**
 * Executes permutations of matrix rows and reorders rows in memory. This may allow for better memory access on a
 * permuted matrix. Should not be done on decompositions. This function requires additional memory, but that allows
 * it to swap whole rows at once, which makes it more cache-friendly and potentially faster for large matrices.
 * @param mtx pointer to the memory where the matrix is stored
 * @param aux_row memory that can be used by the function to store an intermediate matrix row
 */
void JMTX_NAME_TYPED(matrix_drm_commit_permutations2)(JMTX_NAME_TYPED(matrix_drm) * mtx, JMTX_SCALAR_T *aux_row);

/**
 * Computes the result of a Givens rotation being applied on a (square) matrix as multiplication on the left.
 * The rotation is characterized by a rotation angle theta and two indices, which indicate which two rows the rotation
 * is applied to.
 *
 * This function takes cos(theta) and sin(theta) instead of just the angle directly, because in some cases sine and
 * cosine may be computed directly without computing the angle. In that case it would be redundant to convert those into
 * an angle, then convert them back to sine and cosine.
 *
 * @param mtx input matrix, to which the Givens rotation is applied to
 * @param r1 first integer characterizing the rotation
 * @param r2 second integer characterizing the rotation
 * @param ct value of cos(theta), which is used for the rotation
 * @param st value of sin(theta), which is used for the rotation
 */
void JMTX_NAME_TYPED(matrix_drm_givens_rotation_left)(JMTX_NAME_TYPED(matrix_drm) * mtx, unsigned r1, unsigned r2,
                                                      JMTX_SCALAR_T ct, JMTX_SCALAR_T st);

/**
 * Computes the result of a Givens rotation being applied on a (square) matrix as multiplication on the right.
 * The rotation is characterized by a rotation angle theta and two indices, which indicate which two columns the
 * rotation is applied to.
 *
 * This function takes cos(theta) and sin(theta) instead of just the angle directly, because in some cases sine and
 * cosine may be computed directly without computing the angle. In that case it would be redundant to convert those into
 * an angle, then convert them back to sine and cosine.
 *
 * @param mtx input matrix, to which the Givens rotation is applied to
 * @param c1 first integer characterizing the rotation
 * @param c2 second integer characterizing the rotation
 * @param ct value of cos(theta), which is used for the rotation
 * @param st value of sin(theta), which is used for the rotation
 */
void JMTX_NAME_TYPED(matrix_drm_givens_rotation_right)(JMTX_NAME_TYPED(matrix_drm) * mtx, unsigned c1, unsigned c2,
                                                       JMTX_SCALAR_T ct, JMTX_SCALAR_T st);

/**
 * Computes the matrix product of two matrices A and B as C = A B. This is can be written as:
 * $$
 *  C_{i,j} = \sum\limits_{k=0}^{N-1} A_{i,k} B_{k,j}
 * $$
 *
 * This version of the function uses a pre-exsiting matrix as output.
 *
 * @param a matrix A
 * @param b matrix B
 * @param out Where the output should be returned. Should have all permutations committed.
 * @return JMTX_RESULT_SUCCES if successful, JMTX_RESULT_DIMS_MISMATCH if the dimensions of a and b don't allow
 * multiplication, or if c does not have correct dimensions
 */
jmtx_result JMTX_NAME_TYPED(matrix_drm_multiply_matrix)(JMTX_NAME_TYPED(matrix_drm) * a,
                                                        JMTX_NAME_TYPED(matrix_drm) * b,
                                                        JMTX_NAME_TYPED(matrix_drm) * out);

/**
 * Shifts the diagonal of the matrix mtx by the value v, so that: A_{i,i} = A_{i,i} + v for all i
 * @param mtx matrix which should have its diagonal shifted
 * @param v value by which to shift the diagonal
 */
void JMTX_NAME_TYPED(matrix_drm_shift_diagonal)(JMTX_NAME_TYPED(matrix_drm) * mtx, JMTX_SCALAR_T v);

#endif // JMTX_DENSE_ROW_MAJOR_H
