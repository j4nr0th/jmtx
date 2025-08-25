//
// Created by jan on 21.7.2023.
//

#ifndef JMTXZ_DENSE_ROW_MAJOR_INTERNAL_H
#define JMTXZ_DENSE_ROW_MAJOR_INTERNAL_H
#ifndef JMTX_MATRIX_BASE_INTERNAL_H
#    include "../../matrix_base_internal.h"
#endif
#ifndef JMTXC_DENSE_ROW_MAJOR_H
#    include "../../../include/jmtx/cdouble/matrices/dense_row_major.h"
#endif
struct jmtxz_matrix_drm_struct
{
    jmtx_matrix_base base;
    //  Length this->base.rows, contains the order in which the rows are permuted (what row is read from where)
    uint32_t *permutations;
    //  Length this->base.rows, contains the reverse order in which the rows are permuted (what row is written to where)
    uint32_t *rperm;
    //  Length this->base.rows * this->base.cols, contains every entry in row-major ordering.
    _Complex double *restrict values;
};

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
jmtxz_matrix_drm jmtxz_matrix_drm_from_data(unsigned rows, unsigned cols,
                                            _Complex double values[JMTX_ARRAY_ATTRIB(static rows * cols)]);

#endif // JMTXZ_DENSE_ROW_MAJOR_INTERNAL_H
