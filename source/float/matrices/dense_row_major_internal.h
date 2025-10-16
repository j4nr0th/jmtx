//
// Created by jan on 21.7.2023.
//

#ifndef JMTX_DENSE_ROW_MAJOR_INTERNAL_H
#define JMTX_DENSE_ROW_MAJOR_INTERNAL_H
#ifndef JMTX_MATRIX_BASE_INTERNAL_H
#    include "../../matrix_base_internal.h"
#endif
#ifndef JMTXF_DENSE_ROW_MAJOR_H
#    include "../../../include/jmtx/float/matrices/dense_row_major.h"
#endif
struct jmtxf_matrix_drm_struct
{
    jmtx_matrix_base base;
    //  Length this->base.rows, contains the order in which the rows are permuted (what row is read from where)
    uint32_t *permutations;
    //  Length this->base.rows, contains the reverse order in which the rows are permuted (what row is written to where)
    uint32_t *rperm;
    //  Length this->base.rows * this->base.cols, contains every entry in row-major ordering.
    float *restrict values;
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
jmtxf_matrix_drm jmtxf_matrix_drm_from_data(unsigned rows, unsigned cols,
                                            float values[JMTX_ARRAY_ATTRIB(static rows * cols)]);

jmtx_result jmtxf_matrix_drm_set_default_permutations(jmtxf_matrix_drm *this);

inline void jmtxf_matrix_drm_clear_permutations(jmtxf_matrix_drm *this)
{
    if (this->permutations)
    {
        this->base.allocator_callbacks.free(this->base.allocator_callbacks.state, this->permutations);
        this->permutations = NULL;
    }
    if (this->rperm)
    {
        this->base.allocator_callbacks.free(this->base.allocator_callbacks.state, this->rperm);
        this->rperm = NULL;
    }
}

inline float *MATRIX_DRM_ENTRY_READ(const jmtxf_matrix_drm *this, unsigned row, unsigned col)
{
    if (this->permutations)
    {
        return this->values + this->permutations[row] * this->base.cols + col;
    }
    return this->values + row * this->base.cols + col;
}

inline float *MATRIX_DRM_ENTRY_WRITE(const jmtxf_matrix_drm *this, unsigned row, unsigned col)
{
    if (this->rperm)
    {
        return this->values + this->rperm[row] * this->base.cols + col;
    }
    return this->values + row * this->base.cols + col;
}

#endif // JMTX_DENSE_ROW_MAJOR_INTERNAL_H
