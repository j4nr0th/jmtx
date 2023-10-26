//
// Created by jan on 29.9.2023.
//

#ifndef JMTX_F_INTERFACE_H
#define JMTX_F_INTERFACE_H
#include <stdint.h>
#include <ISO_Fortran_binding.h>

int32_t f_bin_search(int32_t v, int32_t len, const int32_t* values);

float f_crs_mul_vec_row(int32_t n_elements, int32_t dim, const int32_t* indices, const float* values, const float* vec);

int32_t f_crs_mul_vec(int32_t n_elements, int32_t dims, const int32_t* elements_before, const int32_t* indices, const float* values, const float* vector_in, float* vector_out);

void f_pt_jacobi(int32_t n_elements, int32_t dims, const int32_t* end_of_row_offsets, const int32_t* indices,
                 const float* values, float tolerance, int32_t max_iterations, float* x, const float* y);

#endif //JMTX_F_INTERFACE_H
