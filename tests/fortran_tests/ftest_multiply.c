//
// Created by jan on 29.9.2023.
//
#include <stdio.h>
#include "../test_common.h"
#include "../../source/f_exper/f_interface.h"
#include "../../source/matrices/sparse_row_compressed.h"
#include "../../source/matrices/sparse_row_compressed_internal.h"
#include "../../source/matrices/sparse_row_compressed_safe.h"

enum
{
    ARRAY_LEN = 32,
    MTX_DIMS = 4,
    N_MTX_ELM = 6,
};

int main()
{
    int32_t array_of_values[ARRAY_LEN];
    for (int32_t i = 0; i < ARRAY_LEN; ++i)
    {
        array_of_values[i] = i;
    }

    const int32_t len = ARRAY_LEN;
    int32_t ret_v;
    for (int32_t i = 0; i < ARRAY_LEN; ++i)
    {
        int32_t v_i = i;
        ASSERT((ret_v = f_bin_search(v_i, len, array_of_values) - 1) == i);
    }
    for (int32_t i = 0; i < ARRAY_LEN; ++i)
    {
        array_of_values[i] = 2 * i + 3;
    }
    for (int32_t i = 0; i < ARRAY_LEN; ++i)
    {
        int32_t v_i = 2 * i + 3;
        ASSERT((ret_v = f_bin_search(v_i, len, array_of_values) - 1) == i);
    }
    for (int32_t i = 0; i < ARRAY_LEN; ++i)
    {
        int32_t v_i = 2 * i;
        ASSERT((ret_v = f_bin_search(v_i, len, array_of_values) - 1) == -1);
    }
    
    jmtx_matrix_crs* mtx;
    jmtx_result mtx_res;

    MATRIX_TEST_CALL(jmtxs_matrix_crs_new(&mtx, MTX_DIMS, MTX_DIMS, N_MTX_ELM, NULL));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);

    const uint32_t col_indices[N_MTX_ELM] =
            {
                0,    2,
                0,
                   1, 2,
                         3,
            };
    const uint32_t row_indices[N_MTX_ELM] =
            {
                0,    0,
                1,
                   2, 2,
                         3,
            };
    const float values[N_MTX_ELM] =
            {
                0.3f, 0.1f,
                4.0f,
                1.0f, -2.0f,
                4.2f
            };

    const float vec_in[MTX_DIMS] = {1.0f, 2.0f, 3.0f, 4.0f};

    float vec_out[MTX_DIMS];

    for (unsigned i = 0; i < N_MTX_ELM; ++i)
    {
        MATRIX_TEST_CALL(jmtxs_matrix_crs_set_entry(mtx, row_indices[i], col_indices[i], values[i]));
        ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    }

    print_crs_matrix(mtx);

    const int a = f_crs_mul_vec(
            (int32_t)mtx->n_entries,
            (int32_t)mtx->base.rows,
            (const int32_t*)mtx->end_of_row_offsets,
            (const int32_t*)mtx->indices,
            (const float*)mtx->values,
            vec_in,
            vec_out
            );
    ASSERT(a == 1);
    for (unsigned i = 0; i < MTX_DIMS; ++i)
    {
        ASSERT(vec_out[i] == jmtx_matrix_crs_vector_multiply_row(mtx, vec_in, i));
    }


    jmtxs_matrix_crs_destroy(mtx);
    return 0;
}

