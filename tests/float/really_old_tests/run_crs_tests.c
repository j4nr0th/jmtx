//
// Created by jan on 14.6.2022.
//
#include "../../../include/jmtx/float/matrices/sparse_row_compressed.h"
#include "../test_common.h"
#include <stdio.h>

static int make_the_values_funnier(uint32_t i, uint32_t j, float *x, void *param)
{
    *x += 1.0f;
    return 0;
    uint32_t *p_funny = param;
    switch ((*p_funny) % 3)
    {
    case 0:
        *x += 69.0f;
        break;
    case 1:
        *x += 420.0f;
        break;
    case 2:
        *x += 1337.0f;
        break;
    }
    *p_funny += 1;
    return 0;
}

int main()
{
    jmtxf_matrix_crs *matrix;
    jmtxf_matrix_crs_new(&matrix, 25, 25, 2, NULL);

    const uint32_t indices[3] = {1, 2, 4};
    const float values[3] = {69, 420, 505};
    jmtxf_matrix_crs_set_row(matrix, 2, 3, indices, values);

    const uint32_t indices2[3] = {0, 2};
    const float values2[3] = {1, 20};
    jmtxf_matrix_crs_set_row(matrix, 4, 2, indices2, values2);

    print_crs_matrix(matrix);
    jmtxf_matrix_crs_remove_bellow_magnitude(matrix, 5.0f);
    const uint32_t indices3[5] = {0, 10, 13, 15, 23};
    const float values3[5] = {-1.0f, 1.0f, -2.0f, 1.0f, 20};
    jmtxf_matrix_crs_set_row(matrix, 14, 5, indices3, values3);
    print_crs_matrix(matrix);
    jmtxf_matrix_crs_shrink(matrix);

    uint32_t funny = 0;
    jmtxf_matrix_crs_set_entry(matrix, 1, 9, -100);
    jmtxf_matrix_crs_set_entry(matrix, 6, 9, -200);
    jmtxf_matrix_crs_set_entry(matrix, 2, 24, -300);
    jmtxf_matrix_crs_apply_unary_fn(matrix, make_the_values_funnier, &funny);
    jmtxf_matrix_crs_set_entry(matrix, 2, 9, -400);
    jmtxf_matrix_crs_shrink(matrix);
    jmtxf_matrix_crs_set_entry(matrix, 2, 1, -400);
    jmtxf_matrix_crs_set_entry(matrix, 20, 1, -400);
    printf("\n\nBREAK 1\n");
    print_crs_matrix(matrix);
    printf("Row 2 is:\n");
    {
        uint32_t n;
        uint32_t *ind;
        float *val;
        n = jmtxf_matrix_crs_get_row(matrix, 2, &ind, &val);
        for (uint32_t i = 0; i < n; ++i)
        {
            printf(" (%u, %g)", ind[i], val[i]);
        }
        printf("\n");
    }
    jmtxf_matrix_crs_remove_zeros(matrix);

    print_crs_matrix(matrix);
    jmtxf_matrix_crs *matrix_1;
    jmtxf_matrix_crs *matrix_2;

    jmtxf_matrix_crs_new(&matrix_1, 5, 5, 0, NULL);

    {
        const uint32_t i[] = {0};
        const float v[] = {1.0f};
        jmtxf_matrix_crs_set_row(matrix_1, 0, 1, i, v);
    }

    {
        const uint32_t i[] = {0, 1};
        const float v[] = {0.75f, 0.25f};
        jmtxf_matrix_crs_set_row(matrix_1, 1, 2, i, v);
    }

    {
        const uint32_t i[] = {1, 2, 3};
        const float v[] = {1, 2, 1};
        jmtxf_matrix_crs_set_row(matrix_1, 2, 3, i, v);
    }

    {
        const uint32_t i[] = {3, 4};
        const float v[] = {0.75f, 0.25f};
        jmtxf_matrix_crs_set_row(matrix_1, 3, 2, i, v);
    }

    {
        const uint32_t i[] = {4};
        const float v[] = {1};
        jmtxf_matrix_crs_set_row(matrix_1, 4, 1, i, v);
    }
    jmtxf_matrix_crs_transpose(matrix_1, &matrix_2, NULL);

    printf("\n\nBREAK 2\n");
    print_crs_matrix(matrix_1);
    print_crs_matrix(matrix_2);

    uint32_t err_count = 0;

    jmtxf_matrix_crs *tmp;
    jmtxf_matrix_crs_transpose(matrix, &tmp, NULL);
    for (uint32_t i = 0; i < 25; ++i)
    {
        for (uint32_t j = 0; j < 25; ++j)
        {
            float x1, x2;
            x1 = jmtxf_matrix_crs_get_entry(matrix, i, j);
            x2 = jmtxf_matrix_crs_get_entry(tmp, j, i);
            if (x1 != x2)
            {
                err_count += 1;
            }
        }
    }
    printf("Error count for transpose: %u\n", err_count);
    jmtxf_matrix_crs_destroy(tmp);

    jmtxf_matrix_crs_copy(matrix, &tmp, NULL);
    err_count = 0;
    for (uint32_t i = 0; i < 25; ++i)
    {
        for (uint32_t j = 0; j < 25; ++j)
        {
            float x1, x2;
            x1 = jmtxf_matrix_crs_get_entry(matrix, i, j);
            x2 = jmtxf_matrix_crs_get_entry(tmp, i, j);
            if (x1 != x2)
            {
                err_count += 1;
            }
        }
    }
    printf("Error count for copy: %u\n", err_count);

    jmtxf_matrix_crs_destroy(tmp);

    jmtxf_matrix_crs_destroy(matrix_1);
    jmtxf_matrix_crs_destroy(matrix_2);

    const float v1[25] = {1, 2, 3, 4, 5, [10] = 69.0f};
    float y[25];

    jmtxf_matrix_crs_vector_multiply(matrix, v1, y);
    printf("x:");
    for (uint32_t i = 0; i < 25; ++i)
    {
        printf(" %g", v1[i]);
    }
    printf("\n");
    printf("y:");
    for (uint32_t i = 0; i < 25; ++i)
    {
        printf(" %g", y[i]);
    }

    printf("\nHello world!\n");
    //    printf("Matrix is currently using %zu bytes of memory\n", jmtxf_matrix_crs_memory_usage(matrix));
    jmtxf_matrix_crs_destroy(matrix);
    printf("Long double size: %zu bytes\n", sizeof(long double));
    return 0;
}
