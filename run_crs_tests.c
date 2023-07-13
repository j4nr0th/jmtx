//
// Created by jan on 14.6.2022.
//
#include "source/matrices/sparse_row_compressed.h"
#include <stdio.h>
#include <assert.h>

static void print_matrix(const jmtx_matrix_crs* mtx)
{
    printf("elements:");
    for (uint32_t i = 0, l = 0; i < mtx->n_elements; ++i)
    {
        if (l < mtx->base.rows && mtx->elements_before[l + 1] <= i + 1)
        {
            l += 1;
        }
        printf(" (%u, %u, %+3.3f)", l, mtx->indices[i + 1], mtx->elements[i + 1]);
    }
    printf("\nelement offsets:");
    for (uint32_t i = 0; i < mtx->base.rows; ++i)
    {
        printf(" %u,", mtx->elements_before[i]);
    }
    printf(" %u", mtx->elements_before[mtx->base.rows]);
    printf("\nMatrix:\n[");

    for (uint32_t i = 0; i < mtx->base.rows; ++i)
    {
        if (i != 0)
            printf(" [");
        else
            printf("[");
        float x;
        for (uint32_t j = 0; j < mtx->base.cols - 1; ++j)
        {
            jmtx_matrix_crs_get_element(mtx, i, j, &x);
            printf("%+03.3f ", x);
        }
        jmtx_matrix_crs_get_element(mtx, i, mtx->base.cols - 1, &x);
        printf("%+03.3f] - %u", x, mtx->elements_before[i + 1] - mtx->elements_before[i]);
        if (i != mtx->base.rows - 1) printf("\n");
    }

    printf("]\n");
}

#define BEEF_CHECK(mtx) jmtx_matrix_crs_beef_check(&(mtx), &beef_status), assert(beef_status == 0xBEEF)

static int make_the_values_funnier(uint32_t i, uint32_t j, float* x, void* param)
{
    *x += 1.0f;
    return 0;
    uint32_t* p_funny = param;
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
    int beef_status;
    jmtx_matrix_crs matrix;
    jmtx_matrix_crs_new(&matrix, 25, 25, 2, NULL);
    BEEF_CHECK(matrix);

    const uint32_t indices[3] = {1, 2, 4};
    const float values[3] = { 69, 420, 505};
    jmtx_matrix_crs_set_row(&matrix, 2, 3, indices, values);
    BEEF_CHECK(matrix);
    const uint32_t indices2[3] = {0, 2};
    const float values2[3] = { 1, 20};
    jmtx_matrix_crs_set_row(&matrix, 4, 2, indices2, values2);
    BEEF_CHECK(matrix);
    print_matrix(&matrix);
    jmtx_matrix_crs_remove_bellow(&matrix, 5.0f);
    const uint32_t indices3[5] = {0, 10, 13, 15, 23};
    const float values3[5] = {-1.0f, 1.0f, -2.0f, 1.0f, 20};
    jmtx_matrix_crs_set_row(&matrix, 14, 5, indices3, values3);
    print_matrix(&matrix);
    jmtx_matrix_crs_shrink(&matrix);
    BEEF_CHECK(matrix);


    uint32_t funny = 0;
    jmtx_matrix_crs_set_element(&matrix, 1, 9, -100);
    jmtx_matrix_crs_set_element(&matrix, 6, 9, -200);
    jmtx_matrix_crs_set_element(&matrix, 2, 24, -300);
    jmtx_matrix_crs_apply_unary_fn(&matrix, make_the_values_funnier, &funny);
    jmtx_matrix_crs_set_element(&matrix, 2, 9, -400);
    jmtx_matrix_crs_shrink(&matrix);
    jmtx_matrix_crs_set_element(&matrix, 2, 1, -400);
    jmtx_matrix_crs_set_element(&matrix, 20, 1, -400);
    printf("\n\nBREAK 1\n");
    print_matrix(&matrix);
    printf("Row 2 is:\n");
    {
        uint32_t n;
        uint32_t* ind;
        float* val;
        jmtx_matrix_crs_get_row(&matrix, 2, &n, &ind, &val);
        for (uint32_t i = 0; i < n; ++i)
        {
            printf(" (%u, %g)", ind[i], val[i]);
        }
        printf("\n");
    }
    jmtx_matrix_crs_remove_zeros(&matrix);
    BEEF_CHECK(matrix);
    print_matrix(&matrix);
    jmtx_matrix_crs matrix_1;
    jmtx_matrix_crs matrix_2;

    jmtx_matrix_crs_new(&matrix_1, 5, 5, 0, NULL);

    {
        const uint32_t i[] = { 0 };
        const float v[] = { 1.0f };
        jmtx_matrix_crs_set_row(&matrix_1, 0, 1, i, v);
    }

    {
        const uint32_t i[] = { 0, 1 };
        const float v[] = { 0.75f, 0.25f };
        jmtx_matrix_crs_set_row(&matrix_1, 1, 2, i, v);
    }

    {
        const uint32_t i[] = { 1, 2, 3 };
        const float v[] = { 1, 2, 1};
        jmtx_matrix_crs_set_row(&matrix_1, 2, 3, i, v);
    }

    {
        const uint32_t i[] = { 3, 4 };
        const float v[] = { 0.75f, 0.25f };
        jmtx_matrix_crs_set_row(&matrix_1, 3, 2, i, v);
    }

    {
        const uint32_t i[] = { 4 };
        const float v[] = {1};
        jmtx_matrix_crs_set_row(&matrix_1, 4, 1, i, v);
    }
    jmtx_matrix_crs_transpose(&matrix_1, &matrix_2);

    printf("\n\nBREAK 2\n");
    print_matrix(&matrix_1);
    print_matrix(&matrix_2);

    uint32_t err_count = 0;

    jmtx_matrix_crs tmp;
    jmtx_matrix_crs_transpose(&matrix, &tmp);
    for (uint32_t i = 0; i < 25; ++i)
    {
        for (uint32_t j = 0; j < 25; ++j)
        {
            float x1,x2;
            jmtx_matrix_crs_get_element(&matrix, i, j, &x1);
            jmtx_matrix_crs_get_element(&tmp, j, i, &x2);
            if (x1 != x2)
            {
                err_count += 1;
            }
        }
    }
    printf("Error count for transpose: %u\n", err_count);
    jmtx_matrix_crs_destroy(&tmp);

    jmtx_matrix_crs_copy(&matrix, &tmp);
    err_count = 0;
    for (uint32_t i = 0; i < 25; ++i)
    {
        for (uint32_t j = 0; j < 25; ++j)
        {
            float x1,x2;
            jmtx_matrix_crs_get_element(&matrix, i, j, &x1);
            jmtx_matrix_crs_get_element(&tmp, i, j, &x2);
            if (x1 != x2)
            {
                err_count += 1;
            }
        }
    }
    printf("Error count for copy: %u\n", err_count);

    jmtx_matrix_crs_destroy(&tmp);


    jmtx_matrix_crs_destroy(&matrix_1);
    jmtx_matrix_crs_destroy(&matrix_2);



    const float v1[25] = {1, 2, 3, 4, 5, [10] = 69.0f};
    float y[25];

    jmtx_matrix_crs_vector_multiply(&matrix, v1, y);
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

    BEEF_CHECK(matrix);
    printf("\nHello world!\n");
//    printf("Matrix is currently using %zu bytes of memory\n", jmtx_matrix_crs_memory_usage(matrix));
    jmtx_matrix_crs_destroy(&matrix);
    printf("Beef status: %x\n", beef_status);
    printf("Long double size: %zu bytes\n", sizeof(long double));
    return 0;
}

