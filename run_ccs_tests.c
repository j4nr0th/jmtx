//
// Created by jan on 15.6.2022.
//
#include "source/matrices/sparse_column_compressed.h"
#include <stdio.h>
#include <assert.h>

static void print_matrix(const jmtx_matrix_ccs* mtx)
{
    printf("elements:");
    for (uint32_t i = 0, l = 0; i < mtx->n_elements; ++i)
    {
        if (l < mtx->base.rows && mtx->elements_before[l + 1] <= i + 1)
        {
            l += 1;
        }
        printf(" (%u, %u, %g)", l, mtx->indices[i + 1], mtx->elements[i + 1]);
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
        jmtx_scalar_t x;
        for (uint32_t j = 0; j < mtx->base.cols - 1; ++j)
        {
            matrix_ccs_get_element(mtx, i, j, &x);
            printf("%g ", x);
        }
        matrix_ccs_get_element(mtx, i, mtx->base.cols - 1, &x);
        printf("%g] - %u", x, mtx->elements_before[i + 1] - mtx->elements_before[i]);
        if (i != mtx->base.rows - 1) printf("\n");
    }

    printf("]\n");
}

#define BEEF_CHECK(mtx) matrix_ccs_beef_check(&(mtx), &beef_status), assert(beef_status == 0xBEEF)

static int make_the_values_funnier(uint32_t i, uint32_t j, jmtx_scalar_t* x, void* param)
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
    jmtx_matrix_ccs matrix;
    matrix_ccs_new(&matrix, 25, 25, 2, NULL);
    BEEF_CHECK(matrix);

    const uint32_t indices[3] = {1, 2, 4};
    const jmtx_scalar_t values[3] = { 69, 420, 505};
    matrix_ccs_set_col(&matrix, 2, 3, indices, values);
    BEEF_CHECK(matrix);
    const uint32_t indices2[3] = {0, 2};
    const jmtx_scalar_t values2[3] = { 1, 20};
    matrix_ccs_set_col(&matrix, 4, 2, indices2, values2);
    BEEF_CHECK(matrix);
    print_matrix(&matrix);
    matrix_ccs_remove_bellow(&matrix, 5.0f);
    const uint32_t indices3[5] = {0, 10, 13, 15, 23};
    const jmtx_scalar_t values3[5] = {-1.0f, 1.0f, -2.0f, 1.0f, 20};
    matrix_ccs_set_col(&matrix, 14, 5, indices3, values3);
    print_matrix(&matrix);
    matrix_ccs_shrink(&matrix);
    BEEF_CHECK(matrix);


    uint32_t funny = 0;
    matrix_ccs_set_element(&matrix, 1, 9,  -100);
    matrix_ccs_set_element(&matrix, 6, 9,  -200);
    matrix_ccs_set_element(&matrix, 2, 24, -300);
    matrix_ccs_apply_unary_fn(&matrix, make_the_values_funnier, &funny);
    matrix_ccs_set_element(&matrix, 2, 9, -400);
    matrix_ccs_shrink(&matrix);
    matrix_ccs_set_element(&matrix, 2, 1, -400);
    matrix_ccs_set_element(&matrix, 20, 1, -400);
    printf("\n\nBREAK 1\n");
    print_matrix(&matrix);
    printf("Column 2 is:\n");
    {
        uint32_t n;
        uint32_t* ind;
        jmtx_scalar_t* val;
        matrix_ccs_get_col(&matrix, 2, &n, &ind, &val);
        for (uint32_t i = 0; i < n; ++i)
        {
            printf(" (%u, %g)", ind[i], val[i]);
        }
        printf("\n");
    }
    matrix_ccs_remove_zeros(&matrix);
    BEEF_CHECK(matrix);
    print_matrix(&matrix);
    jmtx_matrix_ccs matrix_1;
    jmtx_matrix_ccs matrix_2;

    matrix_ccs_new(&matrix_1, 5, 5, 0, NULL);

    {
        const uint32_t i[] = { 0 };
        const jmtx_scalar_t v[] = { 1.0f };
        matrix_ccs_set_col(&matrix_1, 0, 1, i, v);
    }

    {
        const uint32_t i[] = { 0, 1 };
        const jmtx_scalar_t v[] = { 0.75f, 0.25f };
        matrix_ccs_set_col(&matrix_1, 1, 2, i, v);
    }

    {
        const uint32_t i[] = { 1, 2, 3 };
        const jmtx_scalar_t v[] = { 1, 2, 1};
        matrix_ccs_set_col(&matrix_1, 2, 3, i, v);
    }

    {
        const uint32_t i[] = { 3, 4 };
        const jmtx_scalar_t v[] = { 0.75f, 0.25f };
        matrix_ccs_set_col(&matrix_1, 3, 2, i, v);
    }

    {
        const uint32_t i[] = { 4 };
        const jmtx_scalar_t v[] = {1};
        matrix_ccs_set_col(&matrix_1, 4, 1, i, v);
    }
    matrix_ccs_transpose(&matrix_1, &matrix_2);

    printf("\n\nBREAK 2\n");
    print_matrix(&matrix_1);
    print_matrix(&matrix_2);

    uint32_t err_count = 0;

    jmtx_matrix_ccs tmp;
    matrix_ccs_transpose(&matrix, &tmp);
    print_matrix(&matrix);
    print_matrix(&tmp);
    for (uint32_t i = 0; i < 25; ++i)
    {
        for (uint32_t j = 0; j < 25; ++j)
        {
            jmtx_scalar_t x1,x2;
            matrix_ccs_get_element(&matrix, i, j, &x1);
            matrix_ccs_get_element(&tmp, j, i, &x2);
            if (x1 != x2)
            {
                err_count += 1;
                matrix_ccs_get_element(&matrix, i, j, &x1);
                matrix_ccs_get_element(&tmp, j, i, &x2);
            }
        }
    }
    printf("Error count for transpose: %u\n", err_count);
    matrix_ccs_destroy(&tmp);

    matrix_ccs_copy(&matrix, &tmp);
    err_count = 0;
    for (uint32_t i = 0; i < 25; ++i)
    {
        for (uint32_t j = 0; j < 25; ++j)
        {
            jmtx_scalar_t x1,x2;
            matrix_ccs_get_element(&matrix, i, j, &x1);
            matrix_ccs_get_element(&tmp, i, j, &x2);
            if (x1 != x2)
            {
                err_count += 1;
            }
        }
    }
    printf("Error count for copy: %u\n", err_count);

    matrix_ccs_destroy(&tmp);


    matrix_ccs_destroy(&matrix_1);
    matrix_ccs_destroy(&matrix_2);



    const jmtx_scalar_t v1[25] = {1, 2, 3, 4, 5, [10] = 69.0f};
    jmtx_scalar_t y[25];

    matrix_ccs_vector_multiply(&matrix, v1, y);
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
//    printf("Matrix is currently using %zu bytes of memory\n", matrix_ccs_memory_usage(matrix));
    matrix_ccs_destroy(&matrix);
    printf("Beef status: %x\n", beef_status);
    printf("Long double size: %zu bytes\n", sizeof(long double));
    return 0;
}

