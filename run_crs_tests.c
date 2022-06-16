//
// Created by jan on 14.6.2022.
//
#include "sparse_row_compressed.h"
#include <stdio.h>
#include <assert.h>

static void print_matrix(const CrsMatrix* mtx)
{
    printf("elements:");
    for (uint i = 0, l = 0; i < mtx->n_elements; ++i)
    {
        if (l < mtx->rows && mtx->elements_before[l + 1] <= i + 1)
        {
            l += 1;
        }
        printf(" (%u, %u, %+3.3f)", l, mtx->indices[i + 1], mtx->elements[i + 1]);
    }
    printf("\nelement offsets:");
    for (uint i = 0; i < mtx->rows; ++i)
    {
        printf(" %u,", mtx->elements_before[i]);
    }
    printf(" %u", mtx->elements_before[mtx->rows]);
    printf("\nMatrix:\n[");

    for (uint i = 0; i < mtx->rows; ++i)
    {
        if (i != 0)
            printf(" [");
        else
            printf("[");
        scalar_t x;
        for (uint j = 0; j < mtx->columns - 1; ++j)
        {
            matrix_crs_get_element(mtx, i, j, &x);
            printf("%+03.3f ", x);
        }
        matrix_crs_get_element(mtx, i, mtx->columns - 1, &x);
        printf("%+03.3f] - %u", x, mtx->elements_before[i + 1] - mtx->elements_before[i]);
        if (i != mtx->rows - 1) printf("\n");
    }

    printf("]\n");
}

#define BEEF_CHECK(mtx) matrix_crs_beef_check(&(mtx), &beef_status), assert(beef_status == 0xBEEF)

static int make_the_values_funnier(uint i, uint j, scalar_t* x, void* param)
{
    *x += 1.0f;
    return 0;
    uint* p_funny = param;
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
    CrsMatrix matrix;
    matrix_crs_new(&matrix, 25, 25, 2);
    BEEF_CHECK(matrix);

    const uint indices[3] = {1, 2, 4};
    const scalar_t values[3] = { 69, 420, 505};
    matrix_crs_set_row(&matrix, 2, 3, indices, values);
    BEEF_CHECK(matrix);
    const uint indices2[3] = {0, 2};
    const scalar_t values2[3] = { 1, 20};
    matrix_crs_set_row(&matrix, 4, 2, indices2, values2);
    BEEF_CHECK(matrix);
    print_matrix(&matrix);
    matrix_crs_remove_bellow(&matrix, 5.0f);
    const uint indices3[5] = {0, 10, 13, 15, 23};
    const scalar_t values3[5] = {-1.0f, 1.0f, -2.0f, 1.0f, 20};
    matrix_crs_set_row(&matrix, 14, 5, indices3, values3);
    print_matrix(&matrix);
    matrix_crs_shrink(&matrix);
    BEEF_CHECK(matrix);


    uint funny = 0;
    matrix_crs_set_element(&matrix, 1, 9,  -100);
    matrix_crs_set_element(&matrix, 6, 9,  -200);
    matrix_crs_set_element(&matrix, 2, 24, -300);
    matrix_crs_apply_unary_fn(&matrix, make_the_values_funnier, &funny);
    matrix_crs_set_element(&matrix, 2, 9, -400);
    matrix_crs_shrink(&matrix);
    matrix_crs_set_element(&matrix, 2, 1, -400);
    matrix_crs_set_element(&matrix, 20, 1, -400);
    printf("\n\nBREAK 1\n");
    print_matrix(&matrix);
    printf("Row 2 is:\n");
    {
        uint n;
        uint* ind;
        scalar_t* val;
        matrix_crs_get_row(&matrix, 2, &n, &ind, &val);
        for (uint i = 0; i < n; ++i)
        {
            printf(" (%u, %g)", ind[i], val[i]);
        }
        printf("\n");
    }
    matrix_crs_remove_zeros(&matrix);
    BEEF_CHECK(matrix);
    print_matrix(&matrix);
    CrsMatrix matrix_1;
    CrsMatrix matrix_2;

    matrix_crs_new(&matrix_1, 5, 5, 0);

    {
        const uint i[] = { 0 };
        const scalar_t v[] = { 1.0f };
        matrix_crs_set_row(&matrix_1, 0, 1, i, v);
    }

    {
        const uint i[] = { 0, 1 };
        const scalar_t v[] = { 0.75f, 0.25f };
        matrix_crs_set_row(&matrix_1, 1, 2, i, v);
    }

    {
        const uint i[] = { 1, 2, 3 };
        const scalar_t v[] = { 1, 2, 1};
        matrix_crs_set_row(&matrix_1, 2, 3, i, v);
    }

    {
        const uint i[] = { 3, 4 };
        const scalar_t v[] = { 0.75f, 0.25f };
        matrix_crs_set_row(&matrix_1, 3, 2, i, v);
    }

    {
        const uint i[] = { 4 };
        const scalar_t v[] = {1};
        matrix_crs_set_row(&matrix_1, 4, 1, i, v);
    }
    matrix_crs_transpose(&matrix_1, &matrix_2);

    printf("\n\nBREAK 2\n");
    print_matrix(&matrix_1);
    print_matrix(&matrix_2);

    uint err_count = 0;

    CrsMatrix tmp;
    matrix_crs_transpose(&matrix, &tmp);
    for (uint i = 0; i < 25; ++i)
    {
        for (uint j = 0; j < 25; ++j)
        {
            scalar_t x1,x2;
            matrix_crs_get_element(&matrix, i, j, &x1);
            matrix_crs_get_element(&tmp, j, i, &x2);
            if (x1 != x2)
            {
                err_count += 1;
            }
        }
    }
    printf("Error count for transpose: %u\n", err_count);
    matrix_crs_destroy(&tmp);

    matrix_crs_copy(&matrix, &tmp);
    err_count = 0;
    for (uint i = 0; i < 25; ++i)
    {
        for (uint j = 0; j < 25; ++j)
        {
            scalar_t x1,x2;
            matrix_crs_get_element(&matrix, i, j, &x1);
            matrix_crs_get_element(&tmp, i, j, &x2);
            if (x1 != x2)
            {
                err_count += 1;
            }
        }
    }
    printf("Error count for copy: %u\n", err_count);

    matrix_crs_destroy(&tmp);


    matrix_crs_destroy(&matrix_1);
    matrix_crs_destroy(&matrix_2);



    const scalar_t v1[25] = {1, 2, 3, 4, 5, [10] = 69.0f};
    scalar_t y[25];

    matrix_crs_vector_multiply(&matrix, v1, y);
    printf("x:");
    for (uint i = 0; i < 25; ++i)
    {
        printf(" %g", v1[i]);
    }
    printf("\n");
    printf("y:");
    for (uint i = 0; i < 25; ++i)
    {
        printf(" %g", y[i]);
    }

    BEEF_CHECK(matrix);
    printf("\nHello world!\n");
    printf("Matrix is currently using %zu bytes of memory\n", matrix_crs_memory_usage(matrix));
    matrix_crs_destroy(&matrix);
    printf("Beef status: %x\n", beef_status);
    printf("Long double size: %zu bytes\n", sizeof(long double));
    return 0;
}

