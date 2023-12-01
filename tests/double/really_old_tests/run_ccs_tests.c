// Automatically generated from tests/float/really_old_tests/run_ccs_tests.c on Fri Dec  1 06:43:09 2023
//
// Created by jan on 15.6.2022.
//
#include "../../../include/jmtx/double/matrices/sparse_column_compressed.h"
#include "../../../source/double/matrices/sparse_column_compressed_internal.h"
#include "../test_common.h"
#include <stdio.h>
#include <assert.h>

#include "../../../include/jmtx/double/matrices/sparse_column_compressed_safe.h"

static int make_the_values_funnier(uint32_t i, uint32_t j, double* x, void* param)
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
    jmtxd_matrix_ccs* matrix;
    jmtxds_matrix_ccs_new(&matrix, 25, 25, 2, NULL);

    const uint32_t indices[3] = {1, 2, 4};
    const double values[3] = { 69, 420, 505};
    jmtxds_matrix_ccs_set_col(matrix, 2, 3, indices, values);
    const uint32_t indices2[3] = {0, 2};
    const double values2[3] = { 1, 20};
    jmtxds_matrix_ccs_set_col(matrix, 4, 2, indices2, values2);
    print_ccsd_matrix(matrix);
    jmtxds_matrix_ccs_remove_bellow(matrix, 5.0f);
    const uint32_t indices3[5] = {0, 10, 13, 15, 23};
    const double values3[5] = {-1.0f, 1.0f, -2.0f, 1.0f, 20};
    jmtxds_matrix_ccs_set_col(matrix, 14, 5, indices3, values3);
    print_ccsd_matrix(matrix);
    jmtxds_matrix_ccs_shrink(matrix);

    uint32_t funny = 0;
    jmtxds_matrix_ccs_set_entry(matrix, 1, 9, -100);
    jmtxds_matrix_ccs_set_entry(matrix, 6, 9, -200);
    jmtxds_matrix_ccs_set_entry(matrix, 2, 24, -300);
    jmtxds_matrix_ccs_apply_unary_fn(matrix, make_the_values_funnier, &funny);
    jmtxds_matrix_ccs_set_entry(matrix, 2, 9, -400);
    jmtxds_matrix_ccs_shrink(matrix);
    jmtxds_matrix_ccs_set_entry(matrix, 2, 1, -400);
    jmtxds_matrix_ccs_set_entry(matrix, 20, 1, -400);
    printf("\n\nBREAK 1\n");
    print_ccsd_matrix(matrix);
    printf("Column 2 is:\n");
    {
        uint32_t n;
        uint32_t* ind;
        double* val;
        jmtxds_matrix_ccs_get_col(matrix, 2, &n, &ind, &val);
        for (uint32_t i = 0; i < n; ++i)
        {
            printf(" (%u, %g)", ind[i], val[i]);
        }
        printf("\n");
    }
    jmtxds_matrix_ccs_remove_zeros(matrix);
    print_ccsd_matrix(matrix);
    jmtxd_matrix_ccs* matrix_1;
    jmtxd_matrix_ccs* matrix_2;

    jmtxds_matrix_ccs_new(&matrix_1, 5, 5, 0, NULL);

    {
        const uint32_t i[] = { 0 };
        const double v[] = { 1.0f };
        jmtxds_matrix_ccs_set_col(matrix_1, 0, 1, i, v);
    }

    {
        const uint32_t i[] = { 0, 1 };
        const double v[] = { 0.75f, 0.25f };
        jmtxds_matrix_ccs_set_col(matrix_1, 1, 2, i, v);
    }

    {
        const uint32_t i[] = { 1, 2, 3 };
        const double v[] = { 1, 2, 1};
        jmtxds_matrix_ccs_set_col(matrix_1, 2, 3, i, v);
    }

    {
        const uint32_t i[] = { 3, 4 };
        const double v[] = { 0.75f, 0.25f };
        jmtxds_matrix_ccs_set_col(matrix_1, 3, 2, i, v);
    }

    {
        const uint32_t i[] = { 4 };
        const double v[] = {1};
        jmtxds_matrix_ccs_set_col(matrix_1, 4, 1, i, v);
    }
    jmtxds_matrix_ccs_transpose(matrix_1, &matrix_2, NULL);

    printf("\n\nBREAK 2\n");
    print_ccsd_matrix(matrix_1);
    print_ccsd_matrix(matrix_2);

    uint32_t err_count = 0;

    jmtxd_matrix_ccs* tmp;
    jmtxds_matrix_ccs_transpose(matrix, &tmp, NULL);
    print_ccsd_matrix(matrix);
    print_ccsd_matrix(tmp);
    for (uint32_t i = 0; i < 25; ++i)
    {
        for (uint32_t j = 0; j < 25; ++j)
        {
            double x1,x2;
            x1 = jmtxds_matrix_ccs_get_entry(matrix, i, j, &x1);
            x2 = jmtxds_matrix_ccs_get_entry(tmp, j, i, &x2);
            if (x1 != x2)
            {
                err_count += 1;
            }
        }
    }
    printf("Error count for transpose: %u\n", err_count);
    jmtxds_matrix_ccs_destroy(tmp);

    jmtxds_matrix_ccs_copy(matrix, &tmp, NULL);
    err_count = 0;
    for (uint32_t i = 0; i < 25; ++i)
    {
        for (uint32_t j = 0; j < 25; ++j)
        {
            double x1,x2;
            jmtxds_matrix_ccs_get_entry(matrix, i, j, &x1);
            jmtxds_matrix_ccs_get_entry(tmp, i, j, &x2);
            if (x1 != x2)
            {
                err_count += 1;
            }
        }
    }
    printf("Error count for copy: %u\n", err_count);

    jmtxds_matrix_ccs_destroy(tmp);


    jmtxds_matrix_ccs_destroy(matrix_1);
    jmtxds_matrix_ccs_destroy(matrix_2);



    const double v1[25] = {1, 2, 3, 4, 5, [10] = 69.0f};
    double y[25];

    jmtxds_matrix_ccs_vector_multiply(matrix, v1, y);
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
//    printf("Matrix is currently using %zu bytes of memory\n", matrix_ccs_memory_usage(matrix));
    jmtxds_matrix_ccs_destroy(matrix);
    printf("Long double size: %zu bytes\n", sizeof(long double));
    return 0;
}

