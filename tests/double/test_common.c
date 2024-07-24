// Automatically generated from tests/float/test_common.c on Fri Dec  1 06:43:09 2023
//
// Created by jan on 15.6.2022.
//

#include "test_common.h"
#include "../../include/jmtx/double/matrices/sparse_row_compressed.h"
#include "../../include/jmtx/double/matrices/sparse_column_compressed.h"
#include <stdio.h>
#include "../../source/double/matrices/sparse_row_compressed_internal.h"
#include "../../source/double/matrices/sparse_column_compressed_internal.h"
#include "../../source/double/matrices/band_row_major_internal.h"
#include "../../source/double/matrices/sparse_diagonal_compressed_internal.h"
#include "../../include/jmtx/double/matrices/sparse_diagonal_compressed_safe.h"
#include "../../source/double/matrices/dense_row_major_internal.h"
#include <inttypes.h>


void print_crsd_matrix(const jmtxd_matrix_crs* mtx)
{
    printf("Entries:");
    for (uint32_t i = 0, l = 0; i < mtx->n_entries; ++i)
    {
        while (l < mtx->base.rows && mtx->end_of_row_offsets[l] <= i)
        {
            l += 1;
        }
        printf(" (%u, %u, %g)", l, mtx->indices[i], mtx->values[i]);
    }
    printf("\nEnd of row offsets:");
    for (uint32_t i = 0; i < mtx->base.rows; ++i)
    {
        printf(" %u,", mtx->end_of_row_offsets[i]);
    }
    printf("\nMatrix:\n[\n");

    for (uint32_t i = 0; i < mtx->base.rows; ++i)
    {
        printf("\t[");
        for (uint32_t j = 0; j < mtx->base.cols; ++j)
        {
            const double x = jmtxd_matrix_crs_get_entry(mtx, i, j);
            printf("%10g ", x);
        }
        printf("] - %u", mtx->end_of_row_offsets[i] - (i ? mtx->end_of_row_offsets[i - 1] : 0));
        printf("\n");
    }

    printf("]\n");
}

void print_ccsd_matrix(const jmtxd_matrix_ccs* mtx)
{
    printf("values:");
    for (uint32_t i = 0, l = 0; i < mtx->n_entries; ++i)
    {
        if (l < mtx->base.rows && mtx->end_of_column_offsets[l] <= i)
        {
            l += 1;
        }
        printf(" (%u, %u, %g)", l, mtx->indices[i], mtx->values[i]);
    }
    printf("\nelement offsets:");
    for (uint32_t i = 0; i < mtx->base.cols; ++i)
    {
        printf(" %u,", mtx->end_of_column_offsets[i]);
    }
    printf("\nMatrix:\n[\n");

    for (uint32_t i = 0; i < mtx->base.rows; ++i)
    {
        printf("\t[");
        for (uint32_t j = 0; j < mtx->base.cols; ++j)
        {
            const double x = jmtxd_matrix_ccs_get_entry(mtx, i, j);
            printf("%10g ", x);
        }
        const uint32_t n_row = jmtxd_matrix_ccs_elements_in_row(mtx, i);
        printf("] - %u", n_row);
        printf("\n");
    }
    printf("\t ");
    for (uint32_t i = 0; i < mtx->base.cols; ++i)
    {
        printf("%5"PRIu32" ", (i == 0 ? mtx->end_of_column_offsets[0] : mtx->end_of_column_offsets[i] - mtx->end_of_column_offsets[i - 1]));
    }

    printf("\n]\n");
}

void print_brmd_matrix(const jmtxd_matrix_brm* mtx)
{
    for (uint_fast32_t i = 0; i < mtx->base.rows; ++i)
    {
        double* p_vals;
        uint_fast32_t first = jmtxd_matrix_brm_first_pos_in_row(mtx, i);
        uint_fast32_t last = jmtxd_matrix_brm_get_row(mtx, i, &p_vals) + first;
        uint_fast32_t j = 0;
        while (j < first)
        {
            printf(" %10g", 0.0f);
            j += 1;
        }
        uint_fast32_t k = 0;
        while (j < last)
        {
            printf(" %10g", p_vals[k++]);
            j += 1;
        }
        while (j < mtx->base.cols)
        {
            printf(" %10g", 0.0f);
            j += 1;
        }
        printf("\t- %u\n", (unsigned)(last - first));
    }
}

void print_cdsd_matrix(const jmtxd_matrix_cds* mtx)
{
    printf("\nMatrix:\n[\n");
    for (uint_fast32_t i = 0; i < mtx->base.rows; ++i)
    {
        printf("\t[");
        for (uint32_t j = 0; j < mtx->base.cols; ++j)
        {
            const double x = jmtxd_matrix_cds_get_entry(mtx, i, j);
            printf("%10g ", x);
        }
        printf("]\n");
    }
    printf("]\n");
}

void print_drmd_matrix(const struct jmtxd_matrix_drm_struct* mtx)
{
    printf("\nMatrix:\n[\n");
    for (uint_fast32_t i = 0; i < mtx->base.rows; ++i)
    {
        printf("\t[");
        for (uint32_t j = 0; j < mtx->base.cols; ++j)
        {
            const float x = mtx->values[j + mtx->base.cols * (mtx->permutations ? mtx->permutations[i] : i)];
            printf("%10g ", x);
        }
        printf("]\n");
    }
    printf("]\n");
}

void print_vecd(unsigned n, const double x[JMTX_ARRAY_ATTRIB(static n)])
{
    printf("\nVector:\n[\n");
    for (uint_fast32_t i = 0; i < n; ++i)
    {
        printf("%10g\n", x[i]);
    }
    printf("]\n");
}

