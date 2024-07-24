// Automatically generated from tests/cfloat/test_common.c on Fri Dec  1 18:48:10 2023
// Automatically generated from tests/cdouble/test_common.c on Fri Dec  1 17:35:45 2023
//
// Created by jan on 15.6.2022.
//

#include "test_common.h"
#include "../../include/jmtx/cdouble/matrices/sparse_row_compressed.h"
#include "../../include/jmtx/cdouble/matrices/sparse_column_compressed.h"
#include <stdio.h>
#include "../../source/cdouble/matrices/sparse_row_compressed_internal.h"
#include "../../source/cdouble/matrices/sparse_column_compressed_internal.h"
#include "../../source/cdouble/matrices/band_row_major_internal.h"
#include "../../source/cdouble/matrices/sparse_diagonal_compressed_internal.h"
#include "../../source/cdouble/matrices/dense_row_major_internal.h"
#include <inttypes.h>
#include <complex.h>

void print_crsz_matrix(const jmtxz_matrix_crs* mtx)
{
    printf("Entries:");
    for (uint32_t i = 0, l = 0; i < mtx->n_entries; ++i)
    {
        while (l < mtx->base.rows && mtx->end_of_row_offsets[l] <= i)
        {
            l += 1;
        }
        printf(" (%u, %u, %g%+g)", l, mtx->indices[i], creal(mtx->values[i]), cimag(mtx->values[i]));
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
            const _Complex double x = jmtxz_matrix_crs_get_entry(mtx, i, j);
            printf("%10g%+10g ", creal(x), cimag(x));
        }
        printf("] - %u", mtx->end_of_row_offsets[i] - (i ? mtx->end_of_row_offsets[i - 1] : 0));
        printf("\n");
    }

    printf("]\n");
}

void print_ccsz_matrix(const jmtxz_matrix_ccs* mtx)
{
    printf("values:");
    for (uint32_t i = 0, l = 0; i < mtx->n_entries; ++i)
    {
        if (l < mtx->base.rows && mtx->end_of_column_offsets[l] <= i)
        {
            l += 1;
        }
        printf(" (%u, %u, %g%+g)", l, mtx->indices[i], creal(mtx->values[i]), cimag(mtx->values[i]));
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
            const _Complex double x = jmtxz_matrix_ccs_get_entry(mtx, i, j);
            printf("%10g%+10g ", creal(x), cimag(x));
        }
        const uint32_t n_row = jmtxz_matrix_ccs_elements_in_row(mtx, i);
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

void print_brmz_matrix(const jmtxz_matrix_brm* mtx)
{
    for (uint_fast32_t i = 0; i < mtx->base.rows; ++i)
    {
        _Complex double* p_vals;
        uint_fast32_t first = jmtxz_matrix_brm_first_pos_in_row(mtx, i);
        uint_fast32_t last = jmtxz_matrix_brm_get_row(mtx, i, &p_vals) + first;
        uint_fast32_t j = 0;
        while (j < first)
        {
            printf(" %10g", 0.0f);
            j += 1;
        }
        uint_fast32_t k = 0;
        while (j < last)
        {
            complex double x = p_vals[k++];
            printf("%10g%+10g ", creal(x), cimag(x));
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

void print_cdsz_matrix(const jmtxz_matrix_cds* mtx)
{
    printf("\nMatrix:\n[\n");
    for (uint_fast32_t i = 0; i < mtx->base.rows; ++i)
    {
        printf("\t[");
        for (uint32_t j = 0; j < mtx->base.cols; ++j)
        {
            const _Complex double x = jmtxz_matrix_cds_get_entry(mtx, i, j);
            printf("%10g%+10g ", creal(x), cimag(x));
        }
        printf("]\n");
    }
    printf("]\n");
}

void print_drmz_matrix(const struct jmtxz_matrix_drm_struct* mtx)
{
    printf("\nMatrix:\n[\n");
    for (uint_fast32_t i = 0; i < mtx->base.rows; ++i)
    {
        printf("\t[");
        for (uint32_t j = 0; j < mtx->base.cols; ++j)
        {
            const _Complex double x = mtx->values[j + mtx->base.cols * (mtx->permutations ? mtx->permutations[i] : i)];
            printf("%10g%+9gj\n", crealf(x), cimagf(x));
        }
        printf("]\n");
    }
    printf("]\n");
}

void print_vecz(unsigned n, const _Complex double x[JMTX_ARRAY_ATTRIB(static n)])
{
    printf("\nVector:\n[\n");
    for (uint_fast32_t i = 0; i < n; ++i)
    {
        printf("%10g%+9gj\n", crealf(x[i]), cimagf(x[i]));
    }
    printf("]\n");
}
