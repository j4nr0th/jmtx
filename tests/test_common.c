#include "test_common.h"
#include <stdio.h>
#include <inttypes.h>
#include <complex.h>

#ifdef JMTX_TYPE_DEFINED_FLOAT
void jmtxf_print_crs_matrix(const jmtxf_matrix_crs *mtx)
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
            const double x = jmtxf_matrix_crs_get_entry(mtx, i, j);
            printf("%10g ", x);
        }
        printf("] - %u", mtx->end_of_row_offsets[i] - (i ? mtx->end_of_row_offsets[i - 1] : 0));
        printf("\n");
    }

    printf("]\n");
}
void jmtxf_print_ccs_matrix(const jmtxf_matrix_ccs *mtx)
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
            const double x = jmtxf_matrix_ccs_get_entry(mtx, i, j);
            printf("%10g ", x);
        }
        const uint32_t n_row = jmtxf_matrix_ccs_elements_in_row(mtx, i);
        printf("] - %u", n_row);
        printf("\n");
    }
    printf("\t ");
    for (uint32_t i = 0; i < mtx->base.cols; ++i)
    {
        printf("%5" PRIu32 " ", (i == 0 ? mtx->end_of_column_offsets[0]
                                        : mtx->end_of_column_offsets[i] - mtx->end_of_column_offsets[i - 1]));
    }

    printf("\n]\n");
}
void jmtxf_print_brm_matrix(const jmtxf_matrix_brm *mtx)
{
    for (uint_fast32_t i = 0; i < mtx->base.rows; ++i)
    {
        float *p_vals;
        uint_fast32_t first = jmtxf_matrix_brm_first_pos_in_row(mtx, i);
        uint_fast32_t last = jmtxf_matrix_brm_get_row(mtx, i, &p_vals) + first;
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
void jmtxf_print_cds_matrix(const jmtxf_matrix_cds *mtx)
{
    printf("\nMatrix:\n[\n");
    for (uint_fast32_t i = 0; i < mtx->base.rows; ++i)
    {
        printf("\t[");
        for (uint32_t j = 0; j < mtx->base.cols; ++j)
        {
            const double x = jmtxf_matrix_cds_get_entry(mtx, i, j);
            printf("%10g ", x);
        }
        printf("]\n");
    }
    printf("]\n");
}
void jmtxf_print_drm_matrix(const jmtxf_matrix_drm *mtx)
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

#endif // JMTX_TYPE_DEFINED_FLOAT
#ifdef JMTX_TYPE_DEFINED_DOUBLE
void jmtxd_print_crs_matrix(const jmtxd_matrix_crs *mtx)
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
void jmtxd_print_ccs_matrix(const jmtxd_matrix_ccs *mtx)
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
        printf("%5" PRIu32 " ", (i == 0 ? mtx->end_of_column_offsets[0]
                                        : mtx->end_of_column_offsets[i] - mtx->end_of_column_offsets[i - 1]));
    }

    printf("\n]\n");
}
void jmtxd_print_brm_matrix(const jmtxd_matrix_brm *mtx)
{
    for (uint_fast32_t i = 0; i < mtx->base.rows; ++i)
    {
        double *p_vals;
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
void jmtxd_print_cds_matrix(const jmtxd_matrix_cds *mtx)
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
void jmtxd_print_drm_matrix(const jmtxd_matrix_drm *mtx)
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
#endif // JMTX_TYPE_DEFINED_DOUBLE
#ifdef JMTX_TYPE_DEFINED_CFLOAT
void jmtxc_print_crs_matrix(const jmtxc_matrix_crs *mtx)
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
            const _Complex float x = jmtxc_matrix_crs_get_entry(mtx, i, j);
            printf("%10g%+10g ", creal(x), cimag(x));
        }
        printf("] - %u", mtx->end_of_row_offsets[i] - (i ? mtx->end_of_row_offsets[i - 1] : 0));
        printf("\n");
    }

    printf("]\n");
}
void jmtxc_print_ccs_matrix(const jmtxc_matrix_ccs *mtx)
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
            const _Complex float x = jmtxc_matrix_ccs_get_entry(mtx, i, j);
            printf("%10g%+10g ", creal(x), cimag(x));
        }
        const uint32_t n_row = jmtxc_matrix_ccs_elements_in_row(mtx, i);
        printf("] - %u", n_row);
        printf("\n");
    }
    printf("\t ");
    for (uint32_t i = 0; i < mtx->base.cols; ++i)
    {
        printf("%5" PRIu32 " ", (i == 0 ? mtx->end_of_column_offsets[0]
                                        : mtx->end_of_column_offsets[i] - mtx->end_of_column_offsets[i - 1]));
    }

    printf("\n]\n");
}
void jmtxc_print_brm_matrix(const jmtxc_matrix_brm *mtx)
{
    for (uint_fast32_t i = 0; i < mtx->base.rows; ++i)
    {
        _Complex float *p_vals;
        uint_fast32_t first = jmtxc_matrix_brm_first_pos_in_row(mtx, i);
        uint_fast32_t last = jmtxc_matrix_brm_get_row(mtx, i, &p_vals) + first;
        uint_fast32_t j = 0;
        while (j < first)
        {
            printf(" %10g", 0.0f);
            j += 1;
        }
        uint_fast32_t k = 0;
        while (j < last)
        {
            printf(" %10g%+10g", creal(p_vals[k]), cimag(p_vals[k]));
            k += 1;
            j += 1;
        }
        while (j < mtx->base.cols)
        {
            printf(" %10g%+10g", 0.0, 0.0);
            j += 1;
        }
        printf("\t- %u\n", (unsigned)(last - first));
    }
}
void jmtxc_print_cds_matrix(const jmtxc_matrix_cds *mtx)
{
    printf("\nMatrix:\n[\n");
    for (uint_fast32_t i = 0; i < mtx->base.rows; ++i)
    {
        printf("\t[");
        for (uint32_t j = 0; j < mtx->base.cols; ++j)
        {
            const _Complex float x = jmtxc_matrix_cds_get_entry(mtx, i, j);
            printf("%10g%+10g ", creal(x), cimag(x));
        }
        printf("]\n");
    }
    printf("]\n");
}
void jmtxc_print_drm_matrix(const jmtxc_matrix_drm *mtx)
{
    printf("\nMatrix:\n[\n");
    for (uint_fast32_t i = 0; i < mtx->base.rows; ++i)
    {
        printf("\t[");
        for (uint32_t j = 0; j < mtx->base.cols; ++j)
        {
            const _Complex float x = mtx->values[j + mtx->base.cols * (mtx->permutations ? mtx->permutations[i] : i)];
            printf("%10g%+10g ", creal(x), cimag(x));
        }
        printf("]\n");
    }
    printf("]\n");
}
#endif // JMTX_TYPE_DEFINED_CFLOAT

#ifdef JMTX_TYPE_DEFINED_CDOUBLE
void jmtxz_print_crs_matrix(const jmtxz_matrix_crs *mtx)
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
void jmtxz_print_ccs_matrix(const jmtxz_matrix_ccs *mtx)
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
        printf("%5" PRIu32 " ", (i == 0 ? mtx->end_of_column_offsets[0]
                                        : mtx->end_of_column_offsets[i] - mtx->end_of_column_offsets[i - 1]));
    }

    printf("\n]\n");
}
void jmtxz_print_brm_matrix(const jmtxz_matrix_brm *mtx)
{
    for (uint_fast32_t i = 0; i < mtx->base.rows; ++i)
    {
        _Complex double *p_vals;
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
            printf(" %10g%+10g", creal(p_vals[k]), cimag(p_vals[k]));
            k += 1;
            j += 1;
        }
        while (j < mtx->base.cols)
        {
            printf(" %10g%+10g", 0.0, 0.0);
            j += 1;
        }
        printf("\t- %u\n", (unsigned)(last - first));
    }
}
void jmtxz_print_cds_matrix(const jmtxz_matrix_cds *mtx)
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
void jmtxz_print_drm_matrix(const jmtxz_matrix_drm *mtx)
{
    printf("\nMatrix:\n[\n");
    for (uint_fast32_t i = 0; i < mtx->base.rows; ++i)
    {
        printf("\t[");
        for (uint32_t j = 0; j < mtx->base.cols; ++j)
        {
            const _Complex double x = mtx->values[j + mtx->base.cols * (mtx->permutations ? mtx->permutations[i] : i)];
            printf("%10g%+10g ", creal(x), cimag(x));
        }
        printf("]\n");
    }
    printf("]\n");
}
#endif // JMTX_TYPE_DEFINED_CDOUBLE

void jmtxf_print_vec(unsigned n, const float x[JMTX_ARRAY_ATTRIB(static n)])
{
    printf("\nVector:\n[\n");
    for (uint_fast32_t i = 0; i < n; ++i)
    {
        printf("%10g\n", (double)x[i]);
    }
    printf("]\n");
}

void jmtxd_print_vec(unsigned n, const double x[JMTX_ARRAY_ATTRIB(static n)])
{
    printf("\nVector:\n[\n");
    for (uint_fast32_t i = 0; i < n; ++i)
    {
        printf("%10g\n", x[i]);
    }
    printf("]\n");
}

void jmtxc_print_vec(unsigned n, const _Complex float x[JMTX_ARRAY_ATTRIB(static n)])
{
    printf("\nVector:\n[\n");
    for (uint_fast32_t i = 0; i < n; ++i)
    {
        printf("%10g%+10g\n", creal(x[i]), cimag(x[i]));
    }
    printf("]\n");
}

void jmtxz_print_vec(unsigned n, const _Complex double x[JMTX_ARRAY_ATTRIB(static n)])
{
    printf("\nVector:\n[\n");
    for (uint_fast32_t i = 0; i < n; ++i)
    {
        printf("%10g%+10g\n", creal(x[i]), cimag(x[i]));
    }
    printf("]\n");
}

int are_close(const JMTX_SCALAR_T v1, const JMTX_SCALAR_T v2, const JMTX_REAL_T relative_tol, const JMTX_REAL_T abs_tol)
{
    if (v1 == v2)
    {
        return 1;
    }
    ASSERT(abs_tol >= 0.0f);
    ASSERT(relative_tol >= 0.0f);
    const double dif = JMTX_ABS(v1 - v2);
    if (dif < abs_tol)
    {
        return 1;
    }

    const JMTX_REAL_T mag1 = JMTX_ABS(v1);
    const JMTX_REAL_T mag2 = JMTX_ABS(v2);

    if (mag1 > mag2)
    {
        return dif < relative_tol * mag2;
    }

    return dif < relative_tol * mag1;
}
