#include "basic_io.h"
#include "sparse_column_compressed.h"
#include "sparse_row_compressed.h"

#include <inttypes.h>
#include <stdio.h>

jmtx_result JMTX_NAME_TYPED(matrix_crs_to_file)(const JMTX_NAME_TYPED(matrix_crs) * mtx, const char *filename)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTXD_TYPE_CRS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (!filename)
    {
        return JMTX_RESULT_NULL_PARAM;
    }

    FILE *const f_out = fopen(filename, "w");
    if (!f_out)
    {
        return JMTX_RESULT_NO_FILE;
    }

    fprintf(f_out, "Matrix of type %s %" PRIu32 " x %" PRIu32 " with %zu entries",
            jmtx_matrix_type_to_str(mtx->base.type), mtx->base.rows, mtx->base.cols, mtx->n_entries);
    for (JMTX_INDEX_T i = 0; i < mtx->base.rows; ++i)
    {
        JMTX_SCALAR_T *elements = NULL;
        JMTX_INDEX_T *indices = NULL;
        const JMTX_INDEX_T count = JMTX_NAME_TYPED(matrix_crs_get_row)(mtx, i, &indices, &elements);
        fprintf(f_out, "\nRow %" PRIu32 ":", i);
        for (JMTX_INDEX_T j = 0; j < count; ++j)
        {
            fprintf(f_out, "\t(%" PRIu32 ", %+.8e)", indices[j], (double)elements[j]);
        }
    }
    fprintf(f_out, "\n");

    fclose(f_out);

    return JMTX_RESULT_SUCCESS;
}

jmtx_result JMTX_NAME_TYPED(matrix_crs_to_file_explicit)(const JMTX_NAME_TYPED(matrix_crs) * mtx, const char *filename)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTXD_TYPE_CRS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (!filename)
    {
        return JMTX_RESULT_NULL_PARAM;
    }

    FILE *const f_out = fopen(filename, "w");
    if (!f_out)
    {
        return JMTX_RESULT_NO_FILE;
    }

    for (JMTX_INDEX_T i = 0; i < mtx->base.rows; ++i)
    {
        JMTX_SCALAR_T *elements = NULL;
        JMTX_INDEX_T *indices = NULL;
        const JMTX_INDEX_T count = JMTX_NAME_TYPED(matrix_crs_get_row)(mtx, i, &indices, &elements);
        for (JMTX_INDEX_T j = 0, k = 0; j < mtx->base.cols; ++j)
        {
            if (k < count && indices[k] == j)
            {
                fprintf(f_out, "%+.8e\t", (double)elements[j]);
            }
            else
            {
                fprintf(f_out, "%+.8e\t", 0.0f);
            }
        }
        fprintf(f_out, "\n");
    }

    fclose(f_out);

    return JMTX_RESULT_SUCCESS;
}

jmtx_result JMTX_NAME_TYPED(matrix_ccs_to_file)(const JMTX_NAME_TYPED(matrix_ccs) * mtx, const char *filename)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTXD_TYPE_CCS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (!filename)
    {
        return JMTX_RESULT_NULL_PARAM;
    }

    FILE *const f_out = fopen(filename, "w");
    if (!f_out)
    {
        return JMTX_RESULT_NO_FILE;
    }

    fprintf(f_out, "Matrix of type %s %" PRIu32 " x %" PRIu32 " with %zu entries",
            jmtx_matrix_type_to_str(mtx->base.type), mtx->base.rows, mtx->base.cols, mtx->n_entries);
    for (JMTX_INDEX_T i = 0; i < mtx->base.cols; ++i)
    {
        JMTX_SCALAR_T *elements = NULL;
        JMTX_INDEX_T *indices = NULL;
        const JMTX_INDEX_T count = JMTX_NAME_TYPED(matrix_ccs_get_col)(mtx, i, &indices, &elements);
        fprintf(f_out, "\nCol %" PRIu32 ":", i);
        for (JMTX_INDEX_T j = 0; j < count; ++j)
        {
            fprintf(f_out, "\t(%" PRIu32 ", %+.8e)", indices[j], (double)elements[j]);
        }
    }
    fprintf(f_out, "\n");

    fclose(f_out);

    return JMTX_RESULT_SUCCESS;
}

jmtx_result JMTX_NAME_TYPED(matrix_ccs_to_file_explicit)(const JMTX_NAME_TYPED(matrix_ccs) * mtx, const char *filename)

{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTXD_TYPE_CRS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (!filename)
    {
        return JMTX_RESULT_NULL_PARAM;
    }

    FILE *const f_out = fopen(filename, "w");
    if (!f_out)
    {
        return JMTX_RESULT_NO_FILE;
    }

    for (JMTX_INDEX_T i = 0; i < mtx->base.rows; ++i)
    {
        for (JMTX_INDEX_T j = 0; j < mtx->base.cols; ++j)
        {
            fprintf(f_out, "%+.8e\t", (double)JMTX_NAME_TYPED(matrix_ccs_get_entry)(mtx, i, j));
        }
        fprintf(f_out, "\n");
    }

    fclose(f_out);

    return JMTX_RESULT_SUCCESS;
}
