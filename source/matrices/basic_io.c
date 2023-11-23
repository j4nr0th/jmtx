//
// Created by jan on 23.11.2023.
//

#include "basic_io.h"
#include "sparse_row_compressed_internal.h"
#include "sparse_column_compressed_internal.h"

#include <stdio.h>
#include <inttypes.h>

jmtx_result jmtx_matrix_crs_to_file(const jmtx_matrix_crs* mtx, const char* filename)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTX_TYPE_CRS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (!filename)
    {
        return JMTX_RESULT_NULL_PARAM;
    }

    FILE* const f_out = fopen(filename, "w");
    if (!f_out)
    {
        return JMTX_RESULT_NO_FILE;
    }

    fprintf(f_out, "Matrix of type %s %"PRIu32" x %"PRIu32" with %"PRIu32" entries",
            jmtx_matrix_type_to_str(mtx->base.type), mtx->base.rows, mtx->base.cols, mtx->n_entries);
    for (uint32_t i = 0; i < mtx->base.rows; ++i)
    {
        float* elements = NULL;
        uint32_t* indices = NULL;
        const uint32_t count = jmtx_matrix_crs_get_row(mtx, i, &indices, &elements);
        fprintf(f_out, "\nRow %"PRIu32":", i);
        for (uint32_t j = 0; j < count; ++j)
        {
            fprintf(f_out, "\t(%"PRIu32", %+.8e)", indices[j], elements[j]);
        }
    }
    fprintf(f_out, "\n");

    fclose(f_out);

    return JMTX_RESULT_SUCCESS;
}

jmtx_result jmtx_matrix_crs_to_file_explicit(const jmtx_matrix_crs* mtx, const char* filename)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTX_TYPE_CRS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (!filename)
    {
        return JMTX_RESULT_NULL_PARAM;
    }

    FILE* const f_out = fopen(filename, "w");
    if (!f_out)
    {
        return JMTX_RESULT_NO_FILE;
    }

    for (uint32_t i = 0; i < mtx->base.rows; ++i)
    {
        float* elements = NULL;
        uint32_t* indices = NULL;
        const uint32_t count = jmtx_matrix_crs_get_row(mtx, i, &indices, &elements);
        for (uint32_t j = 0, k = 0; j < mtx->base.cols; ++j)
        {
            if (k < count && indices[k] == j)
            {
                fprintf(f_out, "%+.8e\t", elements[j]);
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


jmtx_result jmtx_matrix_ccs_to_file(const jmtx_matrix_ccs* mtx, const char* filename)
{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTX_TYPE_CCS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (!filename)
    {
        return JMTX_RESULT_NULL_PARAM;
    }

    FILE* const f_out = fopen(filename, "w");
    if (!f_out)
    {
        return JMTX_RESULT_NO_FILE;
    }

    fprintf(f_out, "Matrix of type %s %"PRIu32" x %"PRIu32" with %"PRIu32" entries",
            jmtx_matrix_type_to_str(mtx->base.type), mtx->base.rows, mtx->base.cols, mtx->n_entries);
    for (uint32_t i = 0; i < mtx->base.cols; ++i)
    {
        float* elements = NULL;
        uint32_t* indices = NULL;
        const uint32_t count = jmtx_matrix_ccs_get_col(mtx, i, &indices, &elements);
        fprintf(f_out, "\nCol %"PRIu32":", i);
        for (uint32_t j = 0; j < count; ++j)
        {
            fprintf(f_out, "\t(%"PRIu32", %+.8e)", indices[j], elements[j]);
        }
    }
    fprintf(f_out, "\n");

    fclose(f_out);

    return JMTX_RESULT_SUCCESS;
}

jmtx_result jmtx_matrix_ccs_to_file_explicit(const jmtx_matrix_ccs* mtx, const char* filename)

{
    if (!mtx)
    {
        return JMTX_RESULT_NULL_PARAM;
    }
    if (mtx->base.type != JMTX_TYPE_CRS)
    {
        return JMTX_RESULT_WRONG_TYPE;
    }
    if (!filename)
    {
        return JMTX_RESULT_NULL_PARAM;
    }

    FILE* const f_out = fopen(filename, "w");
    if (!f_out)
    {
        return JMTX_RESULT_NO_FILE;
    }

    for (uint32_t i = 0; i < mtx->base.rows; ++i)
    {
        for (uint32_t j = 0, k = 0; j < mtx->base.cols; ++j)
        {
            fprintf(f_out, "%+.8e\t", jmtx_matrix_ccs_get_entry(mtx, i, j));
        }
        fprintf(f_out, "\n");
    }

    fclose(f_out);

    return JMTX_RESULT_SUCCESS;
}
