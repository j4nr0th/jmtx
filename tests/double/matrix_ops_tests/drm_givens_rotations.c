//
// Created by jan on 18.7.2024.
//
#include <math.h>
#include "../../../include/jmtx/double/matrices/dense_row_major.h"
#include "../../cfloat/test_common.h"
#include "../test_common.h"
#include <stdlib.h>

enum
{
    N = 10
};

static void eliminate_entry_upper(jmtxd_matrix_drm* mat, unsigned i, unsigned j)
{
    ASSERT(i < N - 1);
    ASSERT(j < N && j > 0);
    const double a11 = jmtxd_matrix_drm_get_entry(mat, i, i);
    const double a12 = jmtxd_matrix_drm_get_entry(mat, i, j);

    const double mag = hypot(a11, a12);
    const double s = -a12 / mag;
    const double c = +a11 / mag;

    jmtxd_matrix_drm_givens_rotation_right(mat, i, j, c, s);
}

static void eliminate_entry_lower(jmtxd_matrix_drm* mat, unsigned i, unsigned j)
{
    ASSERT(j < i);
    ASSERT(i < N && i > 0);
    const double a11 = jmtxd_matrix_drm_get_entry(mat, j, j);
    const double a21 = jmtxd_matrix_drm_get_entry(mat, i, j);

    const double mag = hypot(a11, a21);
    const double s = -a21 / mag;
    const double c = +a11 / mag;

    jmtxd_matrix_drm_givens_rotation_left(mat, j, i, c, s);
}

static void eliminate_entry_upper2(jmtxd_matrix_drm* mat, unsigned i, unsigned j, jmtxd_matrix_drm* q)
{
    ASSERT(i < N - 1);
    ASSERT(j < N && j > 0);
    const double a11 = jmtxd_matrix_drm_get_entry(mat, i, i);
    const double a12 = jmtxd_matrix_drm_get_entry(mat, i, j);

    const double mag = hypot(a11, a12);
    const double s = -a12 / mag;
    const double c = +a11 / mag;

    jmtxd_matrix_drm_givens_rotation_right(mat, i, j, c, s);
    jmtxd_matrix_drm_givens_rotation_left(q, i, j, c, s);
}

static void eliminate_entry_upper3(jmtxd_matrix_drm* mat, unsigned i, unsigned j, jmtxd_matrix_drm* q, const double k)
{
    ASSERT(i < N - 1);
    ASSERT(j < N && j > 0);
    const double a11 = jmtxd_matrix_drm_get_entry(mat, i, i) - k;
    const double a12 = jmtxd_matrix_drm_get_entry(mat, i, j);

    const double mag = hypot(a11, a12);
    const double s = -a12 / mag;
    const double c = +a11 / mag;

    jmtxd_matrix_drm_givens_rotation_right(mat, i, j, c, s);
    jmtxd_matrix_drm_givens_rotation_left(q, i, j, c, s);
}

static inline void eigenvals2x2(const double a11, const double a12, const double a21, const double a22, double out[2])
{
    const double f = (a11 + a22) / 2;
    const double s = sqrt(f * f + a12*a21 - a11*a22);
    out[0] = f + s;
    out[1] = f - s;
}

static void qr_algo(jmtxd_matrix_drm* mat, unsigned max_itr, double out[])
{
    jmtxd_matrix_drm* q, *l, *t;
    jmtx_result mtx_res = jmtxd_matrix_drm_copy(mat, &l, NULL);
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    mtx_res = jmtxd_matrix_drm_new(&q, N, N, NULL, NULL);
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    mtx_res = jmtxd_matrix_drm_new(&t, N, N, NULL, NULL);
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);

    unsigned p = N;
    double ev = jmtxd_matrix_drm_get_entry(l, p-1, p-1);
    unsigned itr;
    double ev_est[2];
    for (itr = 0; itr < max_itr && p > 2; ++itr)
    {
        jmtxd_matrix_drm_set_all_entries(q, 0);
        for (unsigned i = 0; i < N; ++i)
            jmtxd_matrix_drm_set_entry(q, i, i, 1);
        // print_drmd_matrix(l);
        // print_drmd_matrix(q);
        eigenvals2x2(
            jmtxd_matrix_drm_get_entry(l, p-2, p-2),
            jmtxd_matrix_drm_get_entry(l, p-2, p-1),
            jmtxd_matrix_drm_get_entry(l, p-1, p-2),
            ev,
            ev_est);

        // jmtxd_matrix_drm_shift_diagonal(l, -ev_est[1]);
        for (unsigned i = 0; i < p - 1; ++i)
        {
            for (unsigned j = i + 1; j < p; ++j)
            {
                // eliminate_entry_upper2(l, i, j, q);
                eliminate_entry_upper3(l, i, j, q, ev_est[1]);
                // print_drmd_matrix(l);

            }
        }
        // print_drmd_matrix(l);
        // jmtxd_matrix_drm_transpose_inplace(q, out);
        jmtxd_matrix_drm_multiply_matrix(q, l  , t);
        {
            jmtxd_matrix_drm* tmp = t;
            t = l;
            l = tmp;
        }
        // jmtxd_matrix_drm_multiply_matrix2(l, q, t);
        // jmtxd_matrix_drm_shift_diagonal(l, +ev_est[1]);
        const double new_ev = jmtxd_matrix_drm_get_entry(l, p-1, p-1);

        const double de = new_ev - ev;
        if (fabs(de) < 1e-15)
        {
            printf("Found eigenvalue %g on iteration %u with difference of %e\n", new_ev, itr, de);
            p -= 1;
        }
        else
        {
            ev = new_ev;
        }
        // print_drmd_matrix(t);
    }
    // print_drmd_matrix(l);
    const double a = jmtxd_matrix_drm_get_entry(l, 0, 0);
    const double b = jmtxd_matrix_drm_get_entry(l, 0, 1);
    const double c = jmtxd_matrix_drm_get_entry(l, 1, 0);
    const double d = jmtxd_matrix_drm_get_entry(l, 1, 1);
    eigenvals2x2(a, b, c, d, out);

    for (unsigned i = 2; i < N; ++i)
    {
        out[i] = jmtxd_matrix_drm_get_entry(l, i, i);
    }
    printf("Finished on iteration %u\n", itr);
}

int main()
{
    jmtxd_matrix_drm* mtx;
    jmtx_result mtx_res;

    MATRIX_TEST_CALL(jmtxd_matrix_drm_new(&mtx, N, N, (double []){0.0}, NULL));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);

    for (unsigned i = 0; i < N; ++i)
        jmtxd_matrix_drm_set_entry(mtx, i, i, 1.0);

    print_drmd_matrix(mtx);

    MATRIX_TEST_CALL(jmtxds_matrix_drm_givens_rotation_left(mtx, 1, 2, cos(0.4), sin(0.4)));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    print_drmd_matrix(mtx);

    MATRIX_TEST_CALL(jmtxds_matrix_drm_givens_rotation_right(mtx, 1, 2, cos(0.4), sin(0.4)));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);

    print_drmd_matrix(mtx);
    srand(0);
    for (unsigned i = 0; i < N; ++i)
    {
        for (unsigned j = 0; j < N; ++j)
        {
            jmtxd_matrix_drm_set_entry(mtx, i, j, ((double)rand()/(double)RAND_MAX * 2 - 1));
        }
    }
    print_drmd_matrix(mtx);

    jmtxd_matrix_drm* trp;
    jmtxd_matrix_drm_transpose(mtx, &trp, NULL);
    jmtxd_matrix_drm* spdm, *eye, *spdm2;
    MATRIX_TEST_CALL(jmtxd_matrix_drm_new(&spdm, N, N, NULL, NULL));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    MATRIX_TEST_CALL(jmtxd_matrix_drm_multiply_matrix(trp, mtx, spdm));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    jmtxd_matrix_drm_destroy(trp);
    print_drmd_matrix(spdm);
    MATRIX_TEST_CALL(jmtxd_matrix_drm_copy(spdm, &spdm2, NULL));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);

    MATRIX_TEST_CALL(jmtxd_matrix_drm_new(&eye, N, N, (double []){0.0}, NULL));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);

    for (unsigned i = 0; i < N; ++i)
        jmtxd_matrix_drm_set_entry(eye, i, i, 1.0);

    double results[N];
    const double real[N] =
        {
        12.947353812927712, 8.397417766375117, 4.522814669443053, 2.5762574981520463, 1.9548461558510135, 1.4955423031466695, 0.7171290281538076, 0.4648214913182648, 0.026795848501450264, 0.011904133726680758
        };

    // for (unsigned n = 90; n < 91; ++n)
    {
        qr_algo(spdm, 90, results);
        printf("eigvals for n = %u:\n\t", 90);
        for (unsigned i = 0; i < N; ++i)
        {
            printf("%g ", results[i]);
        }
        printf("\n");
    }

    printf("Errors:\n\t");
    for (unsigned i = 0; i < N; ++i)
    {
        double e = INFINITY;
        unsigned k = ~0u;
        for (unsigned j = 0; j < N; ++j)
        {
            const double new_e = fabs(results[i] - real[j]);
            if (new_e < e)
            {
                e = new_e;
                k = j;
            }
        }
        printf("(%u) %e| ", k, e);
    }
    printf("\n");

    const double test_mtx[10][10] = {
        {
            -1.43774796e+01, -1.62717870e+01, -9.82061561e+00, -1.69440328e+01,
            -1.36585536e+01, -1.85730258e+01, -2.14778763e+01, -1.61921978e+01,
            -2.26153246e+01, -1.64976273e+01,
        },
 { 8.87853218e+00,  1.02529999e+01,  2.37271179e+00,  5.17398860e+00,
   1.00826664e+01 , 1.11075517e+01 , 7.13036746e+00 , 1.17493454e+01,
   1.77385227e+01,  1.04822984e+01},
 { 1.78103937e+01 , 1.03459385e+01 , 1.36421969e+01 , 1.12666930e+01,
   1.34104608e+01 , 1.90105495e+01 , 1.62282831e+01 , 1.66403782e+01,
   2.52998740e+01 , 1.54449266e+01},
 {-2.34713027e+00 ,-1.75548765e-01 ,-1.84975549e+00 , 4.89821115e+00,
  -1.94769257e+00, -2.70141981e+00 ,-1.31135359e+00, -1.09114439e+00,
  -1.62483941e+00 ,-3.44355829e+00},
 {-1.28426281e+01, -7.11970557e+00 ,-6.80271675e+00, -9.94719089e+00,
  -5.86475873e+00 ,-1.53803551e+01, -1.19640075e+01, -1.21198316e+01,
  -2.04377617e+01 ,-1.10298175e+01},
 {-2.72166450e+00, -2.71507999e+00,  6.90919570e-02 ,-1.61678507e+00,
  -4.34207038e+00 , 3.10907001e+00 ,-2.72571998e+00 ,-6.61231055e+00,
  -9.13022538e+00 ,-1.14666503e+00},
 {-4.85345791e+00, -1.47750783e+00, -4.87786804e-01 , 1.33286497e-02,
  -5.81594367e+00, -6.36392117e+00,  4.00901832e+00, -8.32049466e+00,
  -1.18981128e+01, -5.34052768e+00},
 { 1.08147050e+00 , 2.96506565e+00 , 2.37352500e+00,  2.94022953e+00,
   6.97990259e-01, -6.06980191e-01,  2.29288467e+00,  4.90968338e+00,
  -3.68061375e+00,  1.75206477e+00,},
 { 1.05196554e+01 , 6.57305058e+00 , 4.59358487e+00 , 7.79203169e+00,
   9.38301381e+00 , 1.14361057e+01 , 8.98885188e+00,  1.20202018e+01,
   2.38855341e+01 , 8.09582798e+00},
 { 8.44920356e+00  ,5.62215091e+00 , 3.81473180e+00 , 5.18012195e+00,
   7.22907146e+00 , 8.41499314e+00 , 7.65680220e+00 , 9.62350611e+00,
   1.43238900e+01 , 1.05355246e+01}};
    for (unsigned i = 0; i < 10; ++i)
    {
        for (unsigned j = 0; j < 10; ++j)
        {
            jmtxd_matrix_drm_set_entry(mtx, i, j, test_mtx[i][j]);
        }
    }

    {
        qr_algo(mtx, 90, results);
        printf("eigvals for n = %u:\n\t", 90);
        for (unsigned i = 0; i < N; ++i)
        {
            printf("%g ", results[i]);
        }
        printf("\n");
    }

    printf("Errors:\n\t");
    for (unsigned i = 0; i < N; ++i)
    {
        double e = INFINITY;
        unsigned k = ~0u;
        for (unsigned j = 0; j < N; ++j)
        {
            const double new_e = fabs(results[i] - (double)(j + 1));
            if (new_e < e)
            {
                e = new_e;
                k = j;
            }
        }
        printf("(%u) %e| ", k, e);
    }
    printf("\n");

    exit(0);

    for (unsigned i = 0; i < N - 1; ++i)
    {
        for (unsigned j = i + 1; j < N; ++j)
        {
            eliminate_entry_upper2(spdm, i, j, eye);
        }
    }

    print_drmd_matrix(spdm);
    print_drmd_matrix(eye);

    jmtxd_matrix_drm* r;
    MATRIX_TEST_CALL(jmtxd_matrix_drm_new(&r, N, N, NULL, NULL));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    MATRIX_TEST_CALL(jmtxd_matrix_drm_multiply_matrix(eye, spdm, r));
    ASSERT(mtx_res == JMTX_RESULT_SUCCESS);
    print_drmd_matrix(r);


    for (unsigned j = 0; j < N - 1; ++j)
    {
        for (unsigned i = j + 1; i < N; ++i)
        {
            eliminate_entry_lower(spdm2, i, j);
        }
    }
    print_drmd_matrix(spdm2);




    jmtxd_matrix_drm_destroy(r);
    jmtxd_matrix_drm_destroy(spdm2);
    jmtxd_matrix_drm_destroy(spdm);
    jmtxd_matrix_drm_destroy(eye);
    jmtxd_matrix_drm_destroy(mtx);

    return 0;
}