#include "mat_lzz_pX_multiply.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* returns trunc( trunc(a, dA+1)*c div x^dA, dB+1 )           */
/* assumes FFT prime and p large enough                       */
/* output may alias input; b does not have to be zero matrix  */
/* does not use Mat<zz_p> matrix multiplication               */
/*------------------------------------------------------------*/
void middle_product_evaluate_dense(Mat<zz_pX> & b, const Mat<zz_pX> & a, const Mat<zz_pX> & c, long dA, long dB)
{
    const long m = a.NumRows();
    const long n = a.NumCols();
    const long p = c.NumCols();

    long ell;

    Mat<zz_p> vA, vB, iv;
    vandermonde(vA, vB, iv, dA, dB);
    transpose(iv, iv);
    transpose(vB, vB);
    const long nb_points = vA.NumRows();

    // evaluations of matrix a
    Mat<zz_p> tmp_mat(INIT_SIZE, dA + 1, m * n);
    ell = 0;
    for (long i = 0; i < m; ++i)
        for (long j = 0; j < n; ++j, ++ell)
            for (long k = 0; k <= deg(a[i][j]); ++k)
                tmp_mat[dA-k][ell] = a[i][j][k];

    Mat<zz_p> valA;
    mul(valA, vA, tmp_mat);

    // evaluations of matrix c
    tmp_mat.SetDims(nb_points, n * p);
    ell = 0;
    for (long i = 0; i < n; ++i)
        for (long j = 0; j < p; ++j, ++ell)
        {
            const long d = deg(c[i][j]);
            for (long k = 0; k <= d; ++k)
                tmp_mat[k][ell] = c[i][j][k];
            for (long k = d+1; k < nb_points; ++k)
                clear(tmp_mat[k][ell]);
        }

    Mat<zz_p> valC;
    mul(valC, iv, tmp_mat);

    Mat<zz_p> valB(INIT_SIZE, nb_points, m * p);
    Mat<zz_p> valBp;

    Mat<zz_p> valAp(INIT_SIZE, m, n);
    Mat<zz_p> valCp(INIT_SIZE, n, p);
    for (long i = 0; i < nb_points; ++i)
    {
        ell = 0;
        for (long u = 0; u < m; ++u)
            for (long v = 0; v < n; ++v, ++ell)
                valAp[u][v] = valA[i][ell];

        ell = 0;
        for (long u = 0; u < n; ++u)
            for (long v = 0; v < p; ++v, ++ell)
                valCp[u][v] = valC[i][ell];

        mul(valBp, valAp, valCp);

        ell = 0;
        for (long u = 0; u < m; ++u)
            for (long v = 0; v < p; ++v, ++ell)
                valB[i][ell] = valBp[u][v];
    }
    mul(tmp_mat, vB, valB);

    b.SetDims(m, p);
    ell = 0;
    for (long u = 0; u < m; ++u)
        for (long v = 0; v < p; ++v, ++ell)
        {
            b[u][v].SetLength(dB + 1);
            for (long i = 0; i < dB + 1; ++i)
                b[u][v][i] = tmp_mat[i][ell];
            b[u][v].normalize();
        }
}

// experimental code, not correct for the moment
#if 0
void middle_product_evaluate_dense2(Mat<zz_pX> & b, const Mat<zz_pX> & a, const Mat<zz_pX> & c, long dA, long dB)
{
    const long m = a.NumRows();
    const long n = a.NumCols();
    const long p = c.NumCols();
    const long sA = dA/2+1;

    long ell;

    Mat<zz_p> vA, vB, iv;
    vandermonde2(vA, vB, iv, dA, dB);
    transpose(iv, iv);
    transpose(vB, vB);
    const long nb_points = vA.NumRows();

    // evaluations of matrix a, even/odd parts
    Mat<zz_p> tmp_mat_even(INIT_SIZE, sA, m * n);
    Mat<zz_p> tmp_mat_odd(INIT_SIZE, sA, m * n);
    ell = 0;
    for (long i = 0; i < m; ++i)
        for (long j = 0; j < n; ++j, ++ell)
        {
            const long d = deg(a[i][j]);
            for (long k = 0; 2*k < d; ++k)
            {
                tmp_mat_even[sA-1-k][ell] = a[i][j][2*k];
                tmp_mat_odd[sA-1-k][ell] = a[i][j][2*k+1];
            }
            if (d % 2 == 0)
                tmp_mat_even[sA-1-d/2][ell] = a[i][j][d];
        }

    Mat<zz_p> valAeven, valAodd;
    mul(valAeven, vA, tmp_mat_even);
    mul(valAodd, vA, tmp_mat_odd);

    // evaluations of matrix c
    tmp_mat_even.SetDims(nb_points, n * p);
    tmp_mat_odd.SetDims(nb_points, n * p);
    ell = 0;
    for (long i = 0; i < n; ++i)
        for (long j = 0; j < p; ++j, ++ell)
        {
            const long d = deg(c[i][j]);
            long k=0;
            for (; 2*k < d; ++k)
            {
                tmp_mat_even[k][ell] = c[i][j][2*k];
                tmp_mat_odd[k][ell] = c[i][j][2*k+1];
            }
            if (2*k == d)
            {
                tmp_mat_even[k][ell] = c[i][j][d];
                clear(tmp_mat_odd[k][ell]);
                ++k;
            }
            for (; k < nb_points; ++k)
            {
                clear(tmp_mat_even[k][ell]);
                clear(tmp_mat_odd[k][ell]);
            }
        }

    Mat<zz_p> valCeven, valCodd;
    mul(valCeven, iv, tmp_mat_even);
    mul(valCodd, iv, tmp_mat_odd);

    Mat<zz_p> valBeven(INIT_SIZE, nb_points, m * p);
    Mat<zz_p> valBodd(INIT_SIZE, nb_points, m * p);
    Mat<zz_p> valBpeven, valBpodd;

    Mat<zz_p> valApeven(INIT_SIZE, m, n);
    Mat<zz_p> valApodd(INIT_SIZE, m, n);
    Mat<zz_p> valCp(INIT_SIZE, n, p);
    Mat<zz_p> valCpeven(INIT_SIZE, n, p);
    Mat<zz_p> valCpodd(INIT_SIZE, n, p);

    for (long i = 0; i < nb_points; ++i)
    {
        zz_p pt(i+1);
        mul(valAodd[i], pt, valAodd[i]);
        mul(valCodd[i], pt, valCodd[i]);
    }

    zz_p inv2; inv(inv2, to_zz_p(2));
    for (long i = 0; i < nb_points; ++i)
    {
        ell = 0;
        for (long u = 0; u < m; ++u)
            for (long v = 0; v < n; ++v, ++ell)
            {
                add(valApeven[u][v], valAeven[i][ell], valAodd[i][ell]);
                sub(valApodd[u][v], valAeven[i][ell], valAodd[i][ell]);
            }

        ell = 0;
        for (long u = 0; u < n; ++u)
            for (long v = 0; v < p; ++v, ++ell)
            {
                add(valCpeven[u][v], valCeven[i][ell], valCodd[i][ell]);
                sub(valCpodd[u][v], valCeven[i][ell], valCodd[i][ell]);
            }

        mul(valBpeven, valApeven, valCpeven);
        mul(valBpodd, valApodd, valCpodd);

        ell = 0;
        zz_p inv2pt; inv(inv2pt, to_zz_p(2*(i+1)));
        for (long u = 0; u < m; ++u)
            for (long v = 0; v < p; ++v, ++ell)
            {
                add(valBeven[i][ell], valBpeven[u][v], valBpodd[u][v]);
                sub(valBodd[i][ell], valBpeven[u][v], valBpodd[u][v]);
            }
        mul(valBeven[i], valBeven[i], inv2);
        mul(valBodd[i], valBodd[i], inv2pt);
    }

    mul(tmp_mat_even, vB, valBeven);
    mul(tmp_mat_odd, vB, valBodd);

    b.SetDims(m, p);
    ell = 0;
    for (long u = 0; u < m; ++u)
        for (long v = 0; v < p; ++v, ++ell)
        {
            b[u][v].SetLength(dB + 1);
            for (long i = 0; i < dB/2; ++i)
            {
                b[u][v][2*i] = tmp_mat_even[i][ell];
                b[u][v][2*i+1] = tmp_mat_odd[i][ell];
            }
            if (dB % 2 == 0)
                b[u][v][dB] = tmp_mat_even[dB/2][ell];
            b[u][v].normalize();
        }
}
#endif

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
