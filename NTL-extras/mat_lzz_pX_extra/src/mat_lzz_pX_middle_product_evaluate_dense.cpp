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
        {
            const zz_p * cAij = a[i][j].rep.elts();
            for (long k = 0; k <= deg(a[i][j]); ++k)
                tmp_mat[dA-k][ell] = cAij[k];
        }

    Mat<zz_p> valA;
    mul(valA, vA, tmp_mat);

    // evaluations of matrix c
    tmp_mat.SetDims(nb_points, n * p);
    ell = 0;
    for (long i = 0; i < n; ++i)
        for (long j = 0; j < p; ++j, ++ell)
        {
            long d = deg(c[i][j]);
            if (d >= 0)
            {
                const zz_p * cCij = c[i][j].rep.elts();
                for (long k = 0; k <= d; ++k)
                    tmp_mat[k][ell] = cCij[k];
            }
            for (long k = d+1; k < nb_points; k++)
                tmp_mat[k][ell] = 0;
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
            zz_p * bb = b[u][v].rep.elts();
            for (long i = 0; i < dB + 1; ++i)
                bb[i] = tmp_mat[i][ell];
            b[u][v].normalize();
        }
}


// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
