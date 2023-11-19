#include "mat_lzz_pX_multiply.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* c = a*b                                                    */
/* Waksman's algorithm                                        */
/* output may alias input; c does not have to be zero matrix  */
/*------------------------------------------------------------*/
void multiply_waksman(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b)
{
    const long n = a.NumCols();

    if (zz_p::modulus() == 2 || n < 2)
    {
        multiply_naive(c, a, b);
        return;
    }

    if (&c == &a || &c == &b)
    {
        Mat<zz_pX> c2;
        multiply_waksman(c2, a, b);
        c.swap(c2);
        return;
    }

    Vec<zz_pX> d, e;
    zz_pX buf0, buf1, buf2, buf3, buf4;

    const long m = a.NumRows();
    const long p = b.NumCols();

    c.SetDims(m, p);

    d.SetLength(p);
    e.SetLength(m);

    // first iteration j==1 (and j2==1)
    // --> sets initial value for each entry of c
    // (thus, no need to clear its values even if nonzero in input)
    for (long k=0; k<p; ++k)
    {
        // c[0][k] += (a[0][j2-1]+b[j2][k]) * (a[0][j2]+b[j2-1][k]);
        // d[k] += (a[0][j2-1]-b[j2][k]) * (a[0][j2]-b[j2-1][k]);
        add(buf0, a[0][0], b[1][k]);
        sub(buf2, a[0][0], b[1][k]);
        add(buf1, a[0][1], b[0][k]);
        sub(buf3, a[0][1], b[0][k]);
        mul(c[0][k], buf0, buf1);
        mul(d[k], buf2, buf3);
    }

    for (long l=1; l<m; ++l)
    {
        // c[l][0] += (a[l][j2-1]+b[j2][0]) * (a[l][j2]+b[j2-1][0]);
        // e[l] += (a[l][j2-1]-b[j2][0]) * (a[l][j2]-b[j2-1][0]);
        add(buf0, a[l][0], b[1][0]);
        sub(buf2, a[l][0], b[1][0]);
        add(buf1, a[l][1], b[0][0]);
        sub(buf3, a[l][1], b[0][0]);
        mul(c[l][0], buf0, buf1);
        mul(e[l], buf2, buf3);
    }

    for (long l=1; l<m; ++l)
        for (long k=1; k<p; ++k)
        {
            add(buf0, a[l][0], b[1][k]);
            add(buf1, a[l][1], b[0][k]);
            mul(c[l][k], buf0, buf1);
        }

    // other iterations, j>1
    for (long j=2; j<=n/2; ++j)
    {
        const long j2=2*j-1;

        for (long k=0; k<p; ++k)
        {
            // c[0][k] += (a[0][j2-1]+b[j2][k]) * (a[0][j2]+b[j2-1][k]);
            // d[k] += (a[0][j2-1]-b[j2][k]) * (a[0][j2]-b[j2-1][k]);
            add(buf0, a[0][j2-1], b[j2][k]);
            sub(buf2, a[0][j2-1], b[j2][k]);
            add(buf1, a[0][j2], b[j2-1][k]);
            sub(buf3, a[0][j2], b[j2-1][k]);
            mul(buf4, buf0, buf1);
            add(c[0][k], c[0][k], buf4);
            mul(buf4, buf2, buf3);
            add(d[k], d[k], buf4);
        }

        for (long l=1; l<m; ++l)
        {
            // c[l][0] += (a[l][j2-1]+b[j2][0]) * (a[l][j2]+b[j2-1][0]);
            // e[l] += (a[l][j2-1]-b[j2][0]) * (a[l][j2]-b[j2-1][0]);
            add(buf0, a[l][j2-1], b[j2][0]);
            sub(buf2, a[l][j2-1], b[j2][0]);
            add(buf1, a[l][j2], b[j2-1][0]);
            sub(buf3, a[l][j2], b[j2-1][0]);
            mul(buf4, buf0, buf1);
            add(c[l][0], c[l][0], buf4);
            mul(buf4, buf2, buf3);
            add(e[l], e[l], buf4);
        }

        for (long l=1; l<m; ++l)
            for (long k=1; k<p; ++k)
            {
                add(buf0, a[l][j2-1], b[j2][k]);
                add(buf1, a[l][j2], b[j2-1][k]);
                mul(buf2, buf0, buf1);
                add(c[l][k], c[l][k], buf2);
            }
    }

    zz_p inv2;
    inv(inv2, to_zz_p(2));
    for (long l=1; l<m; ++l)
    {
        // e[l] = (e[l]+c[l][0])/2;
        add(e[l], e[l], c[l][0]);
        mul(e[l], e[l], inv2);
        sub(c[l][0], c[l][0], e[l]);
    }

    // buf0 = (e[0]+c[0][0])/2;
    add(buf0, d[0], c[0][0]);
    mul(buf0, buf0, inv2);
    sub(c[0][0], c[0][0], buf0);
    for (long k=1; k<p; ++k)
    {
        // buf1 = (d[k]+c[0][k])/2;
        add(buf1, d[k], c[0][k]);
        mul(buf1, buf1, inv2);
        sub(c[0][k], c[0][k], buf1);
        sub(buf1, buf1, buf0);
        for (long l=1; l<m; ++l)
        {
            // c[l][k] -= buf1 + e[l];
            add(buf2, buf1, e[l]);
            sub(c[l][k], c[l][k], buf2);
        }
    }

    if (n%2) // n is odd
        for (long l=0; l<m; ++l)
            for (long k=0; k<p; ++k)
            {
                // c[l][k] += a[l][n-1]*b[n-1][k];
                mul(buf0, a[l][n-1], b[n-1][k]);
                add(c[l][k], c[l][k], buf0);
            }
}


/*------------------------------------------------------------*/
/* c = a*b                                                    */
/* naive algorithm                                            */
/* output may alias input; c does not have to be zero matrix  */
/*------------------------------------------------------------*/
void multiply_naive(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b)
{
    if (&c == &a || &c == &b)
    {
        Mat<zz_pX> c2;
        multiply_naive(c2, a, b);
        c.swap(c2);
        return;
    }

    const long u = a.NumRows();
    const long v = a.NumCols();
    const long w = b.NumCols();

    c.SetDims(u, w);
    zz_pX buf;
    for (long i = 0; i < u; ++i)
        for (long j = 0; j < w; ++j)
        {
            mul(c[i][j], a[i][0], b[0][j]);
            for (long k = 1; k < v; ++k)
            {
                mul(buf, a[i][k], b[k][j]);
                add(c[i][j], c[i][j], buf);
            }
        }
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
