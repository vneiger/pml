#include "lzz_p_extra.h"
#include "thresholds_matrix_multiply.h"
#include "thresholds_matrix_middle_product.h"
#include "lzz_pX_middle_product.h"
#include "mat_lzz_pX_multiply.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* returns trunc( trunc(a, dA+1)*c div x^dA, dB+1 )           */
/* assumes FFT prime and p large enough                       */
/* output may alias input; b does not have to be zero matrix  */
/*------------------------------------------------------------*/
void middle_product(Mat<zz_pX> & b, const Mat<zz_pX> & a, const Mat<zz_pX> & c, long dA, long dB, long is_prime)
{
    long dmax = max(dA, dB);
    long p = zz_p::modulus();
    long sz = (long) cbrt(a.NumRows() * a.NumCols() * c.NumCols());

    long deg_naive = max_degree_mp_naive(sz);

    if (dmax < deg_naive)
    {
        middle_product_naive(b, a, c, dA, dB);
        return;
    }

    if (is_FFT_ready(NextPowerOfTwo(dA + dB + 1)))
    {
        middle_product_evaluate_FFT(b, a, c, dA, dB);
        return;
    }

    long deg_ev = max_degree_evaluate(sz);
    if (is_prime && p > 2 * (dA + dB + 1) && dmax <= deg_ev)
    {
        middle_product_evaluate_dense(b, a, c, dA, dB);
        return;
    }
    else
    {
        middle_product_3_primes(b, a, c, dA, dB);
        return;
    }

}

/*------------------------------------------------------------*/
/* returns trunc( trunc(a, dA+1)*c div x^dA, dB+1 )           */
/* naive algorithm                                            */
/* output may alias input; b does not have to be zero matrix  */
/*------------------------------------------------------------*/
void middle_product_naive(Mat<zz_pX> & b, const Mat<zz_pX> & a, const Mat<zz_pX> & c, long dA, long dB)
{
    if (&b == &a || &b == &c)
    {
        Mat<zz_pX> b2;
        middle_product_naive(b2, a, c, dA, dB);
        b.swap(b2);
        return;
    }

    const long u = a.NumRows();
    const long v = a.NumCols();
    const long w = c.NumCols();

    b.SetDims(u, w);
    zz_pX buf;
    for (long i = 0; i < u; ++i)
        for (long j = 0; j < w; ++j)
        {
            middle_product(b[i][j], a[i][0], c[0][j], dA, dB);
            for (long k = 1; k < v; ++k)
            {
                middle_product(buf, a[i][k], c[k][j], dA, dB);
                add(b[i][j], b[i][j], buf);
            }
        }
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
