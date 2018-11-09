#include <NTL/matrix.h>
#include <NTL/mat_lzz_p.h>
#include <NTL/lzz_pX.h>
#include <NTL/BasicThreadPool.h>

#include "util.h"
#include "lzz_p_extra.h"
#include "mat_lzz_pX_extra.h"
#include "lzz_pX_CRT.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* c = a*b                                                    */
/*------------------------------------------------------------*/
void multiply(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b, long is_prime)
{
    long dA = deg(a);
    long dB = deg(b);
    long dmax = max(dA, dB);

    if (dA < 0 || dB < 0)
    {
        c.SetDims(a.NumRows(), b.NumCols());
        clear(c);
        return;
    }

    long deg_trs = max_degree_transform();

    if (dmax <= deg_trs)
    {
        multiply_transform(c, a, b, dmax + 1);
        return;
    }

    // only calibrated for square matrices; here's a hack
    long sz = (long) cbrt(a.NumRows() * a.NumCols() * b.NumCols());
    long deg_wak = max_degree_waksman(sz);

    if (dmax <= deg_wak)
    {
        multiply_waksman(c, a, b);
        return;
    }

    if (is_FFT_ready(NextPowerOfTwo(dA + dB + 1)))
    {
        multiply_evaluate_FFT(c, a, b);
        return;
    }

    long p = zz_p::modulus();
    long deg_ev = max_degree_evaluate(sz);
    if (is_prime && p > 2 * (dA + dB + 1) && dmax <= deg_ev)
    {
        multiply_evaluate_dense(c, a, b);
        return;
    }
    else
    {
        multiply_3_primes(c, a, b);
        return;
    }
}

/*------------------------------------------------------------*/
/* multiply by a vector                                       */
/*------------------------------------------------------------*/
void multiply(Vec<zz_pX>& c, const Mat<zz_pX>& a, const Vec<zz_pX>& b, long is_prime)
{
    Mat<zz_pX> cmat, bmat;
    bmat.SetDims(b.length(), 1);
    for (long i = 0; i < b.length(); i++)
        bmat[i][0] = b[i];
    multiply(cmat, a, bmat, is_prime);
    c.SetLength(cmat.NumRows());
    for (long i = 0; i < cmat.NumRows(); i++)
        c[i] = cmat[i][0];
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
