#include <NTL/matrix.h>
#include <NTL/mat_lzz_p.h>
#include <NTL/lzz_pX.h>

#include "util.h"
#include "lzz_p_extra.h"
#include "mat_lzz_pX_extra.h"
#include "lzz_pX_CRT.h"

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
    long sz = (long) cbrt(a.NumRows() * a.NumCols() * b.NumCols());

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

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
