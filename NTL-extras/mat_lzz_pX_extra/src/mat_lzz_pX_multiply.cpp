#include <NTL/matrix.h>
#include <NTL/mat_lzz_p.h>
#include <NTL/lzz_pX.h>

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

    long deg_trs = max_degree_transform();

    if (dmax <= deg_trs)
    {
        multiply_transform(c, a, b, dmax + 1);
        return;
    }

    // only calibrated for square matrices; here's a hack
    long sz = (a.NumRows() + a.NumCols() + b.NumCols()) / 3;
    long deg_wak = max_degree_waksman(sz);

    if (dmax <= deg_wak)
    {
        multiply_waksman(c, a, b);
        return;
    }

    long p = zz_p::modulus();
    long deg_ev = max_degree_evaluate(sz);
    if (is_prime && p > 2 * (dA + dB + 1) && dmax <= deg_ev)
    {
        multiply_evaluate(c, a, b);
        return;
    }
    else
    {
        multiply_3_primes(c, a, b);
        return;
    }
}

/*------------------------------------------------------------*/
/* transpose of b mapsto c = a*b. output is                   */
/*    trunc( rev(a, dA)*c div x^dA, dB+1 )                    */
/* a must have degree at most dA                              */
/* c must have degree at most dA + dB                         */
/*------------------------------------------------------------*/
void t_multiply(Mat<zz_pX> & b, const Mat<zz_pX> & a, const Mat<zz_pX> & c, long dA, long dB, long is_prime)
{
    long dmax = max(dA, dB);
    long p = zz_p::modulus();
    long sz = (a.NumRows() + a.NumCols() + c.NumCols()) / 3;
    long deg_ev = max_degree_evaluate(sz);

    if (is_prime && p > 2 * (dA + dB + 1) && dmax <= deg_ev)
    {
        t_multiply_evaluate(b, a, c, dA, dB);
        return;
    }
    else
    {
        t_multiply_3_primes(b, a, c, dA, dB);
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
