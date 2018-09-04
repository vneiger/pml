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
void middle_product_FFT(Mat<zz_pX> & b, const Mat<zz_pX> & a, const Mat<zz_pX> & c, long dA, long dB)
{
    long s = a.NumRows();
    long t = a.NumCols();
    long u = c.NumCols();

    long idxk = NextPowerOfTwo(dA + dB + 1);
    fftRep R1(INIT_SIZE, idxk);

    long n = 1 << idxk;

    Vec<zz_p> mat_valA, mat_valC;
    Vec<Vec<zz_p>> mat_valB;

    mat_valA.SetLength(n * s * t);
    mat_valC.SetLength(n * t * u);

    Vec<zz_p> tmp;
    long st = s*t;
    for (long i = 0; i < s; i++)
    {
        for (long k = 0; k < t; k++)
        {
            TofftRep(R1, a[i][k], idxk);
            long *frept = & R1.tbl[0][0];
            for (long r = 0, rst = 0; r < n; r++, rst += st)
                mat_valA[rst + i*t + k] = frept[r];
        }
    }

    long tu = t*u;
    for (long i = 0; i < t; i++)
    {
        for (long k = 0; k < u; k++)
        {
            TofftRep(R1, c[i][k], idxk);
            long *frept = & R1.tbl[0][0];
            for (long r = 0, rtu = 0; r < n; r++, rtu += tu)
                mat_valC[rtu + i*u + k] = frept[r];
        }
    }

    Mat<zz_p> va, vb, vc;
    va.SetDims(s, t);
    vc.SetDims(t, u);

    mat_valB.SetLength(s * u);
    for (long i = 0; i < s * u; i++)
        mat_valB[i].SetLength(n);

    for (long j = 0, jst = 0, jtu = 0; j < n; j++, jst += st, jtu += tu)
    {
        for (long i = 0; i < s; i++)
        {
            for (long k = 0; k < t; k++)
            {
                va[i][k] = mat_valA[jst + i*t + k];
            }
        }
        for (long i = 0; i < t; i++)
        {
            for (long k = 0; k < u; k++)
            {
                vc[i][k] = mat_valC[jtu + i*u + k];
            }
        }

        vb = va * vc;

        for (long i = 0; i < s; i++)
        {
            for (long k = 0; k < u; k++)
            {
                mat_valB[i*u + k][j] = vb[i][k];
            }
        }
    }

    b.SetDims(s, u);
    for (long i = 0; i < s; i++)
    {
        for (long k = 0; k < u; k++)
        {
            long *frept = & R1.tbl[0][0];
            for (long r = 0; r < n; r++)
            {
                frept[r] = rep(mat_valB[i*u + k][r]);
            }
            FromfftRep(b[i][k], R1, dA, dA + dB);
        }
    }
}


/*------------------------------------------------------------*/
/* returns trunc( trunc(a, dA+1)*c div x^dA, dB+1 )           */
/* chooses the kind of points                                 */
/* assumes that the field is large enough                     */
/* output may alias input; b does not have to be zero matrix  */
/*------------------------------------------------------------*/
void middle_product_evaluate(Mat<zz_pX> & b, const Mat<zz_pX> & a, const Mat<zz_pX> & c, long dA, long dB)
{
	if (is_FFT_ready(NextPowerOfTwo(dA + dB + 1)))
	{
		middle_product_FFT(b, a, c, dA, dB);
	}
	else
	{
		t_multiply_evaluate_geometric(b, reverse(a, dA), c, dA, dB);
	}
}


/*------------------------------------------------------------*/
/* returns trunc( trunc(a, dA+1)*c div x^dA, dB+1 )           */
/* assumes FFT prime and p large enough                       */
/* output may alias input; b does not have to be zero matrix  */
/*------------------------------------------------------------*/
void middle_product(Mat<zz_pX> & b, const Mat<zz_pX> & a, const Mat<zz_pX> & c, long dA, long dB, long is_prime)
{
    long dmax = max(dA, dB);
    long p = zz_p::modulus();

    long sz = (a.NumRows() + a.NumCols() + c.NumCols()) / 3;
    long deg_naive = max_degree_mp_naive(sz);

    if (dmax <= deg_naive)
    {
        middle_product_naive(b, a, c, dA, dB);
        return;
    }


    long deg_ev = max_degree_mp_evaluate(sz);

    if (is_prime && p > 2 * (dA + dB + 1) && dmax <= deg_ev)
    {
        middle_product_evaluate(b, a, c, dA, dB);
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
