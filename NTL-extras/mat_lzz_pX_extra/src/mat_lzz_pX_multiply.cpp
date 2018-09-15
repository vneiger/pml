#include <NTL/matrix.h>
#include <NTL/mat_lzz_p.h>
#include <NTL/lzz_pX.h>
#include <NTL/BasicThreadPool.h>

#include "lzz_p_extra.h"
#include "mat_lzz_pX_extra.h"
#include "lzz_pX_CRT.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* c = a*b                                                    */
/* assumes FFT prime and p large enough                       */
/* output may alias input; c does not have to be zero matrix  */
/*------------------------------------------------------------*/
void multiply_evaluate_FFT(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b)
{
    zz_pContext context;
    long s = a.NumRows();
    long t = a.NumCols();
    long u = b.NumCols();

    long dA = deg(a);
    long dB = deg(b);

    long idxk = NextPowerOfTwo(dA + dB + 1);
    fftRep R1(INIT_SIZE, idxk);
    long n = 1 << idxk;
    

    Vec<zz_p> mat_valA, mat_valB;
    Vec<Vec<zz_p>> mat_valC;

    mat_valA.SetLength(n * s * t);
    mat_valB.SetLength(n * t * u);

    long st = s*t;
    
    context.save(); // to give the zz_p context to each thread

NTL_EXEC_RANGE(s,first,last)
    context.restore(); // now all threads have the right zz_p context
    
    fftRep R(INIT_SIZE, idxk);
    for (long i = first; i < last; i++)
    {
        for (long k = 0; k < t; k++)
        {
            TofftRep(R, a[i][k], idxk);
            long *frept = & R.tbl[0][0];
            for (long r = 0, rst = 0; r < n; r++, rst += st)
                mat_valA[rst + i*t + k] = frept[r];
        }
    }
NTL_EXEC_RANGE_END

long tu = t*u;

NTL_EXEC_RANGE(t,first,last)
    context.restore();

    fftRep R(INIT_SIZE, idxk);
    for (long i = first; i < last; i++)
    {
        for (long k = 0; k < u; k++)
        {
            TofftRep(R, b[i][k], idxk);
            long *frept = & R.tbl[0][0];
            for (long r = 0, rtu = 0; r < n; r++, rtu += tu)
                mat_valB[rtu + i*u + k] = frept[r];
        }
    }
    
    R1 = R;
NTL_EXEC_RANGE_END

    mat_valC.SetLength(s * u);
    for (long i = 0; i < s * u; i++)
        mat_valC[i].SetLength(n);

NTL_EXEC_RANGE(n,first,last)
    context.restore();
    
    Mat<zz_p> va, vb, vc; 
    va.SetDims(s, t);
    vb.SetDims(t, u);

    for (long j = first, jst = st*first, jtu = tu*first; j < last; j++, jst += st, jtu += tu)
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
                vb[i][k] = mat_valB[jtu + i*u + k];
            }
        }

        vc = va * vb;

        for (long i = 0; i < s; i++)
        {
            for (long k = 0; k < u; k++)
            {
                mat_valC[i*u + k][j] = vc[i][k];
            }
        }
    }
NTL_EXEC_RANGE_END

    c.SetDims(s, u);

NTL_EXEC_RANGE(s,first,last)
    context.restore();
    
    Mat<zz_p> vc; 
    
    fftRep R = R1;
    
    for (long i = first; i < last; i++)
    {
        for (long k = 0; k < u; k++)
        {
            long *frept = & R.tbl[0][0];
            for (long r = 0; r < n; r++)
            {
                frept[r] = rep(mat_valC[i*u + k][r]);
            }
            FromfftRep(c[i][k], R, 0, n - 1);
        }
    }
NTL_EXEC_RANGE_END
}

// NEW CODE TO BE BENCHMARKED
/*------------------------------------------------------------*/
/* c = a*b                                                    */
/* assumes FFT prime and p large enough                       */
/* output may alias input; c does not have to be zero matrix  */
/*------------------------------------------------------------*/
void multiply_evaluate_FFT2(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b)
{
    zz_pContext context;
    long s = a.NumRows();
    long t = a.NumCols();
    long u = b.NumCols();

    long dA = deg(a);
    long dB = deg(b);

    long idxk = NextPowerOfTwo(dA + dB + 1);
    fftRep R1(INIT_SIZE, idxk);
    long n = 1 << idxk;
    

    Mat<Vec<zz_p>> mat_valA, mat_valB;
    Vec<Vec<zz_p>> mat_valC;

    mat_valA.SetDims(s,t);
    mat_valB.SetDims(t,u);

    long st = s*t;
    
    context.save(); // to give the zz_p context to each thread

NTL_EXEC_RANGE(s,first,last)
    context.restore(); // now all threads have the right zz_p context

    fftRep R(INIT_SIZE, idxk);
    for (long i = first; i < last; i++)
    {
        for (long k = 0; k < t; k++)
        {
            TofftRep(R, a[i][k], idxk);
            mat_valA[i][k].SetLength(n);
            long *frept = & R.tbl[0][0];
            for (long r = 0; r < n; r++)
                mat_valA[i][k][r] = frept[r];
        }
    }
NTL_EXEC_RANGE_END

long tu = t*u;

NTL_EXEC_RANGE(t,first,last)
    context.restore();

    fftRep R(INIT_SIZE, idxk);
    for (long i = first; i < last; i++)
    {
        for (long k = 0; k < u; k++)
        {
            TofftRep(R, b[i][k], idxk);
            mat_valB[i][k].SetLength(n);
            long *frept = & R.tbl[0][0];
            for (long r = 0; r < n; r++)
                mat_valB[i][k][r] = frept[r];
        }
    }
    
    R1 = R;
NTL_EXEC_RANGE_END

    mat_valC.SetLength(s * u);
    for (long i = 0; i < s * u; i++)
        mat_valC[i].SetLength(n);

NTL_EXEC_RANGE(n,first,last)
    context.restore();

    Mat<zz_p> va, vb, vc; 
    va.SetDims(s, t);
    vb.SetDims(t, u);

    for (long j = first, jst = st*first, jtu = tu*first; j < last; j++, jst += st, jtu += tu)
    {
        for (long i = 0; i < s; i++)
        {
            for (long k = 0; k < t; k++)
            {
                va[i][k] = mat_valA[i][k][j];
            }
        }
        for (long i = 0; i < t; i++)
        {
            for (long k = 0; k < u; k++)
            {
                vb[i][k] = mat_valB[i][k][j];
            }
        }

        vc = va * vb;

        for (long i = 0; i < s; i++)
        {
            for (long k = 0; k < u; k++)
            {
                mat_valC[i*u + k][j] = vc[i][k];
            }
        }
    }
NTL_EXEC_RANGE_END

    c.SetDims(s, u);

NTL_EXEC_RANGE(s,first,last)
    context.restore();
    
    Mat<zz_p> vc; 
    
    fftRep R = R1;
    
    for (long i = first; i < last; i++)
    {
        for (long k = 0; k < u; k++)
        {
            long *frept = & R.tbl[0][0];
            for (long r = 0; r < n; r++)
            {
                frept[r] = rep(mat_valC[i*u + k][r]);
            }
            FromfftRep(c[i][k], R, 0, n - 1);
        }
    }
NTL_EXEC_RANGE_END
}
// END NEW CODE TO BE BENCHMARKED


/*------------------------------------------------------------*/
/* c = a*b                                                    */
/* chooses the kind of points                                 */
/* assumes that the field is large enough                     */
/* output may alias input; c does not have to be zero matrix  */
/*------------------------------------------------------------*/
void multiply_evaluate(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b)
{
    
    if (is_FFT_ready(NextPowerOfTwo(deg(a) + deg(b) + 1)))
    {
        multiply_evaluate_FFT(c, a, b);
    }
    else
    {
        multiply_evaluate_geometric(c, a, b);
    }
}

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




// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
