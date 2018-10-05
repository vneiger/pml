#include <memory>
#include <NTL/matrix.h>
#include <NTL/mat_lzz_p.h>
#include <NTL/lzz_pX.h>

#include "lzz_p_extra.h"
#include "mat_lzz_pX_extra.h"
#include "lzz_pX_CRT.h"
#include "thresholds_lmultiplier.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/*             BOILERPLATE FOR LMULTIPLIERS                   */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* current number of rows                                     */
/*------------------------------------------------------------*/
long mat_lzz_pX_lmultiplier::NumRows() const
{
    return __s;
}

/*------------------------------------------------------------*/
/* current number of columns                                  */
/*------------------------------------------------------------*/
long mat_lzz_pX_lmultiplier::NumCols() const
{
    return __t;
}

/*------------------------------------------------------------*/
/* degree of current matrix                                   */
/*------------------------------------------------------------*/
long mat_lzz_pX_lmultiplier::degA() const
{
    return __dA;
}

/*------------------------------------------------------------*/
/* max degree of rhs argument                                 */
/*------------------------------------------------------------*/
long mat_lzz_pX_lmultiplier::degB() const
{
    return __dB;
}

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/*                   FFT-BASED LMULTIPLIERS                   */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* constructor: computes FFT of a                             */
/*------------------------------------------------------------*/
mat_lzz_pX_lmultiplier_FFT::mat_lzz_pX_lmultiplier_FFT(const Mat<zz_pX> & a, long dB)
{
    __s = a.NumRows();
    __t = a.NumCols();
    __dA = deg(a);
    __dB = dB;
    
    idxk = NextPowerOfTwo(__dA + __dB + 1);
    long n = 1L << idxk;
    fftRep R1(INIT_SIZE, idxk);

    va.SetLength(n);
    for (long i = 0; i < n; i++)
        va[i].SetDims(__s, __t);
    
    
    long st = __s * __t;
    for (long i = 0; i < __s; i++)
    {
        for (long k = 0; k < __t; k++)
        {
            TofftRep(R1, a[i][k], idxk);
            long *frept = & R1.tbl[0][0];
            for (long r = 0, rst = 0; r < n; r++, rst += st)
                va[r][i][k] = frept[r];
        }
    }

}

/*------------------------------------------------------------*/
/* right multiplication                                       */
/*------------------------------------------------------------*/
void mat_lzz_pX_lmultiplier_FFT::multiply(Mat<zz_pX>& c, const Mat<zz_pX>& b)
{
    long s = NumRows();
    long t = NumCols();
    long u = b.NumCols();

    long dB = deg(b);

    if (dB > degB())
    {
        LogicError("Rhs degree too large in multiplier");
    }

    fftRep R1(INIT_SIZE, idxk);
    long n = 1L << idxk;
    Vec<Vec<zz_p>> mat_valC;
    
    Vec<zz_p> mat_valB;
    mat_valB.SetLength(n * t * u);

    long st = s*t;
    long tu = t*u;
    for (long i = 0; i < t; i++)
    {
        for (long k = 0; k < u; k++)
        {
            TofftRep(R1, b[i][k], idxk);
            long *frept = & R1.tbl[0][0];
            for (long r = 0, rtu = 0; r < n; r++, rtu += tu)
                mat_valB[rtu + i*u + k] = frept[r];
        }
    }

    Mat<zz_p> vb, vc;
    vb.SetDims(t, u);
    
    mat_valC.SetLength(s * u);
    for (long i = 0; i < s * u; i++)
        mat_valC[i].SetLength(n);
    
    for (long j = 0, jst = 0, jtu = 0; j < n; j++, jst += st, jtu += tu)
    {
        for (long i = 0; i < t; i++)
            for (long k = 0; k < u; k++)
                vb[i][k] = mat_valB[jtu + i*u + k];
        
        vc = va[j] * vb;

        for (long i = 0; i < s; i++)
            for (long k = 0; k < u; k++)
                mat_valC[i*u + k][j] = vc[i][k];
    }
    
    c.SetDims(s, u);
    for (long i = 0; i < s; i++)
        for (long k = 0; k < u; k++)
        {
            long *frept = & R1.tbl[0][0];
            for (long r = 0; r < n; r++)
                frept[r] = rep(mat_valC[i*u + k][r]);
            FromfftRep(c[i][k], R1, 0, n - 1);
        }
}

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/*              GEOMETRIC POINTS LMULTIPLIERS                 */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* constructor: computes evaluation of a                      */
/*------------------------------------------------------------*/
mat_lzz_pX_lmultiplier_geometric::mat_lzz_pX_lmultiplier_geometric(const Mat<zz_pX> & a, long dB)
{
    __s = a.NumRows();
    __t = a.NumCols();
    __dA = deg(a);
    __dB = dB;

    long n = __dA + __dB + 1;
    ev = get_geometric_points(n);
    ev.prepare_degree(__dB);

    va.SetLength(n);
    for (long i = 0; i < n; i++)
        va[i].SetDims(__s, __t);
    
    long st = __s * __t;
    for (long i = 0; i < __s; i++)
    {
        for (long k = 0; k < __t; k++)
        {
            Vec<zz_p> tmp;
            ev.evaluate(tmp, a[i][k]);
            for (long r = 0, rst = 0; r < n; r++, rst += st)
                va[r][i][k] = tmp[r];
        }
    }
}

/*------------------------------------------------------------*/
/* right multiplication                                       */
/*------------------------------------------------------------*/
void mat_lzz_pX_lmultiplier_geometric::multiply(Mat<zz_pX>& c, const Mat<zz_pX>& b) 
{
    long s = NumRows();
    long t = NumCols();
    long u = b.NumCols();

    long dB = deg(b);

    if (dB > degB())
    {
        LogicError("Rhs degree too large in multiplier");
    }

    long n = ev.length();
    Vec<zz_p> mat_valB;
    Vec<Vec<zz_p>> mat_valC;
    
    mat_valB.SetLength(n * t * u);
    long st = s*t;
    long tu = t*u;
    for (long i = 0; i < t; i++)
    {
        for (long k = 0; k < u; k++)
        {
            Vec<zz_p> tmp;
            ev.evaluate(tmp, b[i][k]);
            for (long r = 0, rtu = 0; r < n; r++, rtu += tu)
                mat_valB[rtu + i*u + k] = tmp[r];
        }
    }

    Mat<zz_p> vb, vc;
    vb.SetDims(t, u);
    
    mat_valC.SetLength(s * u);
    for (long i = 0; i < s * u; i++)
        mat_valC[i].SetLength(n);
    
    for (long j = 0, jst = 0, jtu = 0; j < n; j++, jst += st, jtu += tu)
    {
        for (long i = 0; i < t; i++)
            for (long k = 0; k < u; k++)
                vb[i][k] = mat_valB[jtu + i*u + k];
        
        vc = va[j] * vb;
        
        for (long i = 0; i < s; i++)
            for (long k = 0; k < u; k++)
                mat_valC[i*u + k][j] = vc[i][k];
    }
    
    c.SetDims(s, u);
    for (long i = 0; i < s; i++)
        for (long k = 0; k < u; k++)
            ev.interpolate(c[i][k], mat_valC[i*u + k]);
    
}

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/*                      3 PRIMES LMULTIPLIERS                 */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* constructor: builds up to 3 FFT multipliers                */
/*------------------------------------------------------------*/
mat_lzz_pX_lmultiplier_3_primes::mat_lzz_pX_lmultiplier_3_primes(const Mat<zz_pX> & a, long dB)
{
    primes = lzz_pX_3_primes(a.NumCols(), deg(a), dB);
    long nb = primes.nb();
    FFT_muls.SetLength(nb);
    for (long i = 0; i < nb; i++)
    {
        zz_pPush push;
        zz_p::FFTInit(i);
        Mat<zz_pX> ap = a;
        reduce_mod_p(ap);
        FFT_muls[i] = mat_lzz_pX_lmultiplier_FFT(ap, dB);
    }
}

/*------------------------------------------------------------*/
/* right multiplication                                       */
/*------------------------------------------------------------*/
void mat_lzz_pX_lmultiplier_3_primes::multiply(Mat<zz_pX>& c, const Mat<zz_pX>& b) 
{
    long nb = primes.nb();
    Vec<Mat<zz_pX>> cs;
    cs.SetLength(nb);

    long p = zz_p::modulus();

    for (long i = 0; i < nb; i++)
    {
        zz_pPush push;
        zz_p::FFTInit(i);
        long pi = zz_p::modulus();

        if (p < pi)
            FFT_muls[i].multiply(cs[i], b);
        else
        {
            Mat<zz_pX> bi = b;
            reduce_mod_p(bi);
            FFT_muls[i].multiply(cs[i], bi);
        }        
    }

    if (nb == 1)
    {
        c = cs[0];
        reduce_mod_p(c);
    }
    else
    {
        primes.reconstruct(c, cs);
    }
}

/*------------------------------------------------------------*/
/* returns a multiplier of the right type                     */
/*------------------------------------------------------------*/
std::unique_ptr<mat_lzz_pX_lmultiplier> get_lmultiplier(const Mat<zz_pX> & a, long dB)
{
    long t = type_of_prime();
    
    if (t == TYPE_FFT_PRIME)
        return std::unique_ptr<mat_lzz_pX_lmultiplier>(new mat_lzz_pX_lmultiplier_FFT(a, dB));

    if ( (t == TYPE_LARGE_PRIME && deg(a) <= MATRIX_LMULTIPLIER_DEGREE_THRESHOLD_LARGE)
         || 
         (t == TYPE_SMALL_PRIME && deg(a) <= MATRIX_LMULTIPLIER_DEGREE_THRESHOLD_SMALL) )
        return std::unique_ptr<mat_lzz_pX_lmultiplier>(new mat_lzz_pX_lmultiplier_3_primes(a, dB));
    else
        return std::unique_ptr<mat_lzz_pX_lmultiplier>(new mat_lzz_pX_lmultiplier_geometric(a, dB));
}
