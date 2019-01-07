#include <memory>
#include <NTL/matrix.h>
#include <NTL/mat_lzz_p.h>
#include <NTL/lzz_pX.h>

#include "lzz_p_extra.h"
#include "mat_lzz_pX_extra.h"
#include "lzz_pX_CRT.h"
#include "thresholds_lmultiplier.h"

NTL_CLIENT

#if defined(NTL_HAVE_LL_TYPE) && defined(NTL_HAVE_SP_LL_ROUTINES)

/*------------------------------------------------------------*/
/* pairwise product of two fftReps using long long's          */
/*------------------------------------------------------------*/
static inline void mul(Vec<ll_type>& z, const fftRep& x, const fftRep& y)
{
    const long *xp = &x.tbl[0][0];
    const long *yp = &y.tbl[0][0];
    long len = min(x.len, y.len);
    for (long j = 0; j < len; j++)
        ll_mul(z[j], xp[j], yp[j]);
}

/*------------------------------------------------------------*/
/* pairwise product of two fftReps using long long's          */
/*------------------------------------------------------------*/
static inline void mul_add(Vec<ll_type>& z, const fftRep& x, const fftRep& y)
{
    const long *xp = &x.tbl[0][0];
    const long *yp = &y.tbl[0][0];
    long len = min(x.len, y.len);
    for (long j = 0; j < len; j++)
        ll_mul_add(z[j], xp[j], yp[j]);
}

#endif

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
/*           FFT-BASED LMULTIPLIER (DIRECT)                   */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* constructor: computes FFT of a                             */
/*------------------------------------------------------------*/
mat_lzz_pX_lmultiplier_FFT_direct::mat_lzz_pX_lmultiplier_FFT_direct(const Mat<zz_pX> & a, long dB)
{
    __s = a.NumRows();
    __t = a.NumCols();
    __dA = deg(a);
    __dB = dB;

    long m = __s;
    long n = __t;
    
    K = NextPowerOfTwo(__dA + dB + 1);
    len = 1L << K;

    n0 = (1L << (2*(NTL_BITS_PER_LONG - NTL_SP_NBITS))) - 1;
    first_slice = n % n0;
    nb_slices = n / n0;
    if (first_slice == 0)
    {
        first_slice = n0;
        nb_slices--;
    }

    pr = zz_p::modulus();
    red1 = sp_PrepRem(pr);
    red2 = make_sp_ll_reduce_struct(pr);

    vala.SetLength(m);
    for (long i = 0; i < m; i++)
    {
        vala[i].SetLength(n);
        for (long j = 0; j < n; j++)
            TofftRep(vala[i][j], a[i][j], K);
    }
}

/*------------------------------------------------------------*/
/* right multiplication                                       */
/*------------------------------------------------------------*/
void mat_lzz_pX_lmultiplier_FFT_direct::multiply(Mat<zz_pX>& c, const Mat<zz_pX>& b)
{
    if (&c == &b)
    {
        Mat<zz_pX> c2;
        multiply(c2, b);
        c = c2;
        return;
    }

#if defined(NTL_HAVE_LL_TYPE) && defined(NTL_HAVE_SP_LL_ROUTINES)
    Vec<fftRep> valb;
    long m, n, p;
    Vec<ll_type> tmp;
    fftRep tmp_r;

    m = NumRows();
    n = NumCols();
    p = b.NumCols();

    c.SetDims(m, p);
    tmp.SetLength(len);
    tmp_r = fftRep(INIT_SIZE, K);
    TofftRep(tmp_r, b[0][0], K);

    valb.SetLength(n);
    if (NumBits(pr) == NTL_SP_NBITS) // we can use normalized remainders; may be a bit faster
    {
        for (long i = 0; i < p; i++)
        {
            for (long j = 0; j < n; j++)
                TofftRep(valb[j], b[j][i], K);

            for (long k = 0; k < m; k++)
            {
                fftRep * va = vala[k].elts();

                mul(tmp, valb[0], va[0]);
                for (long j = 1; j < first_slice; j++)
                    mul_add(tmp, valb[j], va[j]);
                
                long start = first_slice;
                for (long jj = 0; jj < nb_slices; jj++)
                {
                    for (long x = 0; x < len; x++)
                    {
                        tmp[x].lo = sp_ll_red_21_normalized(rem(tmp[x].hi, pr, red1), tmp[x].lo, pr, red2);
                        tmp[x].hi = 0;
                    }
                    for (long j = 0; j < n0; j++)
                        mul_add(tmp, valb[j+start], va[j+start]);
                    start += n0;
                }

                long *tmp_ptr = &tmp_r.tbl[0][0];
                for (long x = 0; x < len; x++) 
                    tmp_ptr[x] = sp_ll_red_21_normalized(rem(tmp[x].hi, pr, red1), tmp[x].lo, pr, red2);
                FromfftRep(c[k][i], tmp_r, 0, __dA + __dB);
            }
        }
    }
    else
    {
        for (long i = 0; i < p; i++)
        {
            for (long j = 0; j < n; j++)
                TofftRep(valb[j], b[j][i], K);

            for (long k = 0; k < m; k++)
            {
                fftRep * va = vala[k].elts();

                mul(tmp, valb[0], va[0]);
                for (long j = 1; j < first_slice; j++)
                    mul_add(tmp, valb[j], va[j]);
                
                long start = first_slice;
                for (long jj = 0; jj < nb_slices; jj++)
                {
                    for (long x = 0; x < len; x++)
                    {
                        tmp[x].lo = sp_ll_red_21(rem(tmp[x].hi, pr, red1), tmp[x].lo, pr, red2);
                        tmp[x].hi = 0;
                    }
                    for (long j = 0; j < n0; j++)
                        mul_add(tmp, valb[j+start], va[j+start]);
                    start += n0;
                }

                long *tmp_ptr = &tmp_r.tbl[0][0];
                for (long x = 0; x < len; x++) 
                    tmp_ptr[x] = sp_ll_red_21(rem(tmp[x].hi, pr, red1), tmp[x].lo, pr, red2);
                FromfftRep(c[k][i], tmp_r, 0, __dA + __dB);
            }
        }
    }
#else
    Vec<fftRep> valb;

    long m = NumRows();
    long n = NumCols();
    long p = b.NumCols();

    c.SetDims(m, p);
    valb.SetLength(n);
    fftRep tmp1 = fftRep(INIT_SIZE, K);
    fftRep tmp2 = fftRep(INIT_SIZE, K);
    for (long i = 0; i < p; i++)
    {
        for (long j = 0; j < n; j++)
            TofftRep(valb[j], b[j][i], K);
        for (long k = 0; k < m; k++)
        {
            fftRep * va = vala[k].elts();
            mul(tmp1, valb[0], va[0]);
            for (long j = 1; j < n; j++)
            {
                mul(tmp2, valb[j], va[j]);
                add(tmp1, tmp1, tmp2);
            }
            FromfftRep(c[k][i], tmp1, 0, __dA + __dB);
        }
    }
#endif
}


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/*            FFT-BASED LMULTIPLIERS (MATMUL)                 */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* constructor: computes FFT of a                             */
/*------------------------------------------------------------*/
mat_lzz_pX_lmultiplier_FFT_matmul::mat_lzz_pX_lmultiplier_FFT_matmul(const Mat<zz_pX> & a, long dB)
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
/* uses Mat<zz_p> matrix multiplication                       */
/*------------------------------------------------------------*/
void mat_lzz_pX_lmultiplier_FFT_matmul::multiply(Mat<zz_pX>& c, const Mat<zz_pX>& b)
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
/*                DENSE ALGORITHM LMULTIPLIERS                */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* computes vandermonde-s and evaluates a                     */
/*------------------------------------------------------------*/
mat_lzz_pX_lmultiplier_dense::mat_lzz_pX_lmultiplier_dense(const Mat<zz_pX> & a, long dB)
{
    long ell;

    __s = a.NumRows();
    __t = a.NumCols();
    __dA = deg(a);
    __dB = dB;

    vandermonde(vA, vB, iV, __dA, __dB);
    nb_points = vA.NumRows();

    Mat<zz_p> tmp_mat(INIT_SIZE, __dA + 1, __s * __t);
    ell = 0;
    for (long i = 0; i < __s; i++)
        for (long j = 0; j < __t; j++)
        {
            long d = deg(a[i][j]);
            if (d >= 0)
            {
                const zz_p * cAij = a[i][j].rep.elts();
                long k;
                for (k = 0; k <= d; k++)  // k <= d-2 so k+1 <= d-1
                {
                    tmp_mat[k][ell] = cAij[k];
                }
            }
            for (long k = d+1; k <= __dA; k++)
            {
                tmp_mat[k][ell] = 0;
            }
            
            ell++;
        }
    valA = vA * tmp_mat;
}

/*------------------------------------------------------------*/
/* right multiplication                                       */
/*------------------------------------------------------------*/
void mat_lzz_pX_lmultiplier_dense::multiply(Mat<zz_pX>& c, const Mat<zz_pX>& b) 
{
    Mat<zz_p> valAp, valBp, valCp, valC, valB, tmp_mat;
    long ell, m, n, p;
    
    m = NumRows();
    n = NumCols();
    p = b.NumCols();

    tmp_mat.SetDims(__dB + 1, n * p);
    ell = 0;
    for (long i = 0; i < n; i++)
        for (long j = 0; j < p; j++)
        {
            long d = deg(b[i][j]);
            if (d >= 0)
            {
                const zz_p * cBij = b[i][j].rep.elts();
                long k;
                for (k = 0; k <= d; k++)  // k <= d-2 so k+1 <= d-1
                {
                    tmp_mat[k][ell] = cBij[k];
                }
            }
            for (long k = d+1; k <= __dA; k++)
            {
                tmp_mat[k][ell] = 0;
            }
            ell++;
        }
    valB = vB * tmp_mat;

    valAp.SetDims(m, n);
    valBp.SetDims(n, p);
    valC.SetDims(nb_points, m * p);
    for (long i = 0; i < nb_points; i++)
    {
        long ell;
        ell = 0;
        for (long u = 0; u < m; u++)
            for (long v = 0; v < n; v++)
            {
                valAp[u][v] = valA[i][ell++];
            }
        
        ell = 0;
        for (long u = 0; u < n; u++)
            for (long v = 0; v < p; v++)
            {
                valBp[u][v] = valB[i][ell++];
            }
        
        valCp = valAp * valBp;
        
        ell = 0;
        for (long u = 0; u < m; u++)
            for (long v = 0; v < p; v++)
            {
                valC[i][ell++] = valCp[u][v];
            }
    }
    tmp_mat = iV*valC;
    
    c.SetDims(m, p);
    ell = 0;
    for (long u = 0; u < m; u++)
        for (long v = 0; v < p; v++)
        {
            c[u][v].rep.SetLength(nb_points);
            zz_p * cc = c[u][v].rep.elts();
            for (long i = 0; i < nb_points; i++)
                cc[i] = tmp_mat[i][ell];
            c[u][v].normalize();
            ell++;
        }
}

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/*                    3 PRIMES LMULTIPLIERS                   */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* constructor: builds up to 3 FFT multipliers                */
/*------------------------------------------------------------*/
mat_lzz_pX_lmultiplier_3_primes::mat_lzz_pX_lmultiplier_3_primes(const Mat<zz_pX> & a, long dB):
    primes(a.NumCols(), deg(a), dB)
{
    long nb = primes.nb();
    FFT_muls.SetLength(nb);
    for (long i = 0; i < nb; i++)
    {
        zz_pPush push;
        zz_p::FFTInit(i);
        Mat<zz_pX> ap;
        reduce_mod_p(ap, a);
        FFT_muls[i] = get_lmultiplier(ap, dB);
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
        {
            FFT_muls[i]->multiply(cs[i], b);
        }
        else
        {
            Mat<zz_pX> bi;
            reduce_mod_p(bi, b);
            FFT_muls[i]->multiply(cs[i], bi);
        }        
    }

    if (nb == 1)
    {
        reduce_mod_p(c, cs[0]);
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
    {
        if (a.NumRows() <= 100)
            return std::unique_ptr<mat_lzz_pX_lmultiplier>(new mat_lzz_pX_lmultiplier_FFT_direct(a, dB));
        else
            return std::unique_ptr<mat_lzz_pX_lmultiplier>(new mat_lzz_pX_lmultiplier_FFT_matmul(a, dB));
    }
    if (t == TYPE_LARGE_PRIME)
    {
        if (a.NumRows() <= 30)
            return std::unique_ptr<mat_lzz_pX_lmultiplier>(new mat_lzz_pX_lmultiplier_3_primes(a, dB));
        if (deg(a) <= 150)
            return std::unique_ptr<mat_lzz_pX_lmultiplier>(new mat_lzz_pX_lmultiplier_dense(a, dB));
        if (a.NumRows() <= 200)
            return std::unique_ptr<mat_lzz_pX_lmultiplier>(new mat_lzz_pX_lmultiplier_3_primes(a, dB));
        else
            return std::unique_ptr<mat_lzz_pX_lmultiplier>(new mat_lzz_pX_lmultiplier_geometric(a, dB));
    }
    else 
    {
        if (a.NumRows() <= 30)
            return std::unique_ptr<mat_lzz_pX_lmultiplier>(new mat_lzz_pX_lmultiplier_3_primes(a, dB));
        if (deg(a) <= 150)
            return std::unique_ptr<mat_lzz_pX_lmultiplier>(new mat_lzz_pX_lmultiplier_dense(a, dB));
        else
            return std::unique_ptr<mat_lzz_pX_lmultiplier>(new mat_lzz_pX_lmultiplier_3_primes(a, dB));
    }

}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
