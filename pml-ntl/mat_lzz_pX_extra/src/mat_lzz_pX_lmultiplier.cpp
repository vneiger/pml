#include <NTL/FFT_impl.h>
#include "lzz_p_extra.h" // type_of_prime
#include "mat_lzz_pX_multiply.h"
#include "thresholds_lmultiplier.h"

// FIXME work in progress:
// constants used to reduce cache misses
#define CACHE_LINE_SIZE 8
#define MATRIX_BLOCK_SIZE 16
// right now these are chosen harcoded for L1 cache line 64B (8 long's) and L1
// total cache 32k --> 512 ~ 16*16 cache lines


NTL_CLIENT

#if defined(NTL_HAVE_LL_TYPE) && defined(NTL_HAVE_SP_LL_ROUTINES) \
     && defined(__GNUC__) && (__GNUC__ >= 4) && !defined(__INTEL_COMPILER)  && !defined(__clang__) \
     && defined (__x86_64__) && NTL_BITS_PER_LONG == 64

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
/*------------------------------------------------------------*/
/*           FFT-BASED LMULTIPLIER (DIRECT)                   */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* constructor: computes FFT of a                             */
/*------------------------------------------------------------*/
mat_lzz_pX_lmultiplier_FFT_direct::mat_lzz_pX_lmultiplier_FFT_direct(const Mat<zz_pX> & a, long dB) :
    mat_lzz_pX_lmultiplier(a, dB)
{
    K = NextPowerOfTwo(__dA + dB + 1);

    n0 = (1L << (2*(NTL_BITS_PER_LONG - NTL_SP_NBITS))) - 1;
    first_slice = (__t % n0 != 0) ? (__t%n0) : n0;
    nb_slices = (__t % n0 != 0) ? (__t / n0) : (__t/n0 - 1);

    pr = zz_p::modulus();
    red1 = sp_PrepRem(pr);
    red2 = make_sp_ll_reduce_struct(pr);

    const long len = FFTRoundUp(__dA + dB + 1, K);
    vala.SetLength(__s);
    for (long i = 0; i < __s; ++i)
    {
        vala[i].SetLength(__t);
        for (long j = 0; j < __t; ++j)
            TofftRep_trunc(vala[i][j], a[i][j], K, len);
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
        c.swap(c2);
        return;
    }

    const long dB = deg(b);
    if (dB > __dB)
        LogicError("Rhs degree too large in lmultiplier");

    const long cube_dim = NumRows() * NumCols() * b.NumCols();
    const long d = (__dA+dB)/2;
    if (NumBits(zz_p::modulus()) < SMALL_PRIME_SIZE)
    {
        if (cube_dim < 4*4*4 || (cube_dim <= 6*6*6 && d > 100))
            multiply_direct(c, b);
        else
            multiply_direct_ll_type(c, b);
    }
    else
    {
        // if < 4*4*4, or if <= 8*8*8 under some conditions, use direct
        if (cube_dim < 4*4*4 || (cube_dim == 4*4*4 && d>=90) || (cube_dim <= 8*8*8 && d >= 65))
            multiply_direct(c, b);
        // otherwise use direct_ll
        else
            multiply_direct_ll_type(c, b);
    }
}

// note: c cannot alias b
void mat_lzz_pX_lmultiplier_FFT_direct::multiply_direct_ll_type(Mat<zz_pX>& c, const Mat<zz_pX>& b)
{
#if defined(NTL_HAVE_LL_TYPE) && defined(NTL_HAVE_SP_LL_ROUTINES) \
     && defined(__GNUC__) && (__GNUC__ >= 4) && !defined(__INTEL_COMPILER)  && !defined(__clang__) \
     && defined (__x86_64__) && NTL_BITS_PER_LONG == 64
    const long m = NumRows();
    const long n = NumCols();
    const long p = b.NumCols();

    // actual wanted length for truncated FFT
    const long len_actual = __dA+deg(b)+1;
    const long len = FFTRoundUp(len_actual, K);

    c.SetDims(m, p);
    Vec<fftRep> valb(INIT_SIZE, n);

    fftRep tmp_r(INIT_SIZE, K);
    TofftRep_trunc(tmp_r, b[0][0], K, len);

    Vec<ll_type> tmp(INIT_SIZE, len);
    if (NumBits(pr) == NTL_SP_NBITS) // we can use normalized remainders; may be a bit faster
    {
        for (long i = 0; i < p; i++)
        {
            for (long j = 0; j < n; j++)
                TofftRep_trunc(valb[j], b[j][i], K, len);

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
                FromfftRep(c[k][i], tmp_r, 0, __dA + deg(b));
            }
        }
    }
    else
    {
        for (long i = 0; i < p; i++)
        {
            for (long j = 0; j < n; j++)
                TofftRep_trunc(valb[j], b[j][i], K, len);

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
                FromfftRep(c[k][i], tmp_r, 0, len_actual-1);
            }
        }
    }
#else
    multiply_direct(c, b);
#endif
}

void mat_lzz_pX_lmultiplier_FFT_direct::multiply_direct(Mat<zz_pX>& c, const Mat<zz_pX>& b)
{
    // dimensions
    const long m = NumRows();
    const long n = NumCols();
    const long p = b.NumCols();

    // actual wanted length for truncated FFT
    const long len_actual = __dA+deg(b)+1;
    const long len = FFTRoundUp(len_actual, K);

    c.SetDims(m, p);
    Vec<fftRep> valb(INIT_SIZE, n);
    fftRep tmp1(INIT_SIZE, K);
    fftRep tmp2(INIT_SIZE, K);
    for (long i = 0; i < p; ++i)
    {
        for (long j = 0; j < n; ++j)
            TofftRep_trunc(valb[j], b[j][i], K, len);
        for (long k = 0; k < m; ++k)
        {
            fftRep * va = vala[k].elts();
            mul(tmp1, valb[0], va[0]);
            for (long j = 1; j < n; ++j)
            {
                mul(tmp2, valb[j], va[j]);
                add(tmp1, tmp1, tmp2);
            }
            FromfftRep(c[k][i], tmp1, 0, len_actual-1);
        }
    }
}


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/*            FFT-BASED LMULTIPLIERS (MATMUL)                 */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* constructor: computes FFT of a                             */
/*------------------------------------------------------------*/
mat_lzz_pX_lmultiplier_FFT_matmul::mat_lzz_pX_lmultiplier_FFT_matmul(const Mat<zz_pX> & a, long dB) :
    mat_lzz_pX_lmultiplier(a, dB)
{
    idxk = NextPowerOfTwo(__dA + __dB + 1);

    // initialize fftRep tools
    fftRep R(INIT_SIZE, idxk);
    Vec<UniqueArray<long>> RR(INIT_SIZE, MATRIX_BLOCK_SIZE);

    const long n = 1<<idxk;
    const long len = FFTRoundUp(__dA + __dB + 1, idxk);
    va.SetLength(len);
    for (long i = 0; i < len; ++i)
        va[i].SetDims(__s, __t);

    for (long i = 0; i < __s; ++i)
        for (long k = 0; k < __t; k+=MATRIX_BLOCK_SIZE)
        {
            const long kk_bnd = std::min((long)MATRIX_BLOCK_SIZE, __t-k);
            for (long kk=0; kk<kk_bnd; ++kk)
            {
                R.tbl[0].SetLength(n);
                TofftRep_trunc(R, a[i][k+kk], idxk, len);
                RR[kk].swap(R.tbl[0]);
            }
            for (long r = 0; r < len; ++r)
                for (long kk=0; kk<kk_bnd; ++kk)
                    va[r][i][k+kk].LoopHole() = RR[kk][r];
        }
}

/*------------------------------------------------------------*/
/* right multiplication                                       */
/* uses Mat<zz_p> matrix multiplication                       */
/*------------------------------------------------------------*/
void mat_lzz_pX_lmultiplier_FFT_matmul::multiply(Mat<zz_pX>& c, const Mat<zz_pX>& b)
{
    const long s = NumRows();
    const long t = NumCols();
    const long u = b.NumCols();

    const long dB = deg(b);

    if (dB > degB())
        LogicError("Rhs degree too large in multiplier");

    // actual wanted length for truncated FFT
    const long len = FFTRoundUp(__dA+dB+1, idxk);
    const long n = 1<<idxk;

    // fft representation
    fftRep R(INIT_SIZE, idxk);
    Vec<UniqueArray<long>> RR(INIT_SIZE, CACHE_LINE_SIZE);

    // stores evaluations in a single vector for a, same for b
    // mat_valB[r*t*u + i*t + k] is b[i][k] evaluated at the r-th point
    const long tu = t*u;
    Vec<long> mat_valB(INIT_SIZE, len * tu);
    for (long i = 0; i < t; ++i)
    {
        for (long k = 0; k < u; k+=CACHE_LINE_SIZE)
        {
            const long kk_bnd = std::min((long)CACHE_LINE_SIZE, u-k);
            for (long kk = 0; kk < kk_bnd; ++kk)
            {
                R.tbl[0].SetLength(n);
                TofftRep_trunc(R, b[i][k+kk], idxk, len);
                R.tbl[0].swap(RR[kk]);
            }
            for (long r = 0, rtu = 0; r < len; ++r, rtu += tu)
                for (long kk = 0; kk < kk_bnd; ++kk)
                    mat_valB[rtu + i*u + k+kk] = RR[kk][r];
        }
    }

    // will store a evaluated at the r-th point, same for b
    Mat<zz_p> vb(INIT_SIZE, t, u);
    // will store the product evaluated at the r-th point, i.e. va*vb
    Vec<Mat<zz_p>> vc(INIT_SIZE, CACHE_LINE_SIZE);

    // stores all evaluations of product c = a*b
    Vec<UniqueArray<long>> mat_valC(INIT_SIZE, s*u);
    for (long i = 0; i < s * u; ++i)
        mat_valC[i].SetLength(n);

    // compute pairwise products
    for (long j = 0, jtu = 0; j < len; j+=CACHE_LINE_SIZE)
    {
        const long jj_bnd = std::min((long)CACHE_LINE_SIZE, len-j);
        for (long jj = 0; jj < jj_bnd; ++jj, jtu += tu)
        {
            for (long i = 0; i < t; ++i)
                for (long k = 0; k < u; ++k)
                    vb[i][k].LoopHole() = mat_valB[jtu + i*u + k];
            mul(vc[jj], va[j+jj], vb);
        }

        for (long i = 0; i < s; i++)
            for (long k = 0; k < u; k+=MATRIX_BLOCK_SIZE)
            {
                const long kk_bnd = std::min((long)MATRIX_BLOCK_SIZE, u-k);
                for (long jj = 0; jj < jj_bnd; ++jj)
                    for (long kk = 0; kk < kk_bnd; ++kk)
                        mat_valC[i*u + k+kk][j+jj] = vc[jj][i][k+kk]._zz_p__rep;
            }
    }

    // interpolate the evaluations stored in mat_valC back into c
    c.SetDims(s, u);
    for (long i = 0; i < s; ++i)
        for (long k = 0; k < u; ++k)
        {
            R.tbl[0].swap(mat_valC[i*u+k]);
            FromfftRep(c[i][k], R, 0, __dA+dB);
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
mat_lzz_pX_lmultiplier_geometric::mat_lzz_pX_lmultiplier_geometric(const Mat<zz_pX> & a, long dB) :
    mat_lzz_pX_lmultiplier(a, dB)
{
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
        LogicError("Rhs degree too large in multiplier");

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
mat_lzz_pX_lmultiplier_dense::mat_lzz_pX_lmultiplier_dense(const Mat<zz_pX> & a, long dB) :
    mat_lzz_pX_lmultiplier(a, dB)
{
    Mat<zz_p> vA;
    vandermonde(vA, vB, iV, __dA, __dB);
    nb_points = vA.NumRows();

    // evaluate a: reorganize data
    Mat<zz_p> tmp_mat(INIT_SIZE, __dA + 1, __s * __t);
    long ell = 0;
    for (long i = 0; i < __s; ++i)
        for (long j = 0; j < __t; ++j, ++ell)
            for (long k = 0; k <= deg(a[i][j]); ++k)
                tmp_mat[k][ell] = a[i][j][k];

    // evaluate a: use Mat<zz_p> multiplication
    Mat<zz_p> valA;
    mul(valA, vA, tmp_mat);

    // reorganize data into valAp
    valAp.SetLength(nb_points);
    for (long i = 0; i < nb_points; ++i)
    {
        valAp[i].SetDims(__s, __t);
        ell = 0;
        for (long u = 0; u < __s; ++u)
            for (long v = 0; v < __t; ++v, ++ell)
                valAp[i][u][v] = valA[i][ell];
    }
}

/*------------------------------------------------------------*/
/* right multiplication                                       */
/*------------------------------------------------------------*/
void mat_lzz_pX_lmultiplier_dense::multiply(Mat<zz_pX>& c, const Mat<zz_pX>& b)
{
    const long m = NumRows();
    const long n = NumCols();
    const long p = b.NumCols();

    Mat<zz_p> tmp_mat(INIT_SIZE, __dB + 1, n * p);
    long ell = 0;
    for (long i = 0; i < n; ++i)
        for (long j = 0; j < p; ++j, ++ell)
            for (long k = 0; k <= deg(b[i][j]); ++k)
                tmp_mat[k][ell] = b[i][j][k];
    Mat<zz_p> valB;
    mul(valB, vB, tmp_mat);

    Mat<zz_p> valBp(INIT_SIZE, n, p);
    Mat<zz_p> valC(INIT_SIZE, nb_points, m * p);
    Mat<zz_p> valCp;
    for (long i = 0; i < nb_points; ++i)
    {
        ell = 0;
        for (long u = 0; u < n; ++u)
            for (long v = 0; v < p; ++v, ++ell)
                valBp[u][v] = valB[i][ell];

        mul(valCp, valAp[i], valBp);

        ell = 0;
        for (long u = 0; u < m; ++u)
            for (long v = 0; v < p; ++v, ++ell)
                valC[i][ell] = valCp[u][v];
    }
    mul(tmp_mat, iV, valC);

    c.SetDims(m, p);
    ell = 0;
    for (long u = 0; u < m; ++u)
        for (long v = 0; v < p; ++v, ++ell)
        {
            c[u][v].SetLength(nb_points);
            for (long i = 0; i < nb_points; ++i)
                c[u][v][i] = tmp_mat[i][ell];
            c[u][v].normalize();
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
    mat_lzz_pX_lmultiplier(a, dB), primes(a.NumCols(), deg(a), dB) // TODO improve: this computes twice deg(a)...
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
    long dA = deg(a);

    if (t == TYPE_FFT_PRIME)
    {
        if (a.NumRows() <= 100)
            return std::unique_ptr<mat_lzz_pX_lmultiplier>(new mat_lzz_pX_lmultiplier_FFT_direct(a, dB));
        else
            return std::unique_ptr<mat_lzz_pX_lmultiplier>(new mat_lzz_pX_lmultiplier_FFT_matmul(a, dB));
    }
    if (t == TYPE_LARGE_PRIME)
    {
        if ((a.NumRows() <= 20 && dA <= 30) || (a.NumRows() <= 40 && dA <= 90) || dA <= 150)
            return std::unique_ptr<mat_lzz_pX_lmultiplier>(new mat_lzz_pX_lmultiplier_dense(a, dB));
        if (a.NumRows() <= 200)
            return std::unique_ptr<mat_lzz_pX_lmultiplier>(new mat_lzz_pX_lmultiplier_3_primes(a, dB));
        else
            return std::unique_ptr<mat_lzz_pX_lmultiplier>(new mat_lzz_pX_lmultiplier_geometric(a, dB));
    }
    else
    {
        if ((a.NumRows() <= 20 && dA <= 30) || (a.NumRows() <= 40 && dA <= 90) || dA <= 150)
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
