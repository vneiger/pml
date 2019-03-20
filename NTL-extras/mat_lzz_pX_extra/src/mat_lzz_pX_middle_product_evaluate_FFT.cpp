#include <NTL/FFT_impl.h>
#include "mat_lzz_pX_multiply.h"

// FIXME work in progress:
// constants used to improve cache efficiency
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
    const long len = min(x.len, y.len);
    for (long j = 0; j < len; ++j)
        ll_mul(z[j], xp[j], yp[j]);
}

/*------------------------------------------------------------*/
/* pairwise product of two fftReps using long long's          */
/*------------------------------------------------------------*/
static inline void mul_add(Vec<ll_type>& z, const fftRep& x, const fftRep& y)
{
    const long *xp = &x.tbl[0][0];
    const long *yp = &y.tbl[0][0];
    const long len = min(x.len, y.len);
    for (long j = 0; j < len; ++j)
        ll_mul_add(z[j], xp[j], yp[j]);
}

#endif


/*------------------------------------------------------------*/
/* returns trunc( trunc(a, dA+1)*c div x^dA, dB+1 )           */
/* assumes FFT prime and p large enough                       */
/* output may alias input; b does not have to be zero matrix  */
/* does not use Mat<zz_p> matrix multiplication               */
/*------------------------------------------------------------*/
void middle_product_evaluate_FFT_direct_ll_type(Mat<zz_pX> & b, const Mat<zz_pX> & a, const Mat<zz_pX> & c, long dA, long dB)
{
#if defined(NTL_HAVE_LL_TYPE) && defined(NTL_HAVE_SP_LL_ROUTINES) \
     && defined(__GNUC__) && (__GNUC__ >= 4) && !defined(__INTEL_COMPILER)  && !defined(__clang__) \
     && defined (__x86_64__) && NTL_BITS_PER_LONG == 64
    if (&b == &a || &b == &c)
    {
        Mat<zz_pX> b2;
        middle_product_evaluate_FFT_direct_ll_type(b2, a, c, dA, dB);
        b.swap(b2);
        return;
    }

    // for computations mod pr
    const long pr = zz_p::modulus();
    const sp_reduce_struct red1 = sp_PrepRem(pr);
    const sp_ll_reduce_struct red2 = make_sp_ll_reduce_struct(pr);

    // dimensions
    const long m = a.NumRows();
    const long n = a.NumCols();
    const long p = c.NumCols();

    // sum of degree, length of FFT
    const long K = NextPowerOfTwo(dA+dB+1);
    const long len = (1<<K);

    const long n0 = (1L << (2*(NTL_BITS_PER_LONG - NTL_SP_NBITS))) - 1;
    const long first_slice = (n % n0 != 0) ? (n%n0) : n0;
    const long nb_slices = (n % n0 != 0) ? (n / n0) : (n/n0 - 1);

    Mat<fftRep> valc(INIT_SIZE, p, n);
    for (long i = 0; i < p; ++i)
        for (long j = 0; j < n; j++)
            TofftRep_trunc(valc[i][j], c[j][i], K, len);

    b.SetDims(m, p);
    Vec<fftRep> vala(INIT_SIZE, n);

    fftRep tmp_r(INIT_SIZE, K);
    TofftRep_trunc(tmp_r, a[0][0], K, len);

    Vec<ll_type> tmp(INIT_SIZE, len);
    if (NumBits(pr) == NTL_SP_NBITS) // we can use normalized remainders; may be a bit faster
    {
        for (long i = 0; i < m; ++i)
        {
            for (long j = 0; j < n; ++j)
                TofftRep_trunc(vala[j], a[i][j], K, len);

            for (long k = 0; k < p; ++k)
            {
                fftRep * vc = valc[k].elts();

                mul(tmp, vala[0], vc[0]);
                for (long j = 1; j < first_slice; ++j)
                    mul_add(tmp, vala[j], vc[j]);

                long start = first_slice;
                for (long jj = 0; jj < nb_slices; ++jj)
                {
                    for (long x = 0; x < len; ++x)
                    {
                        tmp[x].lo = sp_ll_red_21_normalized(rem(tmp[x].hi, pr, red1), tmp[x].lo, pr, red2);
                        tmp[x].hi = 0;
                    }
                    for (long j = 0; j < n0; ++j)
                        mul_add(tmp, vala[j+start], vc[j+start]);
                    start += n0;
                }

                long *tmp_ptr = &tmp_r.tbl[0][0];
                for (long x = 0; x < len; ++x)
                    tmp_ptr[x] = sp_ll_red_21_normalized(rem(tmp[x].hi, pr, red1), tmp[x].lo, pr, red2);
                FromfftRep(b[i][k], tmp_r, dA, dA+dB);
            }
        }
    }
    else
    {
        for (long i = 0; i < m; ++i)
        {
            for (long j = 0; j < n; ++j)
                TofftRep_trunc(vala[j], a[i][j], K, len);

            for (long k = 0; k < p; ++k)
            {
                fftRep * vc = valc[k].elts();

                mul(tmp, vala[0], vc[0]);
                for (long j = 1; j < first_slice; ++j)
                    mul_add(tmp, vala[j], vc[j]);

                long start = first_slice;
                for (long jj = 0; jj < nb_slices; ++jj)
                {
                    for (long x = 0; x < len; ++x)
                    {
                        tmp[x].lo = sp_ll_red_21(rem(tmp[x].hi, pr, red1), tmp[x].lo, pr, red2);
                        tmp[x].hi = 0;
                    }
                    for (long j = 0; j < n0; ++j)
                        mul_add(tmp, vala[j+start], vc[j+start]);
                    start += n0;
                }

                long *tmp_ptr = &tmp_r.tbl[0][0];
                for (long x = 0; x < len; ++x)
                    tmp_ptr[x] = sp_ll_red_21(rem(tmp[x].hi, pr, red1), tmp[x].lo, pr, red2);
                FromfftRep(b[i][k], tmp_r, dA, dA+dB);
            }
        }
    }
#else
    middle_product_evaluate_FFT_direct(b, a, c, dA, dB);
#endif
}

/*------------------------------------------------------------*/
/* returns trunc( trunc(a, dA+1)*c div x^dA, dB+1 )           */
/* assumes FFT prime and p large enough                       */
/* output may alias input; b does not have to be zero matrix  */
/* does not use Mat<zz_p> matrix multiplication               */
/*------------------------------------------------------------*/
void middle_product_evaluate_FFT_direct(Mat<zz_pX> & b, const Mat<zz_pX> & a, const Mat<zz_pX> & c, long dA, long dB)
{
    if (&b == &a || &b == &c)
    {
        Mat<zz_pX> b2;
        middle_product_evaluate_FFT_direct(b2, a, c, dA, dB);
        b.swap(b2);
        return;
    }

    // dimensions
    const long m = a.NumRows();
    const long n = a.NumCols();
    const long p = c.NumCols();
    // number of points for FFT
    const long K = NextPowerOfTwo(dA + dB + 1);
    const long len = (1<<K);

    Mat<fftRep> valc(INIT_SIZE, p, n);
    for (long i = 0; i < p; ++i)
        for (long j = 0; j < n; ++j)
            TofftRep_trunc(valc[i][j], c[j][i], K, len);

    b.SetDims(m, p);
    Vec<fftRep> vala(INIT_SIZE, n);
    fftRep tmp1(INIT_SIZE, K);
    fftRep tmp2(INIT_SIZE, K);
    for (long i = 0; i < m; ++i)
    {
        for (long j = 0; j < n; ++j)
            TofftRep_trunc(vala[j], a[i][j], K, len);
        for (long k = 0; k < p; ++k)
        {
            mul(tmp1, vala[0], valc[k][0]);
            for (long j = 1; j < n; ++j)
            {
                mul(tmp2, vala[j], valc[k][j]);
                add(tmp1, tmp1, tmp2);
            }
            FromfftRep(b[i][k], tmp1, dA, dA + dB);
        }
    }
}


/*------------------------------------------------------------*/
/* returns trunc( trunc(a, dA+1)*c div x^dA, dB+1 )           */
/* assumes FFT prime and p large enough                       */
/* output may alias input; b does not have to be zero matrix  */
/* uses Mat<zz_p> matrix multiplication                       */
/*------------------------------------------------------------*/
void middle_product_evaluate_FFT_matmul(Mat<zz_pX> & b, const Mat<zz_pX> & a, const Mat<zz_pX> & c, long dA, long dB)
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

void middle_product_evaluate_FFT_matmul1(Mat<zz_pX> & b, const Mat<zz_pX> & a, const Mat<zz_pX> & c, long dA, long dB)
{
    // dimensions
    const long s = a.NumRows();
    const long t = a.NumCols();
    const long u = c.NumCols();

    // number of points for FFT
    const long idxk = NextPowerOfTwo(dA + dB + 1);
    const long len = 1 << idxk;

    fftRep R(INIT_SIZE, idxk);
    Vec<UniqueArray<long>> RR(INIT_SIZE, CACHE_LINE_SIZE);

    // stores evaluations in a single vector for a, same for b
    const long st = s*t;
    const long tu = t*u;
    Vec<long> mat_valA(INIT_SIZE, len * st);
    Vec<long> mat_valC(INIT_SIZE, len * tu);

    // evals of a
    for (long i = 0; i < s; ++i)
    {
        for (long k = 0; k < t; k+=CACHE_LINE_SIZE)
        {
            const long kk_bnd = std::min((long)CACHE_LINE_SIZE, t-k);
            for (long kk = 0; kk < kk_bnd; ++kk)
            {
                R.tbl[0].SetLength(len);
                TofftRep(R, a[i][k+kk], idxk);
                R.tbl[0].swap(RR[kk]);
            }
            for (long r = 0, rst = 0; r < len; ++r, rst += st)
                for (long kk=0; kk < kk_bnd; ++kk)
                    mat_valA[rst + i*t + k+kk] = RR[kk][r];
        }
    }

    // evals of c
    for (long i = 0; i < t; ++i)
    {
        for (long k = 0; k < u; k+=CACHE_LINE_SIZE)
        {
            const long kk_bnd = std::min((long)CACHE_LINE_SIZE, u-k);
            for (long kk = 0; kk < kk_bnd; ++kk)
            {
                R.tbl[0].SetLength(len);
                TofftRep(R, c[i][k+kk], idxk);
                R.tbl[0].swap(RR[kk]);
            }
            for (long r = 0, rtu = 0; r < len; r++, rtu += tu)
                for (long kk = 0; kk < kk_bnd; ++kk)
                    mat_valC[rtu + i*u + k+kk] = RR[kk][r];
        }
    }

    // will store a evaluated at the r-th point, same for c
    Mat<zz_p> va(INIT_SIZE, s, t);
    Mat<zz_p> vc(INIT_SIZE, t, u);
    // store va*vc
    Vec<Mat<zz_p>> vb(INIT_SIZE, CACHE_LINE_SIZE);

    // stores the evaluations for b
    Vec<UniqueArray<long>> mat_valB(INIT_SIZE, s * u);
    for (long i = 0; i < s * u; ++i)
        mat_valB[i].SetLength(len);

    // compute pairwise products
    for (long j = 0, jst = 0, jtu = 0; j < len; j+=CACHE_LINE_SIZE)
    {
        const long jj_bnd = std::min((long)CACHE_LINE_SIZE, len-j);
        for (long jj = 0; jj < jj_bnd; ++jj, jst += st, jtu += tu)
        {
            for (long i = 0; i < s; ++i)
                for (long k = 0; k < t; ++k)
                    va[i][k].LoopHole() = mat_valA[jst + i*t + k];
            for (long i = 0; i < t; ++i)
                for (long k = 0; k < u; ++k)
                    vc[i][k].LoopHole() = mat_valC[jtu + i*u + k];
            mul(vb[jj], va, vc);
        }

        for (long i = 0; i < s; ++i)
            for (long k = 0; k < u; k+=MATRIX_BLOCK_SIZE)
            {
                const long kk_bnd = std::min((long)MATRIX_BLOCK_SIZE, u-k);
                for (long jj = 0; jj < jj_bnd; ++jj)
                    for (long kk = 0; kk < kk_bnd; ++kk)
                        mat_valB[i*u + k+kk][j+jj] = vb[jj][i][k+kk]._zz_p__rep;
            }
    }

    // interpolate
    b.SetDims(s, u);
    for (long i = 0; i < s; ++i)
        for (long k = 0; k < u; ++k)
        {
            R.tbl[0].swap(mat_valB[i*u+k]);
            FromfftRep(b[i][k], R, dA, dA + dB);
        }
}

void middle_product_evaluate_FFT_matmul2(Mat<zz_pX> & b, const Mat<zz_pX> & a, const Mat<zz_pX> & c, long dA, long dB)
{
    // dimensions
    const long s = a.NumRows();
    const long t = a.NumCols();
    const long u = c.NumCols();

    // initialize FFT representation
    const long idxk = NextPowerOfTwo(dA + dB + 1);
    const long len = 1 << idxk;

    fftRep R(INIT_SIZE, idxk);

    Vec<UniqueArray<long>> mat_valA(INIT_SIZE, s*t);
    for (long i = 0; i < s; ++i)
        for (long k = 0; k < t; ++k)
        {
            R.tbl[0].SetLength(len);
            TofftRep(R, a[i][k], idxk);
            R.tbl[0].swap(mat_valA[i*t + k]);
        }

    Vec<UniqueArray<long>> mat_valC(INIT_SIZE, t*u);
    for (long i = 0; i < t; ++i)
        for (long k = 0; k < u; ++k)
        {
            R.tbl[0].SetLength(len);
            TofftRep(R, c[i][k], idxk);
            R.tbl[0].swap(mat_valC[i*u + k]);
        }

    // evaluated matrices for a and c
    Vec<Mat<zz_p>> va(INIT_SIZE, CACHE_LINE_SIZE);
    for (long jj = 0; jj < CACHE_LINE_SIZE; ++jj)
         va[jj].SetDims(s, t);
    Vec<Mat<zz_p>> vc(INIT_SIZE, CACHE_LINE_SIZE);
    for (long jj = 0; jj < CACHE_LINE_SIZE; ++jj)
        vc[jj].SetDims(t, u);

    // will store vb[jj] = va[jj] * vc[jj]
    Vec<Mat<zz_p>> vb(INIT_SIZE, CACHE_LINE_SIZE);

    // vector containing evaluations of result
    Vec<UniqueArray<long>> mat_valB(INIT_SIZE, s * u);
    for (long i = 0; i < s * u; ++i)
        mat_valB[i].SetLength(len);

    // pointwise products
    for (long j = 0; j < len; j+=CACHE_LINE_SIZE)
    {
        const long jj_bnd = std::min((long)CACHE_LINE_SIZE, len-j);
        for (long i = 0; i < s; ++i)
            for (long k = 0; k < t; k+=MATRIX_BLOCK_SIZE)
            {
                const long kk_bnd = std::min((long)MATRIX_BLOCK_SIZE, t-k);
                for (long jj=0; jj<jj_bnd; ++jj)
                    for (long kk = 0; kk < kk_bnd; ++kk)
                        va[jj][i][k+kk].LoopHole() = mat_valA[i*t + k+kk][j+jj];
            }
        for (long i = 0; i < t; i++)
            for (long k = 0; k < u; k+=MATRIX_BLOCK_SIZE)
            {
                const long kk_bnd = std::min((long)MATRIX_BLOCK_SIZE, u-k);
                for (long jj=0; jj<jj_bnd; ++jj)
                    for (long kk = 0; kk < kk_bnd; ++kk)
                        vc[jj][i][k+kk].LoopHole() = mat_valC[i*u + k+kk][j+jj];
            }

        for (long jj=0; jj<jj_bnd; ++jj)
            mul(vb[jj], va[jj], vc[jj]);

        for (long i = 0; i < s; ++i)
            for (long k = 0; k < u; k+=MATRIX_BLOCK_SIZE)
            {
                const long kk_bnd = std::min((long)MATRIX_BLOCK_SIZE, u-k);
                for (long kk = 0; kk < kk_bnd; ++kk)
                    for (long jj=0; jj<jj_bnd; ++jj)
                        mat_valB[i*u + k+kk][j+jj] = vb[jj][i][k+kk]._zz_p__rep;
            }
    }

    // interpolate
    b.SetDims(s, u);
    for (long i = 0; i < s; ++i)
        for (long k = 0; k < u; ++k)
        {
            R.tbl[0].swap(mat_valB[i*u + k]);
            FromfftRep(b[i][k], R, dA, dA + dB);
        }
}

void middle_product_evaluate_FFT_matmul3(Mat<zz_pX> & b, const Mat<zz_pX> & a, const Mat<zz_pX> & c, long dA, long dB)
{
    // dimensions
    const long s = a.NumRows();
    const long t = a.NumCols();
    const long u = c.NumCols();

    // number of points in FFT representation
    const long idxk = NextPowerOfTwo(dA + dB + 1);
    const long len = 1 << idxk;

    fftRep R(INIT_SIZE, idxk);
    Vec<UniqueArray<long>> RR(INIT_SIZE, MATRIX_BLOCK_SIZE);

    // matrix of evaluations of a: mat_valA[j] contains
    // the evaluation of a at the j-th point
    Vec<Mat<zz_p>> mat_valA(INIT_SIZE, len);
    for (long j = 0; j < len; ++j)
        mat_valA[j].SetDims(s,t);

    for (long i = 0; i < s; ++i)
        for (long k = 0; k < t; k+=MATRIX_BLOCK_SIZE)
        {
            const long kk_bnd = std::min((long)MATRIX_BLOCK_SIZE, t-k);
            for (long kk=0; kk<kk_bnd; ++kk)
            {
                R.tbl[0].SetLength(len);
                TofftRep(R, a[i][k+kk], idxk);
                RR[kk].swap(R.tbl[0]);
            }
            for (long r = 0; r < len; ++r)
                for (long kk=0; kk<kk_bnd; ++kk)
                    mat_valA[r][i][k+kk].LoopHole() = RR[kk][r];
        }

    // matrix of evaluations of c, similar
    Vec<Mat<zz_p>> mat_valC(INIT_SIZE, len);
    for (long j = 0; j < len; ++j)
        mat_valC[j].SetDims(t,u);

    for (long i = 0; i < t; ++i)
        for (long k = 0; k < u; k+=MATRIX_BLOCK_SIZE)
        {
            const long kk_bnd = std::min((long)MATRIX_BLOCK_SIZE, u-k);
            for (long kk=0; kk<kk_bnd; ++kk)
            {
                R.tbl[0].SetLength(len);
                TofftRep(R, c[i][k+kk], idxk);
                RR[kk].swap(R.tbl[0]);
            }
            for (long r = 0; r < len; ++r)
                for (long kk=0; kk<kk_bnd; ++kk)
                    mat_valC[r][i][k+kk].LoopHole() = RR[kk][r];
        }

    Mat<zz_p> tmp;
    for (long j = 0; j < len; ++j)
    {
        mul(tmp, mat_valA[j], mat_valC[j]);
        tmp.swap(mat_valA[j]);
    }

    b.SetDims(s, u);
    for (long i = 0; i < s; ++i)
        for (long k = 0; k < u; k+=MATRIX_BLOCK_SIZE)
        {
            const long kk_bnd = std::min((long)MATRIX_BLOCK_SIZE, u-k);
            for (long r = 0; r < len; ++r)
                for (long kk=0; kk<kk_bnd; ++kk)
                    RR[kk][r] = mat_valA[r][i][k+kk]._zz_p__rep;
            for (long kk=0; kk<kk_bnd; ++kk)
            {
                R.tbl[0].swap(RR[kk]);
                FromfftRep(b[i][k+kk], R, dA, dA + dB);
            }
        }
}


/*------------------------------------------------------------*/
/* returns trunc( trunc(a, dA+1)*c div x^dA, dB+1 )           */
/* assumes FFT prime and p large enough                       */
/* output may alias input; b does not have to be zero matrix  */
/* uses Mat<zz_p> matrix multiplication                       */
/*------------------------------------------------------------*/
void middle_product_evaluate_FFT(Mat<zz_pX> & b, const Mat<zz_pX> & a, const Mat<zz_pX> & c, long dA, long dB)
{
    long s = a.NumRows();
    long t = a.NumCols();
    long u = c.NumCols();

    long thresh;
    if (NumBits(zz_p::modulus()) < 30)
        thresh = 20 * 20 * 20;
    else
        thresh = 45 * 45 * 45;

    if ((s * t * u) < thresh)  // fine for close-to-square matrices
        middle_product_evaluate_FFT_direct(b, a, c, dA, dB);
    else
        middle_product_evaluate_FFT_matmul(b, a, c, dA, dB);
}

void middle_product_evaluate_FFT_new(Mat<zz_pX> & b, const Mat<zz_pX> & a, const Mat<zz_pX> & c, long dA, long dB)
{
    //const long s = a.NumRows();
    //const long t = a.NumCols();
    //const long u = c.NumCols();
    const long cube_dim = a.NumRows() * a.NumCols() * c.NumCols();
    const long d = (dA+dB)/2;

    // TODO refine for degrees 1...50 (around 50)
    // --> evaluate dense is sometimes nice, then use it

    // could do automatic tuning?
    // (these thresholds should be reasonable on most recent machines;
    // note they were mostly designed by using close-to-square matrices, with
    // deg(a) and deg(b) similar)
    if (NumBits(zz_p::modulus()) < SMALL_PRIME_SIZE)
    {
        if (cube_dim < 4*4*4 || (cube_dim <= 6*6*6 && d > 100))
            middle_product_evaluate_FFT_direct(b, a, c, dA, dB);
        else if (cube_dim <= 8*8*8)
            middle_product_evaluate_FFT_direct_ll_type(b, a, c, dA, dB);
        else if (d<32)
            middle_product_evaluate_dense(b, a, c, dA, dB);
        else if (cube_dim <= 16*16*16)
        {
            if (d>50)
                middle_product_evaluate_FFT_matmul2(b, a, c, dA, dB);
            else // 32...50
                //middle_product_evaluate_dense2(b, a, c, dA, dB);
                middle_product_evaluate_dense(b, a, c, dA, dB);
        }
        else if (d > 100 || (cube_dim < 64*64*64 && d > 80))
            middle_product_evaluate_FFT_matmul1(b, a, c, dA, dB);
        else // dim = 16..63, d = 32..80  ||  dim = 65..., d=32..100
            //middle_product_evaluate_dense2(b, a, c, dA, dB);
            middle_product_evaluate_dense(b, a, c, dA, dB);
    }
    else
    {
        // if < 4*4*4, or if <= 8*8*8 under some conditions, use direct
        if (cube_dim < 4*4*4 || (cube_dim == 4*4*4 && d>=90) || (cube_dim <= 8*8*8 && d >= 65))
            middle_product_evaluate_FFT_direct(b, a, c, dA, dB);

        // if <= 8*8*8 and direct was not used,
        // or if < 20*20*20,
        // or if < 25*25*25 and d>=1300,
        // use direct_ll
        else if (cube_dim < 20*20*20 || (cube_dim < 25*25*25 && d>=1300))
            middle_product_evaluate_FFT_direct_ll_type(b, a, c, dA, dB);

        // if between 20*20*20 and 25*25*25 and direct_ll was not used,
        // or if degree not too small,
        // use matmul1
        else if (cube_dim < 25*25*25 || d>512 ||
                 (cube_dim<48*48*48 && d>32) ||
                 (cube_dim<128*128*128 && d>256))
            middle_product_evaluate_FFT_matmul1(b, a, c, dA, dB);

        // else, not small dimension and small degree: use matmul3
        else
            middle_product_evaluate_FFT_matmul3(b, a, c, dA, dB);
    }
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
