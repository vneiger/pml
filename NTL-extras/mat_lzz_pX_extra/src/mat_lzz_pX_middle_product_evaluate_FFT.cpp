#include <NTL/FFT_impl.h>
#include "mat_lzz_pX_multiply.h"

NTL_CLIENT

#if defined(NTL_HAVE_LL_TYPE) && defined(NTL_HAVE_SP_LL_ROUTINES)

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
#if defined(NTL_HAVE_LL_TYPE) && defined(NTL_HAVE_SP_LL_ROUTINES)
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

    // stores evaluations in a single vector for a, same for b
    const long st = s*t;
    const long tu = t*u;
    Vec<long> mat_valA(INIT_SIZE, len * st);
    Vec<long> mat_valC(INIT_SIZE, len * tu);

    // evals of a
    for (long i = 0; i < s; ++i)
    {
        for (long k = 0; k < t; ++k)
        {
            TofftRep(R, a[i][k], idxk);
            long *frept = & R.tbl[0][0];
            for (long r = 0, rst = 0; r < len; r++, rst += st)
                mat_valA[rst + i*t + k] = frept[r];
        }
    }

    // evals of c
    for (long i = 0; i < t; ++i)
    {
        for (long k = 0; k < u; ++k)
        {
            TofftRep(R, c[i][k], idxk);
            long *frept = & R.tbl[0][0];
            for (long r = 0, rtu = 0; r < len; r++, rtu += tu)
                mat_valC[rtu + i*u + k] = frept[r];
        }
    }

    // will store a evaluated at the r-th point, same for c
    Mat<zz_p> va(INIT_SIZE, s, t);
    Mat<zz_p> vc(INIT_SIZE, t, u);
    // store va*vc
    Mat<zz_p> vb;

    // stores the evaluations for b
    Vec<UniqueArray<long>> mat_valB(INIT_SIZE, s * u);
    for (long i = 0; i < s * u; ++i)
        mat_valB[i].SetLength(len);

    // compute pairwise products
    for (long j = 0, jst = 0, jtu = 0; j < len; ++j, jst += st, jtu += tu)
    {
        for (long i = 0; i < s; ++i)
            for (long k = 0; k < t; ++k)
                va[i][k].LoopHole() = mat_valA[jst + i*t + k];
        for (long i = 0; i < t; ++i)
            for (long k = 0; k < u; ++k)
                vc[i][k].LoopHole() = mat_valC[jtu + i*u + k];

        mul(vb, va, vc);

        for (long i = 0; i < s; ++i)
            for (long k = 0; k < u; ++k)
                mat_valB[i*u + k][j] = vb[i][k]._zz_p__rep;
    }

    b.SetDims(s, u);
    for (long i = 0; i < s; ++i)
        for (long k = 0; k < u; ++k)
        {
            long *frept = & R.tbl[0][0];
            for (long r = 0; r < len; ++r)
                frept[r] = mat_valB[i*u + k][r];
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
    Mat<zz_p> va(INIT_SIZE, s, t);
    Mat<zz_p> vc(INIT_SIZE, t, u);
    // vb = va * vc
    Mat<zz_p> vb;

    Vec<UniqueArray<long>> mat_valB(INIT_SIZE, s * u);
    for (long i = 0; i < s * u; ++i)
        mat_valB[i].SetLength(len);

    // pointwise products
    for (long j = 0; j < len; ++j)
    {
        for (long i = 0; i < s; ++i)
            for (long k = 0; k < t; ++k)
                va[i][k].LoopHole() = mat_valA[i*t + k][j];
        for (long i = 0; i < t; i++)
            for (long k = 0; k < u; ++k)
                vc[i][k].LoopHole() = mat_valC[i*u + k][j];

        mul(vb, va, vc);

        for (long i = 0; i < s; ++i)
            for (long k = 0; k < u; ++k)
                mat_valB[i*u + k][j] = vb[i][k]._zz_p__rep;
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

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
