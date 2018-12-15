#include <NTL/matrix.h>
#include <NTL/mat_lzz_p.h>
#include <NTL/lzz_pX.h>
#include <NTL/BasicThreadPool.h>

#include "lzz_p_extra.h"
#include "mat_lzz_pX_extra.h"
#include "lzz_pX_CRT.h"

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
/* pairwise mul-add of two fftReps using long long's          */
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
/* c = a*b                                                    */
/* output may alias input; c does not have to be zero matrix  */
/* does not use Mat<zz_p> matrix multiplication               */
/*------------------------------------------------------------*/
void multiply_evaluate_FFT_direct(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b)
{
    if (&c == &a || &c == &b)
    {
        Mat<zz_pX> c2;
        multiply_evaluate_FFT_direct(c2, a, b);
        c = c2;
        return;
    }

#if defined(NTL_HAVE_LL_TYPE) && defined(NTL_HAVE_SP_LL_ROUTINES)
    Vec<Vec<fftRep>> valb;
    Vec<fftRep> vala;
    sp_reduce_struct red1;
    sp_ll_reduce_struct red2;
    long len, m, n, p, n0, dA, dB, K, pr, nb_slices, first_slice;
    Vec<ll_type> tmp;
    fftRep tmp_r;

    pr = zz_p::modulus();
    red1 = sp_PrepRem(pr);
    red2 = make_sp_ll_reduce_struct(pr);

    m = a.NumRows();
    n = a.NumCols();
    p = b.NumCols();

    dA = deg(a);
    dB = deg(b);
    K = NextPowerOfTwo(dA + dB + 1);

    valb.SetLength(p);
    for (long i = 0; i < p; i++)
    {
        valb[i].SetLength(n);
        for (long j = 0; j < n; j++)
            TofftRep(valb[i][j], b[j][i], K);
    }

    len = 1L << K;
    c.SetDims(m, p);
    vala.SetLength(n);

    tmp.SetLength(len);
    tmp_r = fftRep(INIT_SIZE, K);
    TofftRep(tmp_r, a[0][0], K);

    n0 = (1L << (2*(NTL_BITS_PER_LONG - NTL_SP_NBITS))) - 1;
    first_slice = n % n0;
    nb_slices = n / n0;
    if (first_slice == 0)
    {
        first_slice = n0;
        nb_slices--;
    }

    if (NumBits(pr) == NTL_SP_NBITS) // we can use normalized remainders; may be a bit faster
    {
        for (long i = 0; i < m; i++)
        {
            for (long j = 0; j < n; j++)
                TofftRep(vala[j], a[i][j], K);

            for (long k = 0; k < p; k++)
            {
                fftRep * vb = valb[k].elts();

                mul(tmp, vala[0], vb[0]);
                for (long j = 1; j < first_slice; j++)
                    mul_add(tmp, vala[j], vb[j]);

                long start = first_slice;
                for (long jj = 0; jj < nb_slices; jj++)
                {
                    for (long x = 0; x < len; x++)
                    {
                        tmp[x].lo = sp_ll_red_21_normalized(rem(tmp[x].hi, pr, red1), tmp[x].lo, pr, red2);
                        tmp[x].hi = 0;
                    }
                    for (long j = 0; j < n0; j++)
                        mul_add(tmp, vala[j+start], vb[j+start]);
                    start += n0;
                }

                long *tmp_ptr = &tmp_r.tbl[0][0];
                for (long x = 0; x < len; x++) 
                    tmp_ptr[x] = sp_ll_red_21_normalized(rem(tmp[x].hi, pr, red1), tmp[x].lo, pr, red2);
                FromfftRep(c[i][k], tmp_r, 0, dA + dB);
            }
        }
    }
    else
    {
        for (long i = 0; i < m; i++)
        {
            for (long j = 0; j < n; j++)
                TofftRep(vala[j], a[i][j], K);

            for (long k = 0; k < p; k++)
            {
                fftRep * vb = valb[k].elts();

                mul(tmp, vala[0], vb[0]);
                for (long j = 1; j < first_slice; j++)
                    mul_add(tmp, vala[j], vb[j]);

                long start = first_slice;
                for (long jj = 0; jj < nb_slices; jj++)
                {
                    for (long x = 0; x < len; x++)
                    {
                        tmp[x].lo = sp_ll_red_21(rem(tmp[x].hi, pr, red1), tmp[x].lo, pr, red2);
                        tmp[x].hi = 0;
                    }
                    for (long j = 0; j < n0; j++)
                        mul_add(tmp, vala[j+start], vb[j+start]);
                    start += n0;
                }

                long *tmp_ptr = &tmp_r.tbl[0][0];
                for (long x = 0; x < len; x++) 
                    tmp_ptr[x] = sp_ll_red_21(rem(tmp[x].hi, pr, red1), tmp[x].lo, pr, red2);
                FromfftRep(c[i][k], tmp_r, 0, dA + dB);
            }
        }
    }
#else
    // dimensions
    long m = a.NumRows();
    long n = a.NumCols();
    long p = b.NumCols();

    // degrees
    long d = deg(a) + deg(b);
    long K = NextPowerOfTwo(d + 1);

    // evaluate matrix b
    // valb[i][j] = evaluations of b[j][i]
    Vec<Vec<fftRep>> valb(INIT_SIZE, p);
    for (long i = 0; i < p; i++)
    {
        valb[i].SetLength(n);
        for (long j = 0; j < n; j++)
            TofftRep(valb[i][j], b[j][i], K);
    }

    c.SetDims(m, p);
    Vec<fftRep> vala(INIT_SIZE, n);
    fftRep tmp1 = fftRep(INIT_SIZE, K);
    fftRep tmp2 = fftRep(INIT_SIZE, K);
    // compute each row of the product c = a*b
    for (long i = 0; i < m; i++)
    {
        // vala[j] contains the evaluations of a[i][j]
        for (long j = 0; j < n; j++)
            TofftRep(vala[j], a[i][j], K);
        // compute c[i][k] = vala * valb[:][j]
        for (long k = 0; k < p; k++)
        {
            fftRep * vb = valb[k].elts();
            mul(tmp1, vala[0], vb[0]);
            for (long j = 1; j < n; j++)
            {
                mul(tmp2, vala[j], vb[j]);
                add(tmp1, tmp1, tmp2);
            }
            FromfftRep(c[i][k], tmp1, 0, d);
        }
    }
#endif
}

/*------------------------------------------------------------*/
/* c = a*b                                                    */
/* output may alias input; c does not have to be zero matrix  */
/* does not use Mat<zz_p> matrix multiplication               */
/*------------------------------------------------------------*/
void multiply_evaluate_FFT_direct_no_ll(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b)
{
    if (&c == &a || &c == &b)
    {
        Mat<zz_pX> c2;
        multiply_evaluate_FFT_direct(c2, a, b);
        c = c2;
        return;
    }

    // dimensions
    long m = a.NumRows();
    long n = a.NumCols();
    long p = b.NumCols();

    // degrees
    long d = deg(a) + deg(b);
    long K = NextPowerOfTwo(d + 1);

    // evaluate matrix b
    // valb[i][j] = evaluations of b[j][i]
    Vec<Vec<fftRep>> valb(INIT_SIZE, p);
    for (long i = 0; i < p; i++)
    {
        valb[i].SetLength(n);
        for (long j = 0; j < n; j++)
            TofftRep(valb[i][j], b[j][i], K);
    }

    c.SetDims(m, p);
    Vec<fftRep> vala(INIT_SIZE, n);
    fftRep tmp1 = fftRep(INIT_SIZE, K);
    fftRep tmp2 = fftRep(INIT_SIZE, K);
    // compute each row of the product c = a*b
    for (long i = 0; i < m; i++)
    {
        // vala[j] contains the evaluations of a[i][j]
        for (long j = 0; j < n; j++)
            TofftRep(vala[j], a[i][j], K);
        // compute c[i][k] = vala * valb[:][j]
        for (long k = 0; k < p; k++)
        {
            fftRep * vb = valb[k].elts();
            mul(tmp1, vala[0], vb[0]);
            for (long j = 1; j < n; j++)
            {
                mul(tmp2, vala[j], vb[j]);
                add(tmp1, tmp1, tmp2);
            }
            FromfftRep(c[i][k], tmp1, 0, d);
        }
    }
}

/*------------------------------------------------------------*/
/* c = a*b                                                    */
/* assumes FFT prime and p large enough                       */
/* output may alias input; c does not have to be zero matrix  */
/* uses Mat<zz_p> matrix multiplication                       */
/* --> for each matrix (a,b,c), list of all evaluations is    */
/* stored in a single array                                   */
/*------------------------------------------------------------*/
void multiply_evaluate_FFT_matmul1(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b)
{
    // dimensions
    const long s = a.NumRows();
    const long t = a.NumCols();
    const long u = b.NumCols();
    // degree of output
    const long d = deg(a) + deg(b);
    // points for FFT representation
    const long idxk = NextPowerOfTwo(d + 1);
    const long n = 1 << idxk;
    fftRep R(INIT_SIZE, idxk);

    // stores evaluations in a single vector for a, same for b
    const long st = s*t;
    Vec<zz_p> mat_valA(INIT_SIZE, n * st);
    const long tu = t*u;
    Vec<zz_p> mat_valB(INIT_SIZE, n * tu);

    // mat_valA[r*s*t + i*t + k] is a[i][k] evaluated at the r-th point
    for (long i = 0; i < s; ++i)
        for (long k = 0; k < t; ++k)
        {
            TofftRep(R, a[i][k], idxk);
            for (long r = 0, rst = 0; r < n; ++r, rst += st)
                mat_valA[rst + i*t+k].LoopHole() = R.tbl[0][r];
        }

    // mat_valB[r*s*t + i*t + k] is b[i][k] evaluated at the r-th point
    for (long i = 0; i < t; ++i)
        for (long k = 0; k < u; ++k)
        {
            TofftRep(R, b[i][k], idxk);
            for (long r = 0, rtu = 0; r < n; ++r, rtu += tu)
                mat_valB[rtu + i*u+k].LoopHole() = R.tbl[0][r];
        }

    // will store a evaluated at the r-th point, same for b
    Mat<zz_p> va(INIT_SIZE, s, t);
    Mat<zz_p> vb(INIT_SIZE, t, u);
    // will store the product evaluated at the r-th point, i.e. va*vb
    Mat<zz_p> vc(INIT_SIZE, s, u);

    // stores all evaluations of product c = a*b
    Vec<UniqueArray<long>> mat_valC(INIT_SIZE, s * u);
    for (long i = 0; i < s * u; ++i)
        mat_valC[i].SetLength(n);

    // compute pairwise products
    for (long j = 0, jst = 0, jtu = 0; j < n; ++j, jst += st, jtu += tu)
    {
        for (long i = 0; i < s; ++i)
            for (long k = 0; k < t; ++k)
                va[i][k].LoopHole() = mat_valA[jst + i*t + k]._zz_p__rep;
        for (long i = 0; i < t; ++i)
            for (long k = 0; k < u; ++k)
                vb[i][k].LoopHole() = mat_valB[jtu + i*u + k]._zz_p__rep;

        mul(vc, va, vb);

        for (long i = 0; i < s; ++i)
            for (long k = 0; k < u; ++k)
                mat_valC[i*u + k][j] = vc[i][k]._zz_p__rep;
    }

    // interpolate the evaluations stored in mat_valC back into c
    c.SetDims(s, u);

    for (long i = 0; i < s; ++i)
        for (long k = 0; k < u; ++k)
        {
            R.tbl[0].swap(mat_valC[i*u + k]);
            FromfftRep(c[i][k], R, 0, d);
        }
}

/*------------------------------------------------------------*/
/* c = a*b                                                    */
/* assumes FFT prime and p large enough                       */
/* output may alias input; c does not have to be zero matrix  */
/* uses Mat<zz_p> matrix multiplication                       */
/* --> for each matrix (a,b,c), evaluations are stored in an  */
/* array: e.g., for a, the entry i*t+j of the corresponding   */
/* array is itself an array containing all evaluations of the */
/* entry a[i][j], where t is the number of columns of a       */
/*------------------------------------------------------------*/
void multiply_evaluate_FFT_matmul2(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b)
{
    // dimensions
    long s = a.NumRows();
    long t = a.NumCols();
    long u = b.NumCols();
    // degree of output
    const long d = deg(a) + deg(b);
    // points for FFT representation
    long idxk = NextPowerOfTwo(d + 1);
    long n = 1 << idxk;
    fftRep R(INIT_SIZE, idxk);

    // matrix of evaluations of a: mat_valA[i][k][r] contains
    // the evaluation of a[i][k] at the r-th point
    Vec<UniqueArray<long>> mat_valA(INIT_SIZE, s*t);
    for (long i = 0; i < s; ++i)
        for (long k = 0; k < t; ++k)
        {
            R.tbl[0].SetLength(n);
            TofftRep(R, a[i][k], idxk);
            R.tbl[0].swap(mat_valA[i*t+k]);
        }

    // matrix of evaluations of b: mat_valB[i][k][r] contains
    // the evaluation of a[i][k] at the r-th point
    Vec<UniqueArray<long>> mat_valB(INIT_SIZE, t*u);
    for (long i = 0; i < t; ++i)
        for (long k = 0; k < u; ++k)
        {
            R.tbl[0].SetLength(n);
            TofftRep(R, b[i][k], idxk);
            R.tbl[0].swap(mat_valB[i*u+k]);
        }

    // vector containing the evaluations of the product c=a*b:
    // mat_valC[i*u+k][j] contains the evaluation of c[i][k]
    // at the j-th point
    Vec<UniqueArray<long>> mat_valC(INIT_SIZE, s*u);
    for (long i = 0; i < mat_valC.length(); ++i)
        mat_valC[i].SetLength(n);

    Mat<zz_p> va(INIT_SIZE, s, t);
    Mat<zz_p> vb(INIT_SIZE, t, u);
    Mat<zz_p> vc(INIT_SIZE, s, u);

    // for each point, compute the evaluation of c
    for (long j = 0; j < n; ++j)
    {
        for (long i = 0; i < s; ++i)
            for (long k = 0; k < t; ++k)
                va[i][k].LoopHole() = mat_valA[i*t+k][j];
        for (long i = 0; i < t; i++)
            for (long k = 0; k < u; k++)
                vb[i][k].LoopHole() = mat_valB[i*u+k][j];

        mul(vc, va, vb);

        for (long i = 0; i < s; ++i)
            for (long k = 0; k < u; ++k)
                mat_valC[i*u + k][j] = vc[i][k]._zz_p__rep;
    }

    // interpolate the evaluations mat_valC back into c
    c.SetDims(s, u);
    for (long i = 0; i < s; ++i)
        for (long k = 0; k < u; ++k)
        {
            R.tbl[0].swap(mat_valC[i*u + k]);
            FromfftRep(c[i][k], R, 0, d);
        }
}

/*------------------------------------------------------------*/
/* c = a*b                                                    */
/* assumes FFT prime and p large enough                       */
/* output may alias input; c does not have to be zero matrix  */
/* uses Mat<zz_p> matrix multiplication                       */
/* --> for each matrix (a,b), evaluations are stored in an    */
/* array of matrices: e.g., for a, the j-th entry of the      */
/* array contains the matrix of a evaluated at the j-th point */
/*------------------------------------------------------------*/
void multiply_evaluate_FFT_matmul3(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b)
{
    // dimensions
    long s = a.NumRows();
    long t = a.NumCols();
    long u = b.NumCols();
    // degree of output
    const long d = deg(a) + deg(b);
    // points for FFT representation
    long idxk = NextPowerOfTwo(d + 1);
    long n = 1 << idxk;
    fftRep R(INIT_SIZE, idxk);

    // matrix of evaluations of a: mat_valA[j] contains
    // the evaluation of a at the j-th point
    Vec<Mat<zz_p>> mat_valA(INIT_SIZE, n);
    for (long j = 0; j < n; ++j)
        mat_valA[j].SetDims(s,t);

    for (long i = 0; i < s; ++i)
        for (long k = 0; k < t; ++k)
        {
            TofftRep(R, a[i][k], idxk);
            for (long r = 0; r < n; ++r)
                mat_valA[r][i][k].LoopHole() = R.tbl[0][r];
        }

    // matrix of evaluations of b: mat_valB[j] contains
    // the evaluation of b at the j-th point
    Vec<Mat<zz_p>> mat_valB(INIT_SIZE, n);
    for (long j = 0; j < n; ++j)
        mat_valB[j].SetDims(t,u);
    for (long i = 0; i < t; ++i)
        for (long k = 0; k < u; ++k)
        {
            TofftRep(R, b[i][k], idxk);
            for (long r = 0; r < n; ++r)
                mat_valB[r][i][k].LoopHole() = R.tbl[0][r];
        }

    // compute pointwise products and store in mat_valA
    for (long j = 0; j < n; ++j)
        mul(mat_valA[j], mat_valA[j], mat_valB[j]);

    // interpolate c from the values in mat_valA
    c.SetDims(s, u);
    for (long i = 0; i < s; ++i)
        for (long k = 0; k < u; ++k)
        {
            for (long r = 0; r < n; ++r)
                R.tbl[0][r] = mat_valA[r][i][k]._zz_p__rep;
            FromfftRep(c[i][k], R, 0, d);
        }
}


// TODO: multi-threaded version of FFT_matmul: currently not integrated, needs more work, and other versions should
// be done as well
//void multiply_evaluate_FFT_matmul_threads(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b)
//{
//    zz_pContext context;
//    long s = a.NumRows();
//    long t = a.NumCols();
//    long u = b.NumCols();
//
//    long dA = deg(a);
//    long dB = deg(b);
//
//    long idxk = NextPowerOfTwo(dA + dB + 1);
//    fftRep R1(INIT_SIZE, idxk);
//    long n = 1 << idxk;
//
//
//    Vec<zz_p> mat_valA, mat_valB;
//    Vec<Vec<zz_p>> mat_valC;
//
//    mat_valA.SetLength(n * s * t);
//    mat_valB.SetLength(n * t * u);
//
//    long st = s*t;
//
//    context.save(); // to give the zz_p context to each thread
//    NTL_EXEC_RANGE(s,first,last)
//    context.restore(); // now all threads have the right zz_p context
//
//    fftRep R(INIT_SIZE, idxk);
//    for (long i = first; i < last; i++)
//    {
//        for (long k = 0; k < t; k++)
//        {
//            TofftRep(R, a[i][k], idxk);
//            long *frept = & R.tbl[0][0];
//            for (long r = 0, rst = 0; r < n; r++, rst += st)
//                mat_valA[rst + i*t + k] = frept[r];
//        }
//    }
//    NTL_EXEC_RANGE_END
//
//    long tu = t*u;
//    NTL_EXEC_RANGE(t, first, last)
//    context.restore();
//
//    fftRep R(INIT_SIZE, idxk);
//    for (long i = first; i < last; i++)
//    {
//        for (long k = 0; k < u; k++)
//        {
//            TofftRep(R, b[i][k], idxk);
//            long *frept = & R.tbl[0][0];
//            for (long r = 0, rtu = 0; r < n; r++, rtu += tu)
//                mat_valB[rtu + i*u + k] = frept[r];
//        }
//    }
//
//    R1 = R;
//    NTL_EXEC_RANGE_END
//
//    mat_valC.SetLength(s * u);
//    for (long i = 0; i < s * u; i++)
//        mat_valC[i].SetLength(n);
//
//    NTL_EXEC_RANGE(n, first, last)
//    context.restore();
//
//    Mat<zz_p> va, vb, vc; 
//    va.SetDims(s, t);
//    vb.SetDims(t, u);
//
//    for (long j = first, jst = st*first, jtu = tu*first; j < last; j++, jst += st, jtu += tu)
//    {
//        for (long i = 0; i < s; i++)
//        {
//            for (long k = 0; k < t; k++)
//            {
//                va[i][k] = mat_valA[jst + i*t + k];
//            }
//        }
//        for (long i = 0; i < t; i++)
//        {
//            for (long k = 0; k < u; k++)
//            {
//                vb[i][k] = mat_valB[jtu + i*u + k];
//            }
//        }
//
//        vc = va * vb;
//
//        for (long i = 0; i < s; i++)
//        {
//            for (long k = 0; k < u; k++)
//            {
//                mat_valC[i*u + k][j] = vc[i][k];
//            }
//        }
//    }
//    NTL_EXEC_RANGE_END
//
//    c.SetDims(s, u);
//
//    NTL_EXEC_RANGE(s, first, last)
//    context.restore();
//
//    Mat<zz_p> vc; 
//    fftRep R = R1;
//    for (long i = first; i < last; i++)
//    {
//        for (long k = 0; k < u; k++)
//        {
//            long *frept = & R.tbl[0][0];
//            for (long r = 0; r < n; r++)
//            {
//                frept[r] = rep(mat_valC[i*u + k][r]);
//            }
//            FromfftRep(c[i][k], R, 0, n - 1);
//        }
//    }
//    NTL_EXEC_RANGE_END
//}


/*------------------------------------------------------------*/
/* c = a*b                                                    */
/* output may alias input; c does not have to be zero matrix  */
/* chooses one of the two above                               */
/*------------------------------------------------------------*/
void multiply_evaluate_FFT(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b)
{
    long s = a.NumRows();
    long t = a.NumCols();
    long u = b.NumCols();

    long thresh;
    // TODO tune better
    if (NumBits(zz_p::modulus()) < 30)
        thresh = 20 * 20 * 20;
    else
        thresh = 45 * 45 * 45;

    if ((s * t * u) < thresh)  // fine for close-to-square matrices // FIXME tune better?
        multiply_evaluate_FFT_direct(c, a, b);
    else
        multiply_evaluate_FFT_matmul1(c, a, b);
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
