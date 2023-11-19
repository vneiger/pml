#include <NTL/FFT_impl.h>
#include "mat_lzz_pX_multiply.h"

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
    const long len = min(x.len, y.len);
    for (long j = 0; j < len; ++j)
        ll_mul(z[j], xp[j], yp[j]);
}

/*------------------------------------------------------------*/
/* pairwise mul-add of two fftReps using long long's          */
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
/* c = a*b                                                    */
/* output may alias input; c does not have to be zero matrix  */
/* does not use Mat<zz_p> matrix multiplication               */
/*------------------------------------------------------------*/
void multiply_evaluate_FFT_direct_ll_type(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b)
{
#if defined(NTL_HAVE_LL_TYPE) && defined(NTL_HAVE_SP_LL_ROUTINES) \
     && defined(__GNUC__) && (__GNUC__ >= 4) && !defined(__INTEL_COMPILER)  && !defined(__clang__) \
     && defined (__x86_64__) && NTL_BITS_PER_LONG == 64
    if (&c == &a || &c == &b)
    {
        Mat<zz_pX> c2;
        multiply_evaluate_FFT_direct_ll_type(c2, a, b);
        c.swap(c2);
        return;
    }

    // for computations mod pr
    const long pr = zz_p::modulus();
    const sp_reduce_struct red1 = sp_PrepRem(pr);
    const sp_ll_reduce_struct red2 = make_sp_ll_reduce_struct(pr);

    // dimensions
    const long m = a.NumRows();
    const long n = a.NumCols();
    const long p = b.NumCols();

    // actual length of output (degree + 1)
    const long len_actual = deg(a) + deg(b) + 1;
    // length (number of points) for FFT representation
    const long K = NextPowerOfTwo(len_actual);
    const long len = FFTRoundUp(len_actual, K);

    const long n0 = (1L << (2*(NTL_BITS_PER_LONG - NTL_SP_NBITS))) - 1;
    const long first_slice = (n % n0 != 0) ? (n%n0) : n0;
    const long nb_slices = (n % n0 != 0) ? (n / n0) : (n/n0 - 1);

    Mat<fftRep> valb(INIT_SIZE, p, n);
    for (long i = 0; i < p; i++)
        for (long j = 0; j < n; j++)
            TofftRep_trunc(valb[i][j], b[j][i], K, len);

    c.SetDims(m, p);
    Vec<fftRep> vala(INIT_SIZE, n);

    fftRep tmp_r(INIT_SIZE, K);
    TofftRep_trunc(tmp_r, a[0][0], K, len);

    Vec<ll_type> tmp(INIT_SIZE, len);
    if (NumBits(pr) == NTL_SP_NBITS) // we can use normalized remainders; may be a bit faster
    {
        for (long i = 0; i < m; i++)
        {
            for (long j = 0; j < n; j++)
                TofftRep_trunc(vala[j], a[i][j], K, len);

            for (long k = 0; k < p; k++)
            {
                fftRep * vb = valb[k].elts();

                mul(tmp, vala[0], vb[0]);
                for (long j = 1; j < first_slice; ++j)
                    mul_add(tmp, vala[j], vb[j]);

                long start = first_slice;
                for (long jj = 0; jj < nb_slices; ++jj)
                {
                    for (long x = 0; x < len; ++x)
                    {
                        tmp[x].lo = sp_ll_red_21_normalized(rem(tmp[x].hi, pr, red1), tmp[x].lo, pr, red2);
                        tmp[x].hi = 0;
                    }
                    for (long j = 0; j < n0; ++j)
                        mul_add(tmp, vala[j+start], vb[j+start]);
                    start += n0;
                }

                long *tmp_ptr = &tmp_r.tbl[0][0];
                for (long x = 0; x < len; ++x)
                    tmp_ptr[x] = sp_ll_red_21_normalized(rem(tmp[x].hi, pr, red1), tmp[x].lo, pr, red2);
                FromfftRep(c[i][k], tmp_r, 0, len_actual-1);
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
                fftRep * vb = valb[k].elts();

                mul(tmp, vala[0], vb[0]);
                for (long j = 1; j < first_slice; ++j)
                    mul_add(tmp, vala[j], vb[j]);

                long start = first_slice;
                for (long jj = 0; jj < nb_slices; ++jj)
                {
                    for (long x = 0; x < len; ++x)
                    {
                        tmp[x].lo = sp_ll_red_21(rem(tmp[x].hi, pr, red1), tmp[x].lo, pr, red2);
                        tmp[x].hi = 0;
                    }
                    for (long j = 0; j < n0; ++j)
                        mul_add(tmp, vala[j+start], vb[j+start]);
                    start += n0;
                }

                long *tmp_ptr = &tmp_r.tbl[0][0];
                for (long x = 0; x < len; ++x)
                    tmp_ptr[x] = sp_ll_red_21(rem(tmp[x].hi, pr, red1), tmp[x].lo, pr, red2);
                FromfftRep(c[i][k], tmp_r, 0, len_actual-1);
            }
        }
    }
#else
    multiply_evaluate_FFT_direct(c,a,b);
#endif
}

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
        c.swap(c2);
        return;
    }

    // dimensions
    const long m = a.NumRows();
    const long n = a.NumCols();
    const long p = b.NumCols();

    // actual length of output (degree + 1)
    const long len_actual = deg(a) + deg(b) + 1;
    // actual length (number of points) for FFT
    const long K = NextPowerOfTwo(len_actual);
    const long len = FFTRoundUp(len_actual, K);

    // evaluate matrix b
    // valb[i][j] = evaluations of b[j][i]
    Mat<fftRep> valb(INIT_SIZE, p, n);
    for (long i = 0; i < p; ++i)
        for (long j = 0; j < n; ++j)
            TofftRep_trunc(valb[i][j], b[j][i], K, len);

    c.SetDims(m, p);
    Vec<fftRep> vala(INIT_SIZE, n);
    fftRep tmp1(INIT_SIZE, K);
    fftRep tmp2(INIT_SIZE, K);
    // compute each row of the product c = a*b
    for (long i = 0; i < m; ++i)
    {
        // vala[j] contains the evaluations of a[i][j]
        for (long j = 0; j < n; ++j)
            TofftRep_trunc(vala[j], a[i][j], K, len);
        // compute c[i][k] = vala * valb[:][j]
        for (long k = 0; k < p; ++k)
        {
            mul(tmp1, vala[0], valb[k][0]);
            for (long j = 1; j < n; ++j)
            {
                mul(tmp2, vala[j], valb[k][j]);
                add(tmp1, tmp1, tmp2);
            }
            FromfftRep(c[i][k], tmp1, 0, len_actual-1);
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
    // len_actual = actual length of output (degree+1)
    const long len_actual = deg(a) + deg(b) + 1;
    // len = number of points for truncated FFT
    const long idxk = NextPowerOfTwo(len_actual);
    const long len = FFTRoundUp(len_actual, idxk);
    const long n = 1<<idxk;
    fftRep R(INIT_SIZE, idxk);
    Vec<UniqueArray<long>> RR(INIT_SIZE, CACHE_LINE_SIZE);

    // stores evaluations in a single vector for a, same for b
    const long st = s*t;
    const long tu = t*u;
    Vec<long> mat_valA(INIT_SIZE, len * st);
    Vec<long> mat_valB(INIT_SIZE, len * tu);

    // mat_valA[r*s*t + i*t + k] is a[i][k] evaluated at the r-th point
    for (long i = 0; i < s; ++i)
    {
        for (long k = 0; k < t; k+=CACHE_LINE_SIZE)
        {
            const long kk_bnd = std::min((long)CACHE_LINE_SIZE, t-k);
            for (long kk = 0; kk < kk_bnd; ++kk)
            {
                R.tbl[0].SetLength(n);
                TofftRep_trunc(R, a[i][k+kk], idxk, len);
                R.tbl[0].swap(RR[kk]);
            }
            for (long r = 0, rst = 0; r < len; ++r, rst += st)
                for (long kk = 0; kk < kk_bnd; ++kk)
                    mat_valA[rst + i*t+k+kk] = RR[kk][r];
        }
    }

    // mat_valB[r*t*u + i*t + k] is b[i][k] evaluated at the r-th point
    for (long i = 0; i < t; ++i)
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
                    mat_valB[rtu + i*u+k+kk] = RR[kk][r];
        }

    // will store a evaluated at the r-th point, same for b
    Mat<zz_p> va(INIT_SIZE, s, t);
    Mat<zz_p> vb(INIT_SIZE, t, u);
    // will store the product evaluated at the r-th point, i.e. va*vb
    Vec<Mat<zz_p>> vc(INIT_SIZE, CACHE_LINE_SIZE);

    // stores all evaluations of product c = a*b
    Vec<UniqueArray<long>> mat_valC(INIT_SIZE, s * u);
    for (long i = 0; i < s * u; ++i)
        mat_valC[i].SetLength(n);

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
                    vb[i][k].LoopHole() = mat_valB[jtu + i*u + k];
            mul(vc[jj], va, vb);
        }

        for (long i = 0; i < s; ++i)
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
            R.tbl[0].swap(mat_valC[i*u + k]);
            FromfftRep(c[i][k], R, 0, len_actual-1);
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
void multiply_evaluate_FFT_matmul2(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b)
{
    // dimensions
    const long s = a.NumRows();
    const long t = a.NumCols();
    const long u = b.NumCols();
    // degree of output
    const long len_actual = deg(a) + deg(b) + 1;
    // points for FFT representation
    const long idxk = NextPowerOfTwo(len_actual);
    const long len = FFTRoundUp(len_actual, idxk);
    const long n = 1<<idxk;
    fftRep R(INIT_SIZE, idxk);

    // matrix of evaluations of a: mat_valA[i*t+k][r] contains
    // the evaluation of a[i][k] at the r-th point
    Vec<UniqueArray<long>> mat_valA(INIT_SIZE, s*t);
    for (long i = 0; i < s; ++i)
        for (long k = 0; k < t; ++k)
        {
            R.tbl[0].SetLength(n);
            TofftRep_trunc(R, a[i][k], idxk, len);
            R.tbl[0].swap(mat_valA[i*t+k]);
        }

    // matrix of evaluations of b: mat_valB[i*u+k][r] contains
    // the evaluation of a[i][k] at the r-th point
    Vec<UniqueArray<long>> mat_valB(INIT_SIZE, t*u);
    for (long i = 0; i < t; ++i)
        for (long k = 0; k < u; ++k)
        {
            R.tbl[0].SetLength(n);
            TofftRep_trunc(R, b[i][k], idxk, len);
            R.tbl[0].swap(mat_valB[i*u+k]);
        }

    // vector containing the evaluations of the product c=a*b:
    // mat_valC[i*u+k][j] contains the evaluation of c[i][k]
    // at the j-th point
    Vec<UniqueArray<long>> mat_valC(INIT_SIZE, s*u);
    for (long i = 0; i < mat_valC.length(); ++i)
        mat_valC[i].SetLength(n);

    Vec<Mat<zz_p>> va(INIT_SIZE, CACHE_LINE_SIZE);
    for (long jj = 0; jj < CACHE_LINE_SIZE; ++jj)
        va[jj].SetDims(s,t);
    Vec<Mat<zz_p>> vb(INIT_SIZE, CACHE_LINE_SIZE);
    for (long jj = 0; jj < CACHE_LINE_SIZE; ++jj)
        vb[jj].SetDims(t,u);
    Vec<Mat<zz_p>> vc(INIT_SIZE, CACHE_LINE_SIZE);

    // for each point, compute the evaluation of c
    for (long j = 0; j < len; j+=CACHE_LINE_SIZE)
    {
        const long jj_bnd = std::min((long)CACHE_LINE_SIZE, len-j);
        for (long i = 0; i < s; ++i)
            for (long k = 0; k < t; k+=MATRIX_BLOCK_SIZE)
            {
                const long kk_bnd = std::min((long)MATRIX_BLOCK_SIZE, t-k);
                for (long jj=0; jj<jj_bnd; ++jj)
                    for (long kk = 0; kk < kk_bnd; ++kk)
                        va[jj][i][k+kk].LoopHole() = mat_valA[i*t+k+kk][j+jj];
            }
        for (long i = 0; i < t; ++i)
            for (long k = 0; k < u; k+=MATRIX_BLOCK_SIZE)
            {
                const long kk_bnd = std::min((long)MATRIX_BLOCK_SIZE, u-k);
                for (long jj=0; jj<jj_bnd; ++jj)
                    for (long kk = 0; kk < kk_bnd; ++kk)
                        vb[jj][i][k+kk].LoopHole() = mat_valB[i*u+k+kk][j+jj];
            }

        for (long jj=0; jj<jj_bnd; ++jj)
            mul(vc[jj], va[jj], vb[jj]);

        for (long i = 0; i < s; ++i)
            for (long k = 0; k < u; k+=MATRIX_BLOCK_SIZE)
            {
                const long kk_bnd = std::min((long)MATRIX_BLOCK_SIZE, u-k);
                for (long kk = 0; kk < kk_bnd; ++kk)
                    for (long jj=0; jj<jj_bnd; ++jj)
                        mat_valC[i*u + k+kk][j+jj] = vc[jj][i][k+kk]._zz_p__rep;
            }
    }

    // interpolate the evaluations mat_valC back into c
    c.SetDims(s, u);
    for (long i = 0; i < s; ++i)
        for (long k = 0; k < u; ++k)
        {
            R.tbl[0].swap(mat_valC[i*u + k]);
            FromfftRep(c[i][k], R, 0, len_actual-1);
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
    const long s = a.NumRows();
    const long t = a.NumCols();
    const long u = b.NumCols();
    // actual length of output (degree +1)
    const long len_actual = deg(a) + deg(b) + 1;
    // number of points used in FFT representation
    const long idxk = NextPowerOfTwo(len_actual);
    const long len = FFTRoundUp(len_actual, idxk);
    const long n = 1 << idxk;
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
                R.tbl[0].SetLength(n);
                TofftRep_trunc(R, a[i][k+kk], idxk, len);
                RR[kk].swap(R.tbl[0]);
            }
            for (long r = 0; r < len; ++r)
                for (long kk=0; kk<kk_bnd; ++kk)
                    mat_valA[r][i][k+kk].LoopHole() = RR[kk][r];
        }

    // matrix of evaluations of b: mat_valB[j] contains
    // the evaluation of b at the j-th point
    Vec<Mat<zz_p>> mat_valB(INIT_SIZE, len);
    for (long j = 0; j < len; ++j)
        mat_valB[j].SetDims(t,u);

    for (long i = 0; i < t; ++i)
        for (long k = 0; k < u; k+=MATRIX_BLOCK_SIZE)
        {
            const long kk_bnd = std::min((long)MATRIX_BLOCK_SIZE, u-k);
            for (long kk=0; kk<kk_bnd; ++kk)
            {
                R.tbl[0].SetLength(n);
                TofftRep_trunc(R, b[i][k+kk], idxk, len);
                RR[kk].swap(R.tbl[0]);
            }
            for (long r = 0; r < len; ++r)
                for (long kk=0; kk<kk_bnd; ++kk)
                    mat_valB[r][i][k+kk].LoopHole() = RR[kk][r];
        }

    // compute pointwise products and store in mat_valA
    Mat<zz_p> tmp;
    for (long j = 0; j < len; ++j)
    {
        mul(tmp, mat_valA[j], mat_valB[j]);
        tmp.swap(mat_valA[j]);
    }

    // interpolate c from the values in mat_valA
    c.SetDims(s, u);
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
                FromfftRep(c[i][k+kk], R, 0, len_actual-1);
            }
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
    const long cube_dim = a.NumRows() * a.NumCols() * b.NumCols();
    const long d = (deg(a)+deg(b))/2;

    // TODO refine for degrees 1...50 (around 50)
    // --> evaluate dense is sometimes nice, then use it

    // could do automatic tuning?
    // (these thresholds should be reasonable on most recent machines;
    // note they were mostly designed by using close-to-square matrices, with
    // deg(a) and deg(b) similar)
    if (NumBits(zz_p::modulus()) < SMALL_PRIME_SIZE)
    {
        if (cube_dim < 4*4*4 || (cube_dim <= 6*6*6 && d > 100))
            multiply_evaluate_FFT_direct(c, a, b);
        else if (cube_dim <= 8*8*8)
            multiply_evaluate_FFT_direct_ll_type(c, a, b);
        else if (d<32)
            multiply_evaluate_dense(c, a, b);
        else if (cube_dim <= 16*16*16)
        {
            if (d>50)
                multiply_evaluate_FFT_matmul2(c, a, b);
            else // 32...50
                multiply_evaluate_dense2(c, a, b);
        }
        else if (d > 100 || (cube_dim < 64*64*64 && d > 80))
            multiply_evaluate_FFT_matmul1(c, a, b);
        else // dim = 16..63, d = 32..80  ||  dim = 65..., d=32..100
            multiply_evaluate_dense2(c, a, b);
    }
    else
    {
        // if < 4*4*4, or if <= 8*8*8 under some conditions, use direct
        if (cube_dim < 4*4*4 || (cube_dim == 4*4*4 && d>=90) || (cube_dim <= 8*8*8 && d >= 65))
            multiply_evaluate_FFT_direct(c, a, b);

        // if <= 8*8*8 and direct was not used,
        // or if < 20*20*20,
        // or if < 25*25*25 and d>=1300,
        // use direct_ll
        else if (cube_dim < 20*20*20 || (cube_dim < 25*25*25 && d>=1300))
            multiply_evaluate_FFT_direct_ll_type(c, a, b);

        // if between 20*20*20 and 25*25*25 and direct_ll was not used,
        // or if degree not too small,
        // use matmul1
        else if (cube_dim < 25*25*25 || d>512 ||
                 (cube_dim<48*48*48 && d>32) ||
                 (cube_dim<128*128*128 && d>256))
            multiply_evaluate_FFT_matmul1(c, a, b);

        // else, not small dimension and small degree: use matmul3
        else
            multiply_evaluate_FFT_matmul3(c, a, b);
    }
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
