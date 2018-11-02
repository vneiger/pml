#include <NTL/matrix.h>
#include <NTL/mat_lzz_p.h>
#include <NTL/lzz_pX.h>
#include <NTL/BasicThreadPool.h>

#include "util.h"
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
/* c = a*b                                                    */
/* output may alias input; c does not have to be zero matrix  */
/* does not use Mat<zz_p> matrix multiplication               */
/*------------------------------------------------------------*/
void multiply_evaluate_direct_FFT(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b)
{

    if (&c == &a || &c == &b)
    {
        Mat<zz_pX> c2;
        multiply_evaluate_direct_FFT(c2, a, b);
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
    Vec<Vec<fftRep>> valb;
    Vec<fftRep> vala;

    long m = a.NumRows();
    long n = a.NumCols();
    long p = b.NumCols();

    long dA = deg(a);
    long dB = deg(b);
    long K = NextPowerOfTwo(dA + dB + 1);

    valb.SetLength(p);
    for (long i = 0; i < p; i++)
    {
        valb[i].SetLength(n);
        for (long j = 0; j < n; j++)
            TofftRep(valb[i][j], b[j][i], K);
    }

    c.SetDims(m, p);
    vala.SetLength(n);
    fftRep tmp1 = fftRep(INIT_SIZE, K);
    fftRep tmp2 = fftRep(INIT_SIZE, K);
    for (long i = 0; i < m; i++)
    {
        for (long j = 0; j < n; j++)
            TofftRep(vala[j], a[i][j], K);
        for (long k = 0; k < p; k++)
        {
            fftRep * vb = valb[k].elts();
            mul(tmp1, vala[0], vb[0]);
            for (long j = 1; j < n; j++)
            {
                mul(tmp2, vala[j], vb[j]);
                add(tmp1, tmp1, tmp2);
            }
            FromfftRep(c[i][k], tmp1, 0, dA + dB);
        }
    }
#endif
}
