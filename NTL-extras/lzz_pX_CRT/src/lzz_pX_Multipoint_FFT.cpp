#include <NTL/lzz_pX.h>
#include <NTL/vector.h>

#include "util.h"
#include "lzz_pX_CRT.h"

NTL_CLIENT


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* multipoint evaluation over zz_p                            */
/* FFT points                                                 */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

// from NTL's old FFT
static long RevInc(long a, long k)
{
    long j, m;

    j = k; 
    m = 1L << (k-1);

    while (j && (m & a)) {
        a ^= m;
        m >>= 1;
        j--;
    }
    if (j) a ^= m;
    return a;
}

// from NTL's old FFT
static void BitReverse_init(Vec<long> & indices, long k)
{
    long n = (1L << k);
    indices.SetLength(n);
    long *rev = indices.elts();
    long i, j;
    for (i = 0, j = 0; i < n; i++, j = RevInc(j, k))
        rev[i] = j;
}

/*------------------------------------------------------------*/
/* constructor for FFT points                                 */
/* n has to be a power of 2, p must be FFT prime              */
/* checks whether the output is bit-reversed or not           */
/*------------------------------------------------------------*/
zz_pX_Multipoint_FFT::zz_pX_Multipoint_FFT(long n)
{
    k = NextPowerOfTwo(n);
    if (n != (1L << k))
    {
        LogicError("Wrong length for zz_pX_Multipoint_FFT");
    }

    this->n = n;
    if (zz_pInfo->p_info == NULL)
    {
        LogicError("Attempt to init a zz_pX_Multipoint_FFT without zz_p::FFTInit");
    }

    if (k > zz_pInfo->MaxRoot) 
        ResourceError("Too many points for FFT");

    BitReverse_init(indices, k);

    zz_pX X;
    X = 0;
    SetCoeff(X, 1, 1);
    wk = fftRep(INIT_SIZE, k);
    TofftRep(wk, X, k); // makes sure that wk has length n (needed for NTL >= 11)

    zz_p pow1 = to_zz_p(1);
    zz_p pow2 = to_zz_p(1);

    zz_p z1 = to_zz_p(wk.tbl[0][1]);
    zz_p z2 = to_zz_p(wk.tbl[0][n >> 1]);
    do_bit_reverse = 1;

    for (long i = 0; i < n; i++)
    {
        if (! ((wk.tbl[0][i] == pow1) || (wk.tbl[0][indices[i]] == pow2)))
            LogicError("Neither bit-reversed nor direct order in FFT");
        if (! (wk.tbl[0][indices[i]] == pow2))
            do_bit_reverse = 0;
        pow1 *= z1;
        pow2 *= z2;
    }
}

/*------------------------------------------------------------*/
/* does a forward FFT                                         */
/*------------------------------------------------------------*/
void zz_pX_Multipoint_FFT::evaluate(Vec<zz_p>& val, const zz_pX& f) const 
{
    fftRep frep(INIT_SIZE, k);
    TofftRep(frep, f, k);
    long *frept = &frep.tbl[0][0];

    val.SetLength(n);
    for (long i = 0; i < n; i++)
    {
        val[i] = frept[i];
    }
}

/*------------------------------------------------------------*/
/* does an inverse FFT                                        */
/*------------------------------------------------------------*/
void zz_pX_Multipoint_FFT::interpolate(zz_pX& f, const Vec<zz_p>& val) {
    long *frept = &wk.tbl[0][0];

    for (long i = 0; i < n; i++)
    {
        frept[i] = rep(val[i]);
    }

    FromfftRep(f, wk, 0, n-1);

#ifdef __NTL_FIX_SIZE_2_FFT
    if (n == 2) // for n=2, the normalization is different in version 11.1.0
    {
        f = (1 / to_zz_p(2)) * f;
    }
#endif
}

/*------------------------------------------------------------*/
/* transpose forward FFT                                      */
/*------------------------------------------------------------*/
void zz_pX_Multipoint_FFT::t_evaluate(zz_pX& f, const Vec<zz_p>& val, long output_size) const 
{
    f = 0;
    zz_pX tmp;
    tmp.rep.SetLength(n);
    zz_p * coeffs_tmp = tmp.rep.elts();

    if (output_size == -1)
        output_size = n;

    if (! do_bit_reverse)
    {
        for (long i = 0; i < n; i++)
        {
            coeffs_tmp[i] = val[i];
        }
    }
    else
    {
        for (long i = 0; i < n; i++)
        {
            coeffs_tmp[i] = val[indices[i]];
        }
    }
    tmp.normalize();

    fftRep frep(INIT_SIZE, k);
    TofftRep(frep, tmp, k);
    long *frept = &frep.tbl[0][0];

    f.rep.SetLength(output_size);
    zz_p * coeffs_f = f.rep.elts();

    if (! do_bit_reverse)
    {
        for (long i = 0; i < output_size; i++)
        {
            coeffs_f[i] = frept[i];
        }
    }
    else 
    {
        for (long i = 0; i < output_size; i++)
        {
            coeffs_f[i] = frept[indices[i]];
        }
    }
    f.normalize();
}

/*------------------------------------------------------------*/
/* transposed inverse FFT                                     */
/*------------------------------------------------------------*/
void zz_pX_Multipoint_FFT::t_interpolate(Vec<zz_p>& val, const zz_pX& f) {
    val.SetLength(n);

    if (n == 1)
    {
        val[0] = coeff(f, 0);
        return;
    }

    if (n == 2) 
    {
        zz_p half = 1/to_zz_p(2);
        val[0] = (coeff(f,0) + coeff(f,1)) * half;
        val[1] = (coeff(f,0) - coeff(f,1)) * half;
        return;
    }

    zz_pX tmp;
    long *frept = &wk.tbl[0][0];

    if (! do_bit_reverse)
    {
        for (long i = 0; i < n; i++)
        {
            frept[i] = rep(coeff(f, i));
        }
    }
    else
    {
        for (long i = 0; i < n; i++)
        {
            frept[i] = rep(coeff(f, indices[i]));
        }
    }

    FromfftRep(tmp, wk, 0, n-1);

#ifdef __NTL_FIX_SIZE_2_FFT
    if (n == 2) // for n=2, the normalization is different in version 11.1.0
    {
        tmp = (1 / to_zz_p(2)) * tmp;
    }
#endif

    if (! do_bit_reverse)
    {
        for (long i = 0; i < n; i++)
            val[i] = coeff(tmp, i);
    }
    else
    {
        for (long i = 0; i < n; i++)
            val[indices[i]] = coeff(tmp, i);
    }
}

/*------------------------------------------------------------*/
/* returns a zz_pX_Multipoint_FFT with at least n points      */
/*------------------------------------------------------------*/
zz_pX_Multipoint_FFT get_FFT_points(long n)
{
    long k = NextPowerOfTwo(n);
    return zz_pX_Multipoint_FFT(1L << k);
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
