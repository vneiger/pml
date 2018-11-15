#include <NTL/lzz_pX.h>
#include <NTL/mat_lzz_p.h>
#include <limits.h>

#include "util.h"
#include "lzz_p_extra.h"
#include "lzz_pX_middle_product.h"
#include "structured_lzz_p.h"
#include "structured_lzz_pX.h"

NTL_CLIENT


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* square matrices such as                                    */
/* [a0  0  0]                                                 */
/* [a1 a0  0]                                                 */
/* [a2 a1 a0]                                                 */
/* (special case of toeplitz matrix)                          */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* sets dimensions to 0                                       */
/*------------------------------------------------------------*/
lower_triangular_toeplitz_lzz_pX::lower_triangular_toeplitz_lzz_pX()
{
    n = 0;
}

/*------------------------------------------------------------*/
/* input vector is as showed above                            */
/*------------------------------------------------------------*/
lower_triangular_toeplitz_lzz_pX::lower_triangular_toeplitz_lzz_pX(const Vec<zz_pX>& input)
{
    n = input.length();
    data = input;

    dX = -1; // x-deg
    dY = -1; // y-deg
    for (long i = 0; i < n; i++)
    {
        if (deg(input[i]) >= dX)
            dX = deg(input[i]);
        if (input[i] != 0)
            dY = i;
    }
}

/*------------------------------------------------------------*/
/* prepare for multiplication with rhs X-degree d             */
/*------------------------------------------------------------*/
void lower_triangular_toeplitz_lzz_pX::prepare_degree(long d)
{
    if (dY == -1 || d == -1)
    {
        return;
    }

    for (long i = 0; i < rhs_degX.length(); i++)
        if (d == rhs_degX[i])
            return;
    
    rhs_degX.append(d);

    zz_pX dataK;
    to_kronecker(dataK, data, dX + d);
    long K = NextPowerOfTwo( (dY + n)*(dX + d + 1) );
    fftRep fftK(INIT_SIZE, K);
    TofftRep(fftK, dataK, K);
    fft.append(fftK);
}


/*------------------------------------------------------------*/
/* getters                                                    */
/*------------------------------------------------------------*/
long lower_triangular_toeplitz_lzz_pX::NumCols() const
{
    return n;
}

long lower_triangular_toeplitz_lzz_pX::NumRows() const
{
    return n;
}

/*------------------------------------------------------------*/
/* turns M into a dense matrix                                */
/*------------------------------------------------------------*/
void lower_triangular_toeplitz_lzz_pX::to_dense(Mat<zz_pX>& Mdense) const 
{
    Mdense.SetDims(n, n);
    for (long j = 0; j < n; j++)
        for (long i = j; i < n; i++)
            Mdense[i][j] = data[i-j];
}

/*------------------------------------------------------------*/
/* right multiplication                                       */
/*------------------------------------------------------------*/
void lower_triangular_toeplitz_lzz_pX::mul_right(Vec<zz_pX>& res, const Vec<zz_pX>& input) const
{
    if (input.length() != n)
    {
        LogicError("Bad size for lower_triangular_toeplitz_lzz_pX right multiplication.");
    }

    if (&res == &input)
    {
        res = mul_right(input);
        return;
    }

    long d = -1; // x-deg of rhs
    long e = -1; // y-deg of rhs
    for (long i = 0; i < n; i++)
    {
        if (deg(input[i]) >= d)
            d = deg(input[i]);
        if (input[i] != 0)
            e = i;
    }
            
    if (d == -1 || dY == -1)
    {
        res.SetLength(n);
        for (long i = 0; i < n; i++)
            res[i] = 0;
        return;
    }

    long i0 = -1;
    for (long i = 0; i < rhs_degX.length(); i++)
        if ((rhs_degX[i] >= d && i0 == -1) || (rhs_degX[i] >= d && rhs_degX[i] < rhs_degX[i0]))
            i0 = i;
    
    zz_pX inK;
    long K;

    // found
    if (i0 >= 0)
    {
        long d_found = rhs_degX[i0];
        to_kronecker(inK, input, dX + d_found);
        K = NextPowerOfTwo( (dY + n) * (dX + d_found + 1) );

        fftRep fft_in(INIT_SIZE, K);
        TofftRep(fft_in, inK, K);
        mul(fft_in, fft_in, fft[i0]);
        zz_pX tmp;
        FromfftRep(tmp, fft_in, 0, n * (dX + d_found + 1) -1);
#ifdef __NTL_FIX_SIZE_2_FFT
        if (K == 1)
        {
            tmp = (1/to_zz_p(2)) * tmp;
        }
#endif   
        from_kronecker(res, tmp, dX + d_found);
    }
    else
    {
        zz_pX repK;
        to_kronecker(inK, input, dX + d);
        to_kronecker(repK, data, dX + d);
        K = NextPowerOfTwo( (dY + e + 1) * (dX + d + 1) );

        fftRep fft_in(INIT_SIZE, K), fft_this(INIT_SIZE, K);
        TofftRep(fft_in, inK, K);
        TofftRep(fft_this, repK, K);
        mul(fft_in, fft_in, fft_this);
        zz_pX tmp;
        long len = min(dY + e + 1, n);
        FromfftRep(tmp, fft_in, 0, len * (dX + d + 1) - 1);
#ifdef __NTL_FIX_SIZE_2_FFT
        if (K == 1)
        {
            tmp = (1/to_zz_p(2)) * tmp;
        }
#endif   
        from_kronecker(res, tmp, dX + d);
    }
    long ell = res.length();
    res.SetLength(n);
    for (long i = ell; i < n; i++)
        res[i] = 0;
}

/*------------------------------------------------------------*/
/* right multiplication mod x^s                               */
/*------------------------------------------------------------*/
void lower_triangular_toeplitz_lzz_pX::mul_right_trunc(Vec<zz_pX>& res, const Vec<zz_pX>& input, long s) const
{
    if (input.length() != n)
    {
        LogicError("Bad size for lower_triangular_toeplitz_lzz_pX right multiplication.");
    }

    if (&res == &input)
    {
        res = mul_right_trunc(input, s);
        return;
    }

    Vec<zz_pX> data_trunc;
    trunc(data_trunc, data, s);
    Vec<zz_pX> input_trunc;
    trunc(input_trunc, input, s);
    long dAx = deg(data_trunc);
    long d = deg(input_trunc);

    if (dAx == -1 || d == -1)
    {
        res.SetLength(n);
        for (long i = 0; i < n; i++)
            res[i] = 0;
        return;
    }
    zz_pX kro_in, kro_A, kro_res;
    to_kronecker(kro_in, input_trunc, dAx + d);
    to_kronecker(kro_A, data_trunc, dAx + d);
    kro_res = MulTrunc(kro_in, kro_A, n * (dAx + d + 1));
    from_kronecker(res, kro_res, dAx + d);
    long ell = res.length();
    res.SetLength(n);
    for (long i = 0; i < ell; i++)
        trunc(res[i], res[i], s);
    for (long i = ell; i < n; i++)
        res[i] = 0;
}

/*------------------------------------------------------------*/
/* left multiplication                                        */
/*------------------------------------------------------------*/
void lower_triangular_toeplitz_lzz_pX::mul_left(Vec<zz_pX>& res, const Vec<zz_pX>& input) const
{
    if (input.length() != n)
    {
        LogicError("Bad size for lower_triangular_toeplitz_lzz_pX left multiplication.");
    }

    if (&res == &input)
    {
        res = mul_left(input);
        return;
    }

    Vec<zz_pX> in_rev, out_rev;
    in_rev.SetLength(n);
    for (long i = 0; i < n; i++)
        in_rev[i] = input[n - 1 - i];
    mul_right(out_rev, in_rev);
    res.SetLength(n);
    for (long i = 0; i < n; i++)
        res[i] = out_rev[n - 1 - i];

}

/*------------------------------------------------------------*/
/* left multiplication mod x^s                                */
/*------------------------------------------------------------*/
void lower_triangular_toeplitz_lzz_pX::mul_left_trunc(Vec<zz_pX>& res, const Vec<zz_pX>& input, long s) const
{
    if (input.length() != n)
    {
        LogicError("Bad size for lower_triangular_toeplitz_lzz_pX right multiplication.");
    }

    if (&res == &input)
    {
        res = mul_left_trunc(input, s);
        return;
    }

    Vec<zz_pX> in_rev, out_rev;
    in_rev.SetLength(n);
    for (long i = 0; i < n; i++)
        in_rev[i] = input[n - 1 - i];
    mul_right_trunc(out_rev, in_rev, s);
    res.SetLength(n);
    for (long i = 0; i < n; i++)
        res[i] = out_rev[n - 1 - i];
}


// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
