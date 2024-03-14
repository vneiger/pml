#include <NTL/lzz_pX.h>
#include <NTL/mat_lzz_p.h>

#include "lzz_p_extra.h"
#include "lzz_pX_middle_product.h"
#include "structured_lzz_p.h"

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
lower_triangular_toeplitz_lzz_p::lower_triangular_toeplitz_lzz_p()
{
    n = 0;
}

/*------------------------------------------------------------*/
/* input vector is as showed above                            */
/*------------------------------------------------------------*/
lower_triangular_toeplitz_lzz_p::lower_triangular_toeplitz_lzz_p(const Vec<zz_p>& input)
{
    n = input.length();
    data = input;
    for (long i = 0; i < n; i++)
    {
        SetCoeff(dataX, i, input[i]);
    }

    long K = NextPowerOfTwo(n+n-1);
    fft = fftRep(INIT_SIZE, K);
    TofftRep(fft, dataX, K);
}

/*------------------------------------------------------------*/
/* getters                                                    */
/*------------------------------------------------------------*/
long lower_triangular_toeplitz_lzz_p::NumCols() const
{
    return n;
}

long lower_triangular_toeplitz_lzz_p::NumRows() const
{
    return n;
}

/*------------------------------------------------------------*/
/* turns M into a dense matrix                                */
/*------------------------------------------------------------*/
void lower_triangular_toeplitz_lzz_p::to_dense(Mat<zz_p>& Mdense) const 
{
    Mdense.SetDims(n, n);
    for (long j = 0; j < n; j++)
        for (long i = j; i < n; i++)
            Mdense[i][j] = data[i-j];
}

/*------------------------------------------------------------*/
/* right multiplication                                       */
/*------------------------------------------------------------*/
void lower_triangular_toeplitz_lzz_p::mul_right(Vec<zz_p>& res, const Vec<zz_p>& input) const
{
    if (input.length() != n)
    {
        LogicError("Bad size for lower_triangular_toeplitz_lzz_p right multiplication.");
    }

    if (&res == &input)
    {
        res = mul_right(input);
        return;
    }

    res.SetLength(n);

    zz_pX input_X;
    input_X.rep.SetLength(n);
    zz_p *cf = input_X.rep.elts();
    for (long i = 0; i < n; i++)
        cf[i] = input[i];
    input_X.normalize();

    if (n <= NTL_zz_pX_MUL_CROSSOVER/2)
    {
        zz_pX tmp;
        MulTrunc(tmp, dataX, input_X, n);
        for (long i = 0; i < n; i++)
            res[i] = coeff(tmp, i);
    }
    else
    {
        long K = NextPowerOfTwo(n+n-1);
        fftRep fft_input = fftRep(INIT_SIZE, K);
        TofftRep(fft_input, input_X, K);
        mul(fft_input, fft_input, fft);
        FromfftRep(res.elts(), fft_input, 0, n-1);
    }
}


/*------------------------------------------------------------*/
/* right matrix multiplication                                */
/*------------------------------------------------------------*/
void lower_triangular_toeplitz_lzz_p::mul_right(Mat<zz_p>& res, const Mat<zz_p>& input) const
{
    if (input.NumRows() != n)
    {
        LogicError("Bad size for lower_triangular_toeplitz_lzz_p right matrix multiplication.");
    }

    if (&res == &input)
    {
        res = mul_right(input);
        return;
    }
}


/*------------------------------------------------------------*/
/* left multiplication                                        */
/*------------------------------------------------------------*/
void lower_triangular_toeplitz_lzz_p::mul_left(Vec<zz_p>& res, const Vec<zz_p>& input) const
{
    if (input.length() != n)
    {
        LogicError("Bad size for lower_triangular_toeplitz_lzz_p left multiplication.");
    }

    if (&res == &input)
    {
        res = mul_left(input);
        return;
    }


    res.SetLength(n);

    zz_pX input_X;
    input_X.rep.SetLength(n);
    zz_p *cf = input_X.rep.elts();
    for (long i = 0; i < n; i++)
        cf[i] = input[n - 1 - i];
    input_X.normalize();

    if (n <= NTL_zz_pX_MUL_CROSSOVER/2)
    {
        zz_pX tmp;
        MulTrunc(tmp, dataX, input_X, n);
        for (long i = 0; i < n; i++)
            res[i] = coeff(tmp, n - 1 - i);
    }
    else
    {
        Vec<zz_p> tmp;
        tmp.SetLength(n);
        long K = NextPowerOfTwo(n+n-1);
        fftRep fft_input = fftRep(INIT_SIZE, K);
        TofftRep(fft_input, input_X, K);
        mul(fft_input, fft_input, fft);
        FromfftRep(tmp.elts(), fft_input, 0, n-1);
        for (long i = 0 ; i < n; i++)
            res[i] = tmp[n - 1 - i];
    }

}

/*------------------------------------------------------------*/
/* left matrix multiplication                                 */
/*------------------------------------------------------------*/
void lower_triangular_toeplitz_lzz_p::mul_left(Mat<zz_p>& res, const Mat<zz_p>& input) const
{
    if (input.NumCols() != n)
    {
        LogicError("Bad size for lower_triangular_toeplitz_lzz_p left matrix multiplication.");
    }

    if (&res == &input)
    {
        res = mul_left(input);
        return;
    }
}


// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
