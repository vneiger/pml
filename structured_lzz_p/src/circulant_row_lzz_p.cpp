#include <NTL/lzz_pX.h>
#include <NTL/mat_lzz_p.h>

#include "lzz_p_extra.h"
#include "structured_lzz_p.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* matrices such as                                           */
/* [a3 a2 a1 a0]                                              */
/* [a0 a3 a2 a1]                                              */
/* [a1 a0 a3 a2]                                              */
/* (special case of toeplitz matrix)                          */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* sets dimensions to 0                                       */
/*------------------------------------------------------------*/
circulant_row_lzz_p::circulant_row_lzz_p()
{
    n = m = 0;
}

/*------------------------------------------------------------*/
/* input vector is as showed above                            */
/*------------------------------------------------------------*/
circulant_row_lzz_p::circulant_row_lzz_p(const Vec<zz_p>& input, long nrows)
{
    n = nrows;
    m = input.length();

    data = input;

    zz_pX poly1, poly2;

    SetCoeff(poly1, 0, input[m-1]);
    for (long i = 1; i < m; i++)
        SetCoeff(poly1, i, input[i-1]);

    long K = NextPowerOfTwo(m+m-1);
    fft1 = fftRep(INIT_SIZE, K);
    TofftRep(fft1, poly1, K);

    for (long i = 0; i < m; i++)
        SetCoeff(poly2, i, input[m - 1 - i]);

    fft3 = fftRep(INIT_SIZE, K);
    TofftRep(fft3, poly2, K);

    long leftover = n % m;
    if (leftover == 0)
        leftover = m;
    long L = NextPowerOfTwo(m + leftover - 1);
    fft2 = fftRep(INIT_SIZE, L);
    TofftRep(fft2, poly2, L);
}

/*------------------------------------------------------------*/
/* getters                                                    */
/*------------------------------------------------------------*/
long circulant_row_lzz_p::NumCols() const
{
    return m;
}

long circulant_row_lzz_p::NumRows() const
{
    return n;
}

/*------------------------------------------------------------*/
/* turns M into a dense matrix                                */
/*------------------------------------------------------------*/
void circulant_row_lzz_p::to_dense(Mat<zz_p>& Mdense) const 
{
    Mdense.SetDims(n, m);
    for (long i = 0; i < n; i++)
        for (long j = 0; j < m; j++)
        {
            long idx = i - j - 1;
            while (idx < 0)
                idx += m;
            idx = idx % m;
            Mdense[i][j] = data[idx];
        } 
}

/*------------------------------------------------------------*/
/* right multiplication                                       */
/*------------------------------------------------------------*/
void circulant_row_lzz_p::mul_right(Vec<zz_p>& res, const Vec<zz_p>& input) const
{
    if (input.length() != m)
    {
        LogicError("Wrong size for circulant_row_lzz_p right multiplication.");
    }

    if (&res == &input)
    {
        res = mul_right(input);
        return;
    }

    long K = NextPowerOfTwo(m + m - 1);
    fftRep fft_input = fftRep(INIT_SIZE, K);
    
    zz_pX input_X;
    input_X.rep.SetLength(m);
    zz_p *cf = input_X.rep.elts();
    for (long i = 0; i < m; i++)
            cf[i] = input[i];
    input_X.normalize();
    TofftRep(fft_input, input_X, K);
    mul(fft_input, fft1, fft_input);

    Vec<zz_p> res_all;
    res_all.SetLength(m + m);
    FromfftRep(res_all.elts(), fft_input, 0, 2*m - 2);
    res_all[2*m - 1] = 0;

#ifdef __NTL_FIX_SIZE_2_FFT
    if (K == 1)
    {
        res_all[0] /= 2;
        res_all[1] /= 2;
    }
#endif   
    
    res.SetLength(n);
    long todo = min(n, m);
    for (long i = 0; i < todo; i++)
        res[i] = res_all[i] + res_all[i + m];

    if (m < n)
    {
        long done = m;
        long left = n-m;
        do 
        {
            todo = min(left, m);
            for (long i = done; i < done + todo; i++)
                res[i] = res[i-m];
            done += todo;
            left -= todo;
        } 
        while (left >= 1);
    }
}

/*------------------------------------------------------------*/
/* right matrix multiplication                                */
/*------------------------------------------------------------*/
void circulant_row_lzz_p::mul_right(Mat<zz_p>& res, const Mat<zz_p>& input) const
{
    if (input.NumRows() != m)
    {
        LogicError("Wrong size for circulant_row_lzz_p right matrix multiplication.");
    }

    if (&res == &input)
    {
        res = mul_right(input);
        return;
    }
    Error ("not implemented yet."); // todo: move to structured_lzz_p.cpp
}

/*------------------------------------------------------------*/
/* left multiplication                                        */
/*------------------------------------------------------------*/
void circulant_row_lzz_p::mul_left(Vec<zz_p>& res, const Vec<zz_p>& input) const
{
    if (input.length() != n)
    {
        LogicError("Wrong size for circulant_row_lzz_p left multiplication.");
    }

    if (&res == &input)
    {
        res = mul_left(input);
        return;
    }

    res.SetLength(m);
    
    long leftover = n % m;
    long nb = n / m;   // number of size-m chunks

    if (leftover == 0)
    {
        leftover = m;
        nb--;
    }

    zz_pX input_X;
    Vec<zz_p> out;

    // last chunk
    out.SetLength(m + leftover);
    long L = NextPowerOfTwo(m + leftover - 1);
    fftRep fft_l = fftRep(INIT_SIZE, L);

    input_X.rep.SetLength(leftover);
    zz_p *cf = input_X.rep.elts();
    for (long i = 0; i < leftover; i++)
        cf[i] = input[i + nb * m];
    input_X.normalize();
    TofftRep(fft_l, input_X, L);
    mul(fft_l, fft2, fft_l);
    FromfftRep(out.elts(), fft_l, 0, m + leftover - 2);

#ifdef __NTL_FIX_SIZE_2_FFT
    if (L == 1)
    {
        out[0] /= 2;
        out[1] /= 2;
    }
#endif   


    out[m + leftover - 1] = 0;
    
    for (long i = 0; i < leftover; i++)
        res[i] = out[i] + out[i + m];
    for (long i = leftover; i < m; i++)
        res[i] = out[i];

    for (long r = 0; r < nb; r++)
    {
        out.SetLength(2 * m);
        long K = NextPowerOfTwo(2 * m - 1);
        fftRep fft_k = fftRep(INIT_SIZE, K);

        input_X.rep.SetLength(m);
        zz_p *cf = input_X.rep.elts();
        for (long i = 0; i < m; i++)
            cf[i] = input[i + r * m];
        input_X.normalize();
        TofftRep(fft_k, input_X, K);
        mul(fft_k, fft3, fft_k);
        FromfftRep(out.elts(), fft_k, 0, m + m - 2);

#ifdef __NTL_FIX_SIZE_2_FFT
    if (K == 1)
    {
        out[0] /= 2;
        out[1] /= 2;
    }
#endif   
        out[m + m - 1] = 0;
    
        for (long i = 0; i < m; i++)
            res[i] += (out[i] + out[i + m]);
    }

}

/*------------------------------------------------------------*/
/* left matrix multiplication                                 */
/*------------------------------------------------------------*/
void circulant_row_lzz_p::mul_left(Mat<zz_p>& res, const Mat<zz_p>& input) const
{
    if (input.NumCols() != n)
    {
        LogicError("Wrong size for circulant_row_lzz_p left matrix multiplication.");
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
