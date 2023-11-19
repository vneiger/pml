#include <NTL/lzz_pX.h>
#include <NTL/mat_lzz_p.h>

#include "lzz_p_extra.h"
#include "structured_lzz_p.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* matrices such as                                           */
/* [a0 a2 a1 a0]                                              */
/* [a1 a0 a2 a1]                                              */
/* [a2 a1 a0 a2]                                              */
/* (special case of toeplitz matrix)                          */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* sets dimensions to 0                                       */
/*------------------------------------------------------------*/
circulant_column_lzz_p::circulant_column_lzz_p()
{
    n = m = 0;
}

/*------------------------------------------------------------*/
/* input vector is as showed above                            */
/*------------------------------------------------------------*/
circulant_column_lzz_p::circulant_column_lzz_p(const Vec<zz_p>& input, long ncols)
{
    n = input.length();
    m = ncols;

    data = input;

    zz_pX poly1, poly2;

    SetCoeff(poly1, 0, input[0]);
    for (long i = 1; i < n; i++)
        SetCoeff(poly1, i, input[n-1-(i-1)]);

    long K = NextPowerOfTwo(n+n-1);
    fft1 = fftRep(INIT_SIZE, K);
    TofftRep(fft1, poly1, K);

    for (long i = 0; i < n; i++)
        SetCoeff(poly2, i, input[i]);

    fft3 = fftRep(INIT_SIZE, K);
    TofftRep(fft3, poly2, K);

    long leftover = m % n;
    if (leftover == 0)
        leftover = n;
    long L = NextPowerOfTwo(n + leftover - 1);
    fft2 = fftRep(INIT_SIZE, L);
    TofftRep(fft2, poly2, L);
}

/*------------------------------------------------------------*/
/* getters                                                    */
/*------------------------------------------------------------*/
long circulant_column_lzz_p::NumCols() const
{
    return m;
}

long circulant_column_lzz_p::NumRows() const
{
    return n;
}

/*------------------------------------------------------------*/
/* turns M into a dense matrix                                */
/*------------------------------------------------------------*/
void circulant_column_lzz_p::to_dense(Mat<zz_p>& Mdense) const 
{
    Mdense.SetDims(n, m);
    for (long i = 0; i < m; i++)
        for (long j = 0; j < n; j++)
        {
            long idx = i - j - 1;
            while (idx < 0)
                idx += n;
            idx = idx % n;
            Mdense[j][i] = data[n-1-idx];
        } 
}

/*------------------------------------------------------------*/
/* right multiplication                                       */
/*------------------------------------------------------------*/
void circulant_column_lzz_p::mul_right(Vec<zz_p>& res, const Vec<zz_p>& input) const
{
    if (input.length() != m)
    {
        LogicError("Wrong size for circulant_column_lzz_p right multiplication.");
    }

    if (&res == &input)
    {
        res = mul_right(input);
        return;
    }

    res.SetLength(n);
    
    long leftover = m % n;
    long nb = m / n;   // number of size-n chunks

    if (leftover == 0)
    {
        leftover = n;
        nb--;
    }

    zz_pX input_X;
    Vec<zz_p> out;

    // last chunk
    out.SetLength(n + leftover);
    long L = NextPowerOfTwo(n + leftover - 1);
    fftRep fft_l = fftRep(INIT_SIZE, L);

    input_X.rep.SetLength(leftover);
    zz_p *cf = input_X.rep.elts();
    for (long i = 0; i < leftover; i++)
        cf[i] = input[i + nb * n];
    input_X.normalize();
    TofftRep(fft_l, input_X, L);
    mul(fft_l, fft2, fft_l);
    FromfftRep(out.elts(), fft_l, 0, n + leftover - 2);

#ifdef __NTL_FIX_SIZE_2_FFT
    if (L == 1)
    {
        out[0] /= 2;
        out[1] /= 2;
    }
#endif   

    out[n + leftover - 1] = 0;
    
    for (long i = 0; i < leftover; i++)
        res[i] = out[i] + out[i + n];
    for (long i = leftover; i < n; i++)
        res[i] = out[i];

    for (long r = 0; r < nb; r++)
    {
        out.SetLength(2 * n);
        long K = NextPowerOfTwo(2 * n - 1);
        fftRep fft_k = fftRep(INIT_SIZE, K);

        input_X.rep.SetLength(n);
        zz_p *cf = input_X.rep.elts();
        for (long i = 0; i < n; i++)
            cf[i] = input[i + r * n];
        input_X.normalize();
        TofftRep(fft_k, input_X, K);
        mul(fft_k, fft3, fft_k);
        FromfftRep(out.elts(), fft_k, 0, n + n - 2);

#ifdef __NTL_FIX_SIZE_2_FFT
    if (K == 1)
    {
        out[0] /= 2;
        out[1] /= 2;
    }
#endif   
        out[n + n - 1] = 0;
    
        for (long i = 0; i < n; i++)
            res[i] += (out[i] + out[i + n]);
    }
}

/*------------------------------------------------------------*/
/* right matrix multiplication                                */
/*------------------------------------------------------------*/
void circulant_column_lzz_p::mul_right(Mat<zz_p>& res, const Mat<zz_p>& input) const
{
    if (input.NumRows() != m)
    {
        LogicError("Wrong size for circulant_column_lzz_p right matrix multiplication.");
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
void circulant_column_lzz_p::mul_left(Vec<zz_p>& res, const Vec<zz_p>& input) const
{
    if (input.length() != n)
    {
        LogicError("Wrong size for circulant_column_lzz_p left multiplication.");
    }

    if (&res == &input)
    {
        res = mul_left(input);
        return;
    }

    long K = NextPowerOfTwo(n + n - 1);
    fftRep fft_input = fftRep(INIT_SIZE, K);
    
    zz_pX input_X;
    input_X.rep.SetLength(m);
    zz_p *cf = input_X.rep.elts();
    for (long i = 0; i < n; i++)
            cf[i] = input[i];
    input_X.normalize();
    TofftRep(fft_input, input_X, K);
    mul(fft_input, fft1, fft_input);

    Vec<zz_p> res_all;
    res_all.SetLength(n + n);
    FromfftRep(res_all.elts(), fft_input, 0, 2*n - 2);
    res_all[2*n - 1] = 0;

#ifdef __NTL_FIX_SIZE_2_FFT
    if (K == 1)
    {
        res_all[0] /= 2;
        res_all[1] /= 2;
    }
#endif   
    
    res.SetLength(m);
    long todo = min(m, n);
    for (long i = 0; i < todo; i++)
        res[i] = res_all[i] + res_all[i + n];

    if (n < m)
    {
        long done = n;
        long left = m-n;
        do 
        {
            todo = min(left, n);
            for (long i = done; i < done + todo; i++)
                res[i] = res[i-n];
            done += todo;
            left -= todo;
        } 
        while (left >= 1);
    }
}

/*------------------------------------------------------------*/
/* left matrix multiplication                                 */
/*------------------------------------------------------------*/
void circulant_column_lzz_p::mul_left(Mat<zz_p>& res, const Mat<zz_p>& input) const
{
    if (input.NumCols() != n)
    {
        LogicError("Wrong size for circulant_column_lzz_p left matrix multiplication.");
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
