#include <NTL/lzz_pX.h>
#include <NTL/mat_lzz_p.h>

#include "structured_lzz_pX.h"

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
circulant_row_lzz_pX::circulant_row_lzz_pX()
{
    n = m = 0;
}

/*------------------------------------------------------------*/
/* input vector is as showed above                            */
/*------------------------------------------------------------*/
circulant_row_lzz_pX::circulant_row_lzz_pX(const Vec<zz_pX>& input, long nrows)
{
    n = nrows;
    m = input.length();
    data = input;
    dataRev.SetLength(m);
    dataRev[0] = input[m - 1];
    for (long i = 1; i < m; i++)
        dataRev[i] = input[i-1];

    dX = -1; // x-deg
    dY = -1; // y-deg
    for (long i = 0; i < m; i++)
    {
        if (deg(dataRev[i]) >= dX)
            dX = deg(dataRev[i]);
        if (dataRev[i] != 0)
            dY = i;
    }
}

/*------------------------------------------------------------*/
/* getters                                                    */
/*------------------------------------------------------------*/
long circulant_row_lzz_pX::NumCols() const
{
    return m;
}

long circulant_row_lzz_pX::NumRows() const
{
    return n;
}

/*------------------------------------------------------------*/
/* turns M into a dense matrix                                */
/*------------------------------------------------------------*/
void circulant_row_lzz_pX::to_dense(Mat<zz_pX>& Mdense) const 
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
/* prepare for multiplication with rhs X-degree d             */
/*------------------------------------------------------------*/
void circulant_row_lzz_pX::prepare_degree(long d)
{
    if (dY == -1 || d == -1)
    {
        return;
    }

    for (long i = 0; i < rhs_degX.length(); i++)
        if (d == rhs_degX[i])
            return;
    
    rhs_degX.append(d);
    lhs_degX.append(d);

    zz_pX dataK;
    to_kronecker(dataK, dataRev, dX + d);
    long K = NextPowerOfTwo( (dY + m) * (dX + d + 1) );
    fftRep fftK(INIT_SIZE, K);
    TofftRep(fftK, dataK, K);
    fft.append(fftK);

    Vec<zz_pX> tmp;
    zz_pX mat_kro;
    tmp.SetLength(m);
    for (long i = 0; i < m; i++)
        tmp[i] = data[m - 1 - i];
    to_kronecker(mat_kro, tmp, dX + d);

    long leftover = n % m;
    long nb = n / m;   // number of size-m chunks
    if (leftover == 0)
    {
        leftover = m;
        nb--;
    }
    long L = NextPowerOfTwo( (m + leftover - 1) * (dX + d + 1));    
    fftRep fft2_here(INIT_SIZE, L);
    TofftRep(fft2_here, mat_kro, L);
    fft2.append(fft2_here);

    K = NextPowerOfTwo( (2 * m - 1) * (dX + d + 1) );
    fftRep fft3_here(INIT_SIZE, K);
    TofftRep(fft3_here, mat_kro, K);
    fft3.append(fft3_here);
}

/*------------------------------------------------------------*/
/* right multiplication                                       */
/*------------------------------------------------------------*/
void circulant_row_lzz_pX::mul_right(Vec<zz_pX>& res, const Vec<zz_pX>& input) const
{
    if (input.length() != m)
    {
        LogicError("Wrong size for circulant_row_lzz_pX right multiplication.");
    }

    if (&res == &input)
    {
        res = mul_right(input);
        return;
    }

    long d = -1; // x-deg of rhs
    long e = -1; // y-deg of rhs
    for (long i = 0; i < m; i++)
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

    Vec<zz_pX> coeff_res;

    // found
    if (i0 >= 0)
    {
        long d_found = rhs_degX[i0];
        long K = NextPowerOfTwo( (dY + m) * (dX + d_found + 1) );
        fftRep fft_input = fftRep(INIT_SIZE, K);
        zz_pX input_X;
        to_kronecker(input_X, input, dX + d_found);
        TofftRep(fft_input, input_X, K);
        mul(fft_input, fft[i0], fft_input);
        zz_pX res_kro;
        FromfftRep(res_kro, fft_input, 0, (dY + m) * (dX + d_found + 1) - 1 );
#ifdef __NTL_FIX_SIZE_2_FFT
        if (K == 1)
        {
            res_kro = (1/to_zz_p(2)) * res_kro;
        }
#endif   
        from_kronecker(coeff_res, res_kro, dX + d_found);
    }
    else
    {
        zz_pX input_X, repK;
        to_kronecker(input_X, input, dX + d);
        to_kronecker(repK, dataRev, dX + d);
        long K = NextPowerOfTwo( (dY + e + 1) * (dX + d + 1) );

        fftRep fft_in(INIT_SIZE, K), fft_this(INIT_SIZE, K);
        TofftRep(fft_in, input_X, K);
        TofftRep(fft_this, repK, K);
        mul(fft_in, fft_in, fft_this);
        zz_pX res_kro;

        FromfftRep(res_kro, fft_in, 0, (dY + e + 1) * (dX + d + 1) - 1);
#ifdef __NTL_FIX_SIZE_2_FFT
        if (K == 1)
        {
            res_kro = (1/to_zz_p(2)) * res_kro;
        }
#endif   
        from_kronecker(coeff_res, res_kro, dX + d);
    }

    long ell = coeff_res.length();
    coeff_res.SetLength(2 * m);
    for (long i = ell; i < 2 * m; i++)
        coeff_res[i] = 0;
    
    res.SetLength(n);
    long todo = min(n, m);
    for (long i = 0; i < todo; i++)
        res[i] = coeff_res[i] + coeff_res[i + m];

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
/* left multiplication                                        */
/*------------------------------------------------------------*/
void circulant_row_lzz_pX::mul_left(Vec<zz_pX>& res, const Vec<zz_pX>& input) const
{
    if (input.length() != n)
    {
        LogicError("Wrong size for circulant_row_lzz_pX left multiplication.");
    }

    if (&res == &input)
    {
        res = mul_left(input);
        return;
    }

    long d = -1; // x-deg of lhs
    for (long i = 0; i < n; i++)
    {
        if (deg(input[i]) >= d)
            d = deg(input[i]);
    }

    res.SetLength(m);

    if (d == -1 || dY == -1)
    {
        for (long i = 0; i < m; i++)
            res[i] = 0;
        return;
    }

    long leftover = n % m;
    long nb = n / m;   // number of size-m chunks
    if (leftover == 0)
    {
        leftover = m;
        nb--;
    }


    long i0 = -1;
    for (long i = 0; i < lhs_degX.length(); i++)
        if ((lhs_degX[i] >= d && i0 == -1) || (lhs_degX[i] >= d && lhs_degX[i] < lhs_degX[i0]))
            i0 = i;

    // found
    if (i0 >= 0)
    {
        Vec<zz_pX> tmp, out;
        zz_pX res_kro, input_X;
        long ell;

        // last chunk
        long d_found = lhs_degX[i0];
        long L = NextPowerOfTwo( (m + leftover - 1) * (dX + d_found + 1));
        fftRep fft_l = fftRep(INIT_SIZE, L);

        tmp.SetLength(leftover);
        for (long i = 0; i < leftover; i++)
            tmp[i] = input[i + nb*m];
        to_kronecker(input_X, tmp, dX + d_found);
        TofftRep(fft_l, input_X, L);

        mul(fft_l, fft2[i0], fft_l);
        FromfftRep(res_kro, fft_l, 0, (m + leftover - 1) * (dX + d_found + 1) - 1);
#ifdef __NTL_FIX_SIZE_2_FFT
        if (L == 1)
        {
            res_kro = (1/to_zz_p(2)) * res_kro;
        }
#endif   
        from_kronecker(out, res_kro, dX + d_found);
        ell = out.length();
        out.SetLength(m + leftover);
        for (long i = ell; i < out.length(); i++)
            out[i] = 0;
        
        for (long i = 0; i < leftover; i++)
            res[i] = out[i] + out[i + m];
        for (long i = leftover; i < m; i++)
            res[i] = out[i];

        tmp.SetLength(m);
        for (long r = 0; r < nb; r++)
        {
            long K = NextPowerOfTwo( (2 * m - 1) * (dX + d_found + 1) );
            fftRep fft_k = fftRep(INIT_SIZE, K);
            input_X = 0;

            for (long i = 0; i < m; i++)
                tmp[i] = input[i + r*m];
            to_kronecker(input_X, tmp, dX + d_found);
            TofftRep(fft_k, input_X, K);

            mul(fft_k, fft3[i0], fft_k);

            FromfftRep(res_kro, fft_k, 0, (m + m - 1) * (dX + d_found + 1) - 1);
#ifdef __NTL_FIX_SIZE_2_FFT
            if (K == 1)
            {
                res_kro = (1/to_zz_p(2)) * res_kro;
            }
#endif   
            from_kronecker(out, res_kro, dX + d_found);
            ell = out.length();
            out.SetLength(m + m);
            for (long i = ell; i < out.length(); i++)
                out[i] = 0;

            for (long i = 0; i < m; i++)
                res[i] += (out[i] + out[i + m]);
        }
    }
    else
    {
        Vec<zz_pX> out, tmp;
        zz_pX input_X, mat_kro, res_kro;

        // last chunk
        long L = NextPowerOfTwo( (m + leftover - 1) * (dX + d + 1));
        fftRep fft_l = fftRep(INIT_SIZE, L), fft2_here = fftRep(INIT_SIZE, L);

        tmp.SetLength(leftover);
        for (long i = 0; i < leftover; i++)
            tmp[i] = input[i + nb*m];
        to_kronecker(input_X, tmp, dX + d);
        TofftRep(fft_l, input_X, L);

        tmp.SetLength(m);
        for (long i = 0; i < m; i++)
            tmp[i] = data[m - 1 - i];
        to_kronecker(mat_kro, tmp, dX + d);
        TofftRep(fft2_here, mat_kro, L);

        mul(fft_l, fft2_here, fft_l);

        FromfftRep(res_kro, fft_l, 0, (m + leftover - 1) * (dX + d + 1) - 1);
#ifdef __NTL_FIX_SIZE_2_FFT
        if (L == 1)
        {
            res_kro = (1/to_zz_p(2)) * res_kro;
        }
#endif   
        from_kronecker(out, res_kro, dX + d);

        long ell = out.length();
        out.SetLength(m + leftover);
        for (long i = ell; i < out.length(); i++)
            out[i] = 0;
        
        for (long i = 0; i < leftover; i++)
            res[i] = out[i] + out[i + m];
        for (long i = leftover; i < m; i++)
            res[i] = out[i];

        for (long r = 0; r < nb; r++)
        {
            long K = NextPowerOfTwo( (2 * m - 1) * (dX + d + 1) );
            fftRep fft_k = fftRep(INIT_SIZE, K), fft3_here = fftRep(INIT_SIZE, K);

            input_X = 0;
            for (long i = 0; i < m; i++)
                tmp[i] = input[i + r*m];
            to_kronecker(input_X, tmp, dX + d);
            TofftRep(fft_k, input_X, K);
            
            mat_kro = 0;
            for (long i = 0; i < m; i++)
                tmp[i] = data[m - 1 - i];
            to_kronecker(mat_kro, tmp, dX + d);
            TofftRep(fft3_here, mat_kro, K);

            mul(fft_k, fft3_here, fft_k);
            FromfftRep(res_kro, fft_k, 0, (m + m - 1) * (dX + d + 1) - 1);
            
#ifdef __NTL_FIX_SIZE_2_FFT
            if (K == 1)
            {
                res_kro = (1/to_zz_p(2)) * res_kro;
            }
#endif   
            from_kronecker(out, res_kro, dX + d);
            ell = out.length();
            out.SetLength(m + m);
            for (long i = ell; i < out.length(); i++)
                out[i] = 0;

            for (long i = 0; i < m; i++)
                res[i] += (out[i] + out[i + m]);
        }

    }
}

/*------------------------------------------------------------*/
/* right multiplication mod x^s                               */
/*------------------------------------------------------------*/
void circulant_row_lzz_pX::mul_right_trunc(Vec<zz_pX>& out, const Vec<zz_pX>& in, long s) const
{
    if (in.length() != m)
    {
        LogicError("Wrong size for circulant_row_lzz_pX right trunc multiplication.");
    }

    Vec<zz_pX> data_trunc;
    trunc(data_trunc, dataRev, s);
    Vec<zz_pX> input_trunc;
    trunc(input_trunc, in, s);

    long dAx = deg(data_trunc);
    long d = deg(input_trunc);

    if (dAx == -1 || d == -1)
    {
        out.SetLength(n);
        for (long i = 0; i < n; i++)
            out[i] = 0;
        return;
    }


    Vec<zz_pX> res;
    zz_pX kro_in, kro_A, kro_res;
    to_kronecker(kro_in, input_trunc, dAx + d);
    to_kronecker(kro_A, data_trunc, dAx + d);
    kro_res = kro_in * kro_A;
    from_kronecker(res, kro_res, dAx + d);
    long ell = res.length();
    res.SetLength(m + m);
    for (long i = 0; i < ell; i++)
        trunc(res[i], res[i], s);
    for (long i = ell; i < m + m; i++)
        res[i] = 0;

    out.SetLength(n);
    long todo = min(n, m);
    for (long i = 0; i < todo; i++)
        out[i] = res[i] + res[i + m];

    if (m < n)
    {
        long done = m;
        long left = n-m;
        do 
        {
            todo = min(left, m);
            for (long i = done; i < done + todo; i++)
                out[i] = out[i-m];
            done += todo;
            left -= todo;
        } 
        while (left >= 1);
    }
}

/*------------------------------------------------------------*/
/* left multiplication mod x^s                                */
/*------------------------------------------------------------*/
void circulant_row_lzz_pX::mul_left_trunc(Vec<zz_pX>& out, const Vec<zz_pX>& in, long s) const
{
    if (&out == &in)
    {
        out = mul_left_trunc(in, s);
        return;
    }

    if (in.length() != n)
    {
        LogicError("Wrong size for circulant_row_lzz_pX left trunc multiplication.");
    }
    out.SetLength(m);

    Vec<zz_pX> data_trunc; 
    mirror(data_trunc, data);
    trunc(data_trunc, data_trunc, s);
    Vec<zz_pX> input_trunc;
    trunc(input_trunc, in, s);
    long dAx = deg(data_trunc);
    long d = deg(input_trunc);

    if (dAx == -1 || d == -1)
    {
        for (long i = 0; i < m; i++)
            out[i] = 0;
        return;
    }

    long leftover = n % m;
    long nb = n / m;   // number of size-m chunks

    if (leftover == 0)
    {
        leftover = m;
        nb--;
    }


    Vec<zz_pX> res;
    zz_pX kro_in, kro_A, kro_res;
    Vec<zz_pX> input_slice;
    input_slice.SetLength(leftover);
    for(long i = 0; i < leftover; i++)
        input_slice[i] = input_trunc[i + nb * m];

    to_kronecker(kro_in, input_slice, dAx + d);
    to_kronecker(kro_A, data_trunc, dAx + d);
    kro_res = kro_in * kro_A;
    from_kronecker(res, kro_res, dAx + d);
    long ell = res.length();
    res.SetLength(m + leftover);
    for (long i = 0; i < ell; i++)
        trunc(res[i], res[i], s);
    for (long i = ell; i < m + leftover; i++)
        res[i] = 0;

    for (long i = 0; i < leftover; i++)
        out[i] = res[i] + res[i + m];
    for (long i = leftover; i < m; i++)
        out[i] = res[i];

    for (long r = 0; r < nb; r++)
    {
        input_slice.SetLength(m);
        for (long i = 0; i < m; i++)
            input_slice[i] = input_trunc[i + r * m];
        to_kronecker(kro_in, input_slice, dAx + d);
        kro_res = kro_in * kro_A;
        from_kronecker(res, kro_res, dAx + d);
        long ell = res.length();
        res.SetLength(m + m);
        for (long i = 0; i < ell; i++)
            trunc(res[i], res[i], s);
        for (long i = ell; i < m + m; i++)
            res[i] = 0;
    
        for (long i = 0; i < m; i++)
            out[i] += (res[i] + res[i + m]);
    }
}


// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
