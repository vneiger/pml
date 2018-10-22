#include <NTL/lzz_pX.h>
#include <NTL/mat_lzz_p.h>

#include "lzz_p_extra.h"
#include "lzz_pX_middle_product.h"
#include "structured_lzz_p.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* Toeplitz matrices                                          */
/* stored as                                                  */
/*       a3 a2 a1 a0                                          */
/*       a4 a3 a2 a1                                          */
/*       a5 a4 a3 a2                                          */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* sets dimensions to 0                                       */
/*------------------------------------------------------------*/
toeplitz_lzz_p::toeplitz_lzz_p()
{
    n = m = 0;
}

/*------------------------------------------------------------*/
/* input vector is as showed above                            */
/*------------------------------------------------------------*/
toeplitz_lzz_p::toeplitz_lzz_p(const Vec<zz_p>& input, long rows, long cols)
{
    n = rows;
    m = cols;
    data = input;
    data_rev.SetLength(n+m-1);
    zz_pX data_X;
    data_X.rep.SetLength(n+m-1);
    
    for (long i = 0; i < n+m-1; i++){
        data_rev[i] = input[n+m-2-i];
        data_X.rep[i] = data_rev[i];
    }    
    data_X.normalize();
    
    long K = NextPowerOfTwo(n+m-1);
    fft_data = fftRep(INIT_SIZE, K);
    TofftRep(fft_data, data_X, K);
}

/*------------------------------------------------------------*/
/* getters                                                    */
/*------------------------------------------------------------*/
long toeplitz_lzz_p::NumCols() const
{
    return m;
}

long toeplitz_lzz_p::NumRows() const
{
    return n;
}

/*------------------------------------------------------------*/
/* turns M into a dense matrix                                */
/*------------------------------------------------------------*/
void toeplitz_lzz_p::to_dense(Mat<zz_p>& Mdense) const 
{
    Mdense.SetDims(n, m);
    for (long i = 0; i < n; i++)
        for (long j = 0; j < m; j++)
            Mdense[i][j] = data[m+i-1-j];
}

/*------------------------------------------------------------*/
/* right multiplication                                       */
/*------------------------------------------------------------*/
void toeplitz_lzz_p::mul_right(Vec<zz_p>& res, const Vec<zz_p>& input) const
{
    if (input.length() != m)
    {
        LogicError("Wrong size for toeplitz_lzz_p right multiplication.");
    }

    if (&res == &input)
    {
        res = mul_right(input);
        return;
    }

    res.SetLength(n);

    if (min(n, m) <= (2*NTL_zz_pX_MUL_CROSSOVER)/4)
    {
        long sp = Kar_stk_size(max(n, m));
        Vec<zz_p> stk;
        stk.SetLength(sp);
        tKarMul_aux(res.elts(), n, input.elts(), m, data_rev.elts(), n+m-1, stk.elts());
    }
    else
    {
        long K = NextPowerOfTwo(n+m-1);
        fftRep fft_input = fftRep(INIT_SIZE, K);

        zz_pX input_X, output_X;
        input_X.rep.SetLength(m);
        zz_p *cf = input_X.rep.elts();
        for (long i = 0; i < m; i++)
            cf[i] = input[m-1-i];
        input_X.normalize();

        TofftRep(fft_input, input_X, K);
        mul(fft_input, fft_input, fft_data);
        FromfftRep(res.elts(), fft_input, m-1, n+m-2);
    }

    for (long i = 0; (i+i) <= n-1; i++)
    {
        long ri = res[i].LoopHole();
        res[i].LoopHole() = res[n-1-i].LoopHole();
        res[n-1-i].LoopHole() = ri;
    }

}


/*------------------------------------------------------------*/
/* right matrix multiplication                                */
/*------------------------------------------------------------*/
void toeplitz_lzz_p::mul_right(Mat<zz_p>& res, const Mat<zz_p>& input) const
{
    if (input.NumRows() != m)
    {
        LogicError("Wrong size for toeplitz_lzz_p right matrix multiplication.");
    }

    if (&res == &input)
    {
        res = mul_right(input);
        return;
    }

    long p = input.NumCols();

    if ( (min(m, n) >= max(m, n)/2) && // close enough to a square matrix
         p >= min(m, n) / 4 )          // close enough to a matrix-matrix product
    {
        long t = type_of_prime();
        if ( (t == TYPE_FFT_PRIME && min(m, n) <= 100)  ||
             (t == TYPE_SMALL_PRIME && min(m, n) <= 1000) ||
             (t == TYPE_LARGE_PRIME && min(m, n) <= 300) )
        {
            Mat<zz_p> hank = to_dense();
            res = hank * input;
            return;
        }
    }

    res.SetDims(n, p);

    Vec<zz_p> res_vec;
    res_vec.SetLength(n);

    if (min(n, m) <= (2*NTL_zz_pX_MUL_CROSSOVER)/4)
    {
        Vec<zz_p> in_vec;
        in_vec.SetLength(m);
        long sp = Kar_stk_size(max(n, m));
        Vec<zz_p> stk;
        stk.SetLength(sp);

        for (long i = 0; i < p; i++)
        {
            for (long j = 0; j < m; j++)
                in_vec[j] = input[j][i];
            tKarMul_aux(res_vec.elts(), n, in_vec.elts(), m, data_rev.elts(), n+m-1, stk.elts());
            for (long j = 0; j < n; j++)
                res[n-1-j][i] = res_vec[j];
        }
    }
    else
    {
        long K = NextPowerOfTwo(n+m-1);
        fftRep fft_input = fftRep(INIT_SIZE, K);
        zz_pX input_X, output_X;

        for (long j = 0; j < p; j++)
        {
            input_X.SetLength(m);
            zz_p *cf = input_X.rep.elts();
            for (long i = 0; i < m; i++)
                cf[i] = input[m-1-i][j];
            input_X.normalize();
            
            TofftRep(fft_input, input_X, K);
            mul(fft_input, fft_input, fft_data);
            FromfftRep(res_vec.elts(), fft_input, m-1, n+m-2);
            for (long i = 0; i < n; i++)
                res[n-1-i][j] = res_vec[i];
        }
    }
}


/*------------------------------------------------------------*/
/* left multiplication                                        */
/*------------------------------------------------------------*/
void toeplitz_lzz_p::mul_left(Vec<zz_p>& res, const Vec<zz_p>& input) const
{
    if (input.length() != n)
    {
        LogicError("Wrong size for toeplitz_lzz_p left matrix multiplication.");
    }

    if (&res == &input)
    {
        res = mul_left(input);
        return;
    }

    res.SetLength(m);
    if (min(n, m) <= (2*NTL_zz_pX_MUL_CROSSOVER)/4)
    {
        Vec<zz_p> in_rev;
        in_rev.SetLength(n);
        for (long i = 0; i < n; i++)
            in_rev[i] = input[n-1-i];

        long sp = Kar_stk_size(max(n, m));
        Vec<zz_p> stk;
        stk.SetLength(sp);
        tKarMul_aux(res.elts(), m, in_rev.elts(), n, data_rev.elts(), n+m-1, stk.elts());
    }
    else
    {
        long K = NextPowerOfTwo(n+m-1);
        fftRep fft_input = fftRep(INIT_SIZE, K);

        zz_pX input_X, output_X;
        input_X.rep.SetLength(n);
        zz_p *cf = input_X.rep.elts();
        for (long i = 0; i < n; i++)
            cf[i] = input[i];
        input_X.normalize();

        TofftRep(fft_input, input_X, K);
        mul(fft_input, fft_input, fft_data);
        FromfftRep(res.elts(), fft_input, n-1, n+m-2);
    }
}

/*------------------------------------------------------------*/
/* left matrix multiplication                                 */
/*------------------------------------------------------------*/
void toeplitz_lzz_p::mul_left(Mat<zz_p>& res, const Mat<zz_p>& input) const
{
    if (input.NumCols() != n)
    {
        LogicError("Wrong size for toeplitz_lzz_p left matrix multiplication.");
    }

    if (&res == &input)
    {
        res = mul_left(input);
        return;
    }

    long p = input.NumRows();

    if ( (min(m, n) >= max(m, n)/2) && // close enough to a square matrix
         p >= min(m, n) / 4 )          // close enough to a matrix-matrix product
    {
        long t = type_of_prime();
        if ( (t == TYPE_FFT_PRIME && min(m, n) <= 100)  ||
             (t == TYPE_SMALL_PRIME && min(m, n) <= 1000) ||
             (t == TYPE_LARGE_PRIME && min(m, n) <= 300) )
        {
            Mat<zz_p> hank = to_dense();
            res = input * hank;
            return;
        }
    }

    res.SetDims(p, m);

    Vec<zz_p> res_vec;
    res_vec.SetLength(m);

    if (min(n, m) <= (2*NTL_zz_pX_MUL_CROSSOVER)/4)
    {
        Vec<zz_p> in_vec;
        in_vec.SetLength(n);
        long sp = Kar_stk_size(max(n, m));
        Vec<zz_p> stk;
        stk.SetLength(sp);
        for (long i = 0; i < p; i++)
        {
            for (long j = 0; j < n; j++)
                in_vec[j] = input[i][n-1-j];
            tKarMul_aux(res_vec.elts(), m, in_vec.elts(), n, data_rev.elts(), n+m-1, stk.elts());
            for (long j = 0; j < m; j++)
                res[i][j] = res_vec[j];
        }
    }
    else
    {
        long K = NextPowerOfTwo(n+m-1);
        fftRep fft_input = fftRep(INIT_SIZE, K);
        zz_pX input_X, output_X;

        for (long i = 0; i < p; i++)
        {
            input_X.rep.SetLength(n);
            zz_p *cf = input_X.rep.elts();
            for (long j = 0; j < n; j++)
                cf[j] = input[i][j];
            input_X.normalize();

            TofftRep(fft_input, input_X, K);
            mul(fft_input, fft_input, fft_data);
            FromfftRep(res_vec.elts(), fft_input, n-1, n+m-2);
            for (long j = 0; j < m; j++)
                res[i][j] = res_vec[j];
        }
    }
}


// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
