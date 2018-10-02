#include <NTL/lzz_pX.h>
#include <NTL/mat_lzz_p.h>

#include "lzz_pX_middle_product.h"
#include "mosaic_hankel_lzz_p.h"

NTL_CLIENT

/*----------------------------------------------------*/
/*----------------------------------------------------*/
/* Hankel matrices                                    */
/* stored as                                          */
/*       a5 a4 a3 a2                                  */
/*       a4 a3 a2 a1                                  */
/*       a3 a2 a1 a0                                  */
/*----------------------------------------------------*/
/*----------------------------------------------------*/

/*----------------------------------------------------*/
/* sets dimensions to 0                               */
/*----------------------------------------------------*/
hankel_lzz_p::hankel_lzz_p()
{
    data.SetLength(0);
    n = m = 0;
}

/*----------------------------------------------------*/
/* input vector is as showed above                    */
/*----------------------------------------------------*/
hankel_lzz_p::hankel_lzz_p(const Vec<zz_p>& input, long rows, long cols)
{
    n = rows;
    m = cols;
    data = input;
    data_rev.SetLength(n+m-1);
    zz_pX data_X;
    data_X.rep.SetLength(n+m-1);
    
    for (long i = 0;i < n+m-1; i++){
        data_rev[i] = input[n+m-2-i];
        data_X.rep[i] = data_rev[i];
    }    
    data_X.normalize();
    
    long K = NextPowerOfTwo(n+m-1);
    fft_data = fftRep(INIT_SIZE, K);
    TofftRep(fft_data, data_X, K);
}

/*----------------------------------------------------*/
/* getters                                            */
/*----------------------------------------------------*/
long hankel_lzz_p::NumCols() const
{
    return m;
}

long hankel_lzz_p::NumRows() const
{
    return n;
}


/*----------------------------------------------------*/
/* turns M into a dense matrix                        */
/*----------------------------------------------------*/
void hankel_lzz_p::to_dense(Mat<zz_p>& Mdense)
{
    Mdense.SetDims(n, m);
    for (long i = 0; i < n; i++)
        for (long j = 0; j < m; j++)
            Mdense[i][j] = data[m+n-2-i-j];
}

/*----------------------------------------------------*/
/* right multiplication                               */
/*----------------------------------------------------*/
void hankel_lzz_p::mul_right(Vec<zz_p>& res, const Vec<zz_p>& input)
{
    if (input.length() != m)
    {
        LogicError("Wrong size for hankel_lzz_p right multiplication.");
    }

    if (&res == &input)
    {
        res = mul_right(input);
        return;
    }

    res.SetLength(n);

    if (min(n, m) <= NTL_zz_pX_MUL_CROSSOVER)
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
}


/*----------------------------------------------------*/
/* right multiplication                               */
/*----------------------------------------------------*/
void hankel_lzz_p::mul_right(Mat<zz_p>& res, const Mat<zz_p>& input)
{
    if (input.NumRows() != m)
    {
        LogicError("Wrong size for hankel_lzz_p right matrix multiplication.");
    }

    if (&res == &input)
    {
        res = mul_right(input);
        return;
    }

    long p = input.NumCols();
    res.SetDims(n, p);

    Vec<zz_p> res_vec;
    res_vec.SetLength(n);

    if (min(n, m) <= NTL_zz_pX_MUL_CROSSOVER)
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
                res[j][i] = res_vec[j];
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
                res[i][j] = res_vec[i];
        }
    }
}



// /*----------------------------------------------------*/
// /* left multiplication                                */
// /*----------------------------------------------------*/
// void mul_left(Vec<zz_p>& res, const hankel_lzz_p& M, const Vec<zz_p>& input){

//     long nM = M.NumRows();
//     long mM = M.NumCols();

//     res.SetLength(mM);
//     if (min(nM, mM) <= NTL_zz_pX_MUL_CROSSOVER){
//         long sp = Kar_stk_size(max(nM, mM));
//         zz_p *stk = new zz_p[sp];
//         tKarMul_aux(res._vec__rep.rep, mM, input._vec__rep.rep, nM, M.data_rev._vec__rep.rep, nM+mM-1, stk);
//         delete[] stk;
//     }
//     else{
//         long K = NextPowerOfTwo(nM+mM-1);
//         fftRep fft_input = fftRep(INIT_SIZE, K);

//         zz_pX input_X, output_X;
//         input_X.rep.SetLength(nM);
//         zz_p *cf = input_X.rep.elts();
//         for (long i = 0; i < nM; i++)
//             cf[i] = input[nM-1-i];
//         input_X.normalize();

//         TofftRep(fft_input, input_X, K);
//         mul(fft_input, fft_input, M.fft_data);
//         FromfftRep(res._vec__rep.rep, fft_input, nM-1, nM+mM-2);
//     }
// }


// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
