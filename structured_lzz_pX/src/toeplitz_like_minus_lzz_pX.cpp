#include <NTL/lzz_pX.h>
#include <NTL/mat_lzz_p.h>

#include "util.h"
#include "lzz_p_extra.h"
#include "mat_lzz_p_extra.h"
#include "mat_lzz_pX_extra.h"
#include "structured_lzz_p.h"
#include "structured_lzz_pX.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* Toeplitz like matrices, where generators G, H are such that*/
/* Z0 M - M Z1 = G H^t                                        */
/* -> M = sum_i L(-g_i) circ (rev h_i)                        */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* sets dimensions to 0                                       */
/*------------------------------------------------------------*/
toeplitz_like_minus_lzz_pX::toeplitz_like_minus_lzz_pX()
{
    m = 0;
    n = 0;
}

/*------------------------------------------------------------*/
/* input vector is as showed above, with G<-U, H<-V           */
/*------------------------------------------------------------*/
toeplitz_like_minus_lzz_pX::toeplitz_like_minus_lzz_pX(const Mat<zz_pX>& U, const Mat<zz_pX>& V)
{
    G = U;
    H = V;
    m = G.NumRows();
    n = H.NumRows();
    long alpha = G.NumCols();

    toeplitz_G.SetLength(alpha);
    circulant_H.SetLength(alpha);

    if (m == 0 || n == 0)
        return;
    
    Vec<zz_pX> vecG, vecH;
    vecG.SetLength(m);
    vecH.SetLength(n);

    for (long i = 0; i < alpha; i++)
    {
        long idx;
        for (long j = 0; j < m; j++)
            vecG[j] = -G[j][i];
        toeplitz_G[i] = lower_triangular_toeplitz_lzz_pX(vecG);

        idx = n - 2;
        if (idx < 0) // happens when n = 1;
            idx += n;
        for (long j = 0; j < n; j++)
        {
            vecH[j] = H[idx][i];
            idx--;
            if (idx < 0)
                idx += n;
        }
        circulant_H[i] = circulant_row_lzz_pX(vecH, m);
    }
}

/*------------------------------------------------------------*/
/* getters                                                    */
/*------------------------------------------------------------*/
long toeplitz_like_minus_lzz_pX::NumRows() const
{
    return m;
}

long toeplitz_like_minus_lzz_pX::NumCols() const
{
    return n;
}

long toeplitz_like_minus_lzz_pX::NumGens() const
{
    return G.NumCols();
}

/*------------------------------------------------------------*/
/* output = M * input                                         */
/*------------------------------------------------------------*/
void toeplitz_like_minus_lzz_pX::mul_right(Vec<zz_pX>& output, const Vec<zz_pX>& input) const
{
    if (&output == &input)
    {
        output = mul_right(input);
        return;
    }

    output.SetLength(m);
    for (long i = 0; i < m; i++)
        output[i] = 0;

    Vec<zz_pX> tmp1, tmp2;

    for (long i = 0; i < NumGens(); i++)
    {
        circulant_H[i].mul_right(tmp1, input);
        toeplitz_G[i].mul_right(tmp2, tmp1);
        for (long j = 0; j < m; j++)
            output[j] += tmp2[j];
    }
}

/*------------------------------------------------------------*/
/* left multiplication                                        */
/*------------------------------------------------------------*/
void toeplitz_like_minus_lzz_pX::mul_left(Vec<zz_pX>& output, const Vec<zz_pX>& input) const
{
    if (&output == &input)
    {
        output = mul_left(input);
        return;
    }

    output.SetLength(n);
    for (long i = 0; i < n; i++)
        output[i] = 0;

    Vec<zz_pX> tmp1, tmp2;

    for (long i = 0; i < NumGens(); i++)
    {
        toeplitz_G[i].mul_left(tmp1, input);
        circulant_H[i].mul_left(tmp2, tmp1);
        for (long j = 0; j < n; j++)
            output[j] += tmp2[j];
    }
}

/*------------------------------------------------------------*/
/* output = M * input mod x^d                                 */
/*------------------------------------------------------------*/
void toeplitz_like_minus_lzz_pX::mul_right_trunc(Vec<zz_pX>& output, const Vec<zz_pX>& input, long d) const
{
    if (&output == &input)
    {
        output = mul_right_trunc(input, d);
        return;
    }

    output.SetLength(m);
    for (long i = 0; i < m; i++)
        output[i] = 0;

    Vec<zz_pX> tmp1, tmp2;

    for (long i = 0; i < NumGens(); i++)
    {
        circulant_H[i].mul_right_trunc(tmp1, input, d);
        toeplitz_G[i].mul_right_trunc(tmp2, tmp1, d);
        for (long j = 0; j < m; j++)
            output[j] += tmp2[j];
    }
}

/*------------------------------------------------------------*/
/* left multiplication mod x^d                                */
/*------------------------------------------------------------*/
void toeplitz_like_minus_lzz_pX::mul_left_trunc(Vec<zz_pX>& output, const Vec<zz_pX>& input, long d) const
{
    if (&output == &input)
    {
        output = mul_left_trunc(input, d);
        return;
    }

    output.SetLength(n);
    for (long i = 0; i < n; i++)
        output[i] = 0;

    Vec<zz_pX> tmp1, tmp2;

    for (long i = 0; i < NumGens(); i++)
    {
        toeplitz_G[i].mul_left_trunc(tmp1, input, d);
        circulant_H[i].mul_left_trunc(tmp2, tmp1, d);
        for (long j = 0; j < n; j++)
            output[j] += tmp2[j];
    }
}


/*------------------------------------------------------------*/
/* turns M into a dense matrix                                */
/*------------------------------------------------------------*/
void toeplitz_like_minus_lzz_pX::to_dense(Mat<zz_pX>& Mdense) const 
{
    Mdense.SetDims(m, n);
    for (long i = 0; i < NumGens(); i++)
    {
        Mdense += toeplitz_G[i].to_dense() * circulant_H[i].to_dense();
    }
}

/*------------------------------------------------------------*/
/* returns Z0 A - A Z1                                        */
/*------------------------------------------------------------*/
void toeplitz_lzz_pX_phi_minus(Mat<zz_pX> & res, const Mat<zz_pX>& A)
{
    if (&res == &A)
    {
        res = toeplitz_lzz_pX_phi_minus(A);
        return;
    }

    long m = A.NumRows();
    long n = A.NumCols();
    
    Mat<zz_p> z0;
    z0 = Z_lzz_p(m, to_zz_p(0));
    Mat<zz_p> z1;
    z1 = Z_lzz_p(n, to_zz_p(1));
    res = z0 * A - A * z1;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
