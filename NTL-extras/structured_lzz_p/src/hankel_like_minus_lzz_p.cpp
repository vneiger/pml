#include <NTL/lzz_pX.h>
#include <NTL/mat_lzz_p.h>

#include "lzz_p_extra.h"
#include "mat_lzz_p_extra.h"
#include "structured_lzz_p.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* Hankel like matrices, where generators G, H are such that  */
/* Z0^t M - M Z1 = G H^t                                      */
/* -> M = sum_i U(-g_i,m,m) horizontal_circ (shift h_i,-1)    */
/* (dual structure to hankel_like_plus_lzz_p)                 */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* sets dimensions to 0                                       */
/*------------------------------------------------------------*/
hankel_like_minus_lzz_p::hankel_like_minus_lzz_p()
{
    m = 0;
    n = 0;
}

/*------------------------------------------------------------*/
/* input vector is as showed above, with G<-U, H<-V           */
/*------------------------------------------------------------*/
hankel_like_minus_lzz_p::hankel_like_minus_lzz_p(const Mat<zz_p>& U, const Mat<zz_p>& V)
{
    G = U;
    H = V;
    m = G.NumRows();
    n = H.NumRows();
    long alpha = G.NumCols();

    hankel_G.SetLength(alpha);
    toeplitz_H.SetLength(alpha);

    if (m == 0 || n == 0)
        return;
    
    Vec<zz_p> vecG, vecH;
    vecG.SetLength(m + m - 1);
    vecH.SetLength(m + n - 1);

    for (long i = 0; i < alpha; i++)
    {
        long idx;

        idx = 0;
        for (long j = 0; j < m - 1; j++, idx++)
            vecG[j] = 0;
        for (long j = 0; j < m; j++, idx++)
            vecG[idx] = -G[m - 1 - j][i];
        hankel_G[i] = hankel_lzz_p(vecG, m, m);

        idx = n - 2;
        if (idx < 0) // happens when n = 1;
            idx += n;
        for (long j = 0; j < m + n - 1; j++)
        {
            vecH[j] = H[idx][i];
            idx--;
            if (idx < 0)
                idx += n;
        }
        toeplitz_H[i] = toeplitz_lzz_p(vecH, m, n);
    }
}

/*------------------------------------------------------------*/
/* getters                                                    */
/*------------------------------------------------------------*/
long hankel_like_minus_lzz_p::NumRows() const
{
    return m;
}

long hankel_like_minus_lzz_p::NumCols() const
{
    return n;
}

long hankel_like_minus_lzz_p::NumGens() const
{
    return G.NumCols();
}

/*------------------------------------------------------------*/
/* output = M * input                                         */
/*------------------------------------------------------------*/
void hankel_like_minus_lzz_p::mul_right(Vec<zz_p>& output, const Vec<zz_p>& input) const
{
    if (&output == &input)
    {
        output = mul_right(input);
        return;
    }

    output.SetLength(m);
    for (long i = 0; i < m; i++)
        output[i] = 0;

    Vec<zz_p> tmp1, tmp2;
    for (long i = 0; i < NumGens(); i++)
    {
        toeplitz_H[i].mul_right(tmp1, input);
        hankel_G[i].mul_right(tmp2, tmp1);
        for (long j = 0; j < m; j++)
            output[j] += tmp2[j];
    }
}

void hankel_like_minus_lzz_p::mul_right(Mat<zz_p>& out, const Mat<zz_p>& in) const
{
}

/*------------------------------------------------------------*/
/* left multiplication                                        */
/*------------------------------------------------------------*/
void hankel_like_minus_lzz_p::mul_left(Vec<zz_p>& out, const Vec<zz_p>& in) const
{
}

void hankel_like_minus_lzz_p::mul_left(Mat<zz_p>& out, const Mat<zz_p>& in) const
{
}

/*------------------------------------------------------------*/
/* turns M into a dense matrix                                */
/*------------------------------------------------------------*/
void hankel_like_minus_lzz_p::to_dense(Mat<zz_p>& Mdense) const 
{
    Mdense.SetDims(m, n);
    for (long i = 0; i < NumGens(); i++)
    {
        Mdense += hankel_G[i].to_dense() * toeplitz_H[i].to_dense();
    }
}

/*------------------------------------------------------------*/
/* returns Z0^t A - A Z1                                      */
/*------------------------------------------------------------*/
void hankel_lzz_p_phi_minus(Mat<zz_p> & res, const Mat<zz_p>& A)
{
    if (&res == &A)
    {
        res = hankel_lzz_p_phi_minus(A);
        return;
    }

    long m = A.NumRows();
    long n = A.NumCols();
    
    Mat<zz_p> Tz0, z0;
    z0 = Z_lzz_p(m, to_zz_p(0));
    transpose(Tz0, z0);
    res = Tz0 * A - A * Z_lzz_p(n, to_zz_p(1));
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
