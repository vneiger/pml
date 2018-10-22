#include <NTL/lzz_pX.h>
#include <NTL/mat_lzz_p.h>

#include "lzz_p_extra.h"
#include "mat_lzz_p_extra.h"
#include "structured_lzz_p.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* Hankel like matrices, where generators G, H are such that  */
/* Z1 M - M Z0^t = G H^t                                      */
/* -> M = sum_i circ(g_i) L(h_i) J                            */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* sets dimensions to 0                                       */
/*------------------------------------------------------------*/
hankel_like_plus_lzz_p::hankel_like_plus_lzz_p()
{
    m = 0;
    n = 0;
}

/*------------------------------------------------------------*/
/* input vector is as showed above, with G<-U, H<-V           */
/*------------------------------------------------------------*/
hankel_like_plus_lzz_p::hankel_like_plus_lzz_p(const Mat<zz_p>& U, const Mat<zz_p>& V)
{
    G = U;
    H = V;
    m = G.NumRows();
    n = H.NumRows();
    long alpha = G.NumCols();

    toeplitz_G.SetLength(alpha);
    hankel_H.SetLength(alpha);
    
    Vec<zz_p> vecG, vecH;
    vecG.SetLength(m + n - 1);
    vecH.SetLength(n + n - 1);

    for (long i = 0; i < alpha; i++)
    {
        for (long j = n, jm = n % m; j < m + n; j++)
        {
            vecG[j - 1] = G[jm][i];
            jm++;
            if (jm == m)
                jm = 0;
        }

        if (m > 0 && n > 0)
            for (long j = 1, c = (n - 1) % m; j < n; j++)
            {
                vecG[n - 1 - j] = G[c][i];
                c--;
                if (c < 0)
                    c = m - 1;
            }

        toeplitz_G[i] = toeplitz_lzz_p(vecG, m, n);

        for (long j = 0; j < n; j++)
            vecH[j] = H[n - 1 - j][i];
        for (long j = n; j < 2*n - 1; j++)
            vecH[j] = 0;
        
        hankel_H[i] = hankel_lzz_p(vecH, n , n);
    }
}

/*------------------------------------------------------------*/
/* getters                                                    */
/*------------------------------------------------------------*/
long hankel_like_plus_lzz_p::NumRows() const
{
    return m;
}

long hankel_like_plus_lzz_p::NumCols() const
{
    return n;
}

long hankel_like_plus_lzz_p::NumGens() const
{
    return G.NumCols();
}

/*------------------------------------------------------------*/
/* output = M * input                                         */
/*------------------------------------------------------------*/
void hankel_like_plus_lzz_p::mul_right(Vec<zz_p>& output, const Vec<zz_p>& input) const
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
        hankel_H[i].mul_right(tmp1, input);
        toeplitz_G[i].mul_right(tmp2, tmp1);
        for (long j = 0; j < m; j++)
            output[j] += tmp2[j];
    }
}

void hankel_like_plus_lzz_p::mul_right(Mat<zz_p>& out, const Mat<zz_p>& in) const
{
}

/*------------------------------------------------------------*/
/* left multiplication                                        */
/*------------------------------------------------------------*/
void hankel_like_plus_lzz_p::mul_left(Vec<zz_p>& out, const Vec<zz_p>& in) const
{
}

void hankel_like_plus_lzz_p::mul_left(Mat<zz_p>& out, const Mat<zz_p>& in) const
{
}

/*------------------------------------------------------------*/
/* turns M into a dense matrix                                */
/*------------------------------------------------------------*/
void hankel_like_plus_lzz_p::to_dense(Mat<zz_p>& Mdense) const 
{
    Mdense.SetDims(m, n);
    Vec<zz_p> vecG, vecH;
    vecG.SetLength(m);
    vecH.SetLength(n);

    for (long i = 0; i < NumGens(); i++)
    {
        for (long k = 0; k < m; k++)
            vecG[k] = G[(k + n) % m][i];
        for (long k = 0; k < n; k++)
            vecH[k] = H[k][i];
        Mdense += circulant_lzz_p(vecG, m, n) * lower_hankel_lzz_p(vecH, n, n);
    }
}

/*------------------------------------------------------------*/
/* returns Z1 A - A Z0^t                                      */
/*------------------------------------------------------------*/
void hankel_lzz_p_phi_plus(Mat<zz_p> & res, const Mat<zz_p>& A)
{
    if (&res == &A)
    {
        res = hankel_lzz_p_phi_plus(A);
        return;
    }

    long m = A.NumRows();
    long n = A.NumCols();
    
    Mat<zz_p> Tz0, z0;
    z0 = Z_lzz_p(n, to_zz_p(0));
    transpose(Tz0, z0);
    res = Z_lzz_p(m, to_zz_p(1)) * A - A * Tz0;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
