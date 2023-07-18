#include <NTL/lzz_pX.h>
#include <NTL/mat_lzz_p.h>

#include "lzz_p_extra.h"
#include "mat_lzz_p_extra.h"
#include "structured_lzz_p.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* Toeplitz like matrices, where generators G, H are such that*/
/* Z1 M - M Z0 = G H^t                                        */
/* -> M = sum_i circ(g_i) L(rev h_i)                          */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* sets dimensions to 0                                       */
/*------------------------------------------------------------*/
toeplitz_like_plus_lzz_p::toeplitz_like_plus_lzz_p()
{
    m = 0;
    n = 0;
}

/*------------------------------------------------------------*/
/* input vector is as showed above, with G<-U, H<-V           */
/*------------------------------------------------------------*/
toeplitz_like_plus_lzz_p::toeplitz_like_plus_lzz_p(const Mat<zz_p>& U, const Mat<zz_p>& V)
{
    G = U;
    H = V;
    m = G.NumRows();
    n = H.NumRows();
    long alpha = G.NumCols();

    circulant_G.SetLength(alpha);
    lower_triangular_toeplitz_H.SetLength(alpha);
    
    Vec<zz_p> vecG, vecH;
    vecG.SetLength(m);
    vecH.SetLength(n);

    for (long i = 0; i < alpha; i++)
    {
        long idx = n % m;
        for (long j = 0; j < m; j++)
        {
            vecG[j] = G[idx][i];
            idx++;
            if (idx == m)
                idx = 0;
        }
        circulant_G[i] = circulant_column_lzz_p(vecG, n);
        
        for (long j = 0; j < n; j++)
            vecH[j] = H[n - 1 - j][i];
        lower_triangular_toeplitz_H[i] = lower_triangular_toeplitz_lzz_p(vecH);
    }
}

/*------------------------------------------------------------*/
/* getters                                                    */
/*------------------------------------------------------------*/
long toeplitz_like_plus_lzz_p::NumRows() const
{
    return m;
}

long toeplitz_like_plus_lzz_p::NumCols() const
{
    return n;
}

long toeplitz_like_plus_lzz_p::NumGens() const
{
    return G.NumCols();
}

/*------------------------------------------------------------*/
/* output = M * input                                         */
/*------------------------------------------------------------*/
void toeplitz_like_plus_lzz_p::mul_right(Vec<zz_p>& output, const Vec<zz_p>& input) const
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
        lower_triangular_toeplitz_H[i].mul_right(tmp1, input);
        circulant_G[i].mul_right(tmp2, tmp1);
        for (long j = 0; j < m; j++)
            output[j] += tmp2[j];
    }
}

void toeplitz_like_plus_lzz_p::mul_right(Mat<zz_p>& out, const Mat<zz_p>& in) const
{
}

/*------------------------------------------------------------*/
/* left multiplication                                        */
/*------------------------------------------------------------*/
void toeplitz_like_plus_lzz_p::mul_left(Vec<zz_p>& out, const Vec<zz_p>& in) const
{
}

void toeplitz_like_plus_lzz_p::mul_left(Mat<zz_p>& out, const Mat<zz_p>& in) const
{
}

/*------------------------------------------------------------*/
/* turns M into a dense matrix                                */
/*------------------------------------------------------------*/
void toeplitz_like_plus_lzz_p::to_dense(Mat<zz_p>& Mdense) const 
{
    Mdense.SetDims(m, n);
    for (long i = 0; i < NumGens(); i++)
    {
        Mdense += circulant_G[i].to_dense() * lower_triangular_toeplitz_H[i].to_dense();
    }
}

/*------------------------------------------------------------*/
/* returns Z1 A - A Z0                                        */
/*------------------------------------------------------------*/
void toeplitz_lzz_p_phi_plus(Mat<zz_p> & res, const Mat<zz_p>& A)
{
    if (&res == &A)
    {
        res = toeplitz_lzz_p_phi_plus(A);
        return;
    }

    long m = A.NumRows();
    long n = A.NumCols();
    res = Z_lzz_p(m, to_zz_p(1)) * A - A * Z_lzz_p(n, to_zz_p(0));
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
