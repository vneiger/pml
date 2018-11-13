#include <NTL/lzz_pX.h>
#include <NTL/mat_lzz_p.h>

#include "lzz_p_extra.h"
#include "mat_lzz_p_extra.h"
#include "lzz_pX_middle_product.h"
#include "mat_lzz_pX_extra.h"
#include "structured_lzz_p.h"
#include "structured_lzz_pX.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/*                   SYLVESTER MATRICES                       */
/* sylv(A,B)= (square) matrix of (F,G) \mapsto AF + BG        */
/* with deg(F) < deg(B), deg(G) < deg(A)                      */
/* using monomial bases in increasing degrees for rows / cols */
/* particular case of mosaic toeplitz                         */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* sets dimensions to 0                                       */
/*------------------------------------------------------------*/
sylvester_lzz_pX::sylvester_lzz_pX()
{
    n = 0;
}

/*------------------------------------------------------------*/
/* dimension = deg(A) + deg(B), A and B as above              */
/*------------------------------------------------------------*/
sylvester_lzz_pX::sylvester_lzz_pX(const Vec<zz_pX>& A, const Vec<zz_pX>& B)
{
    if (A[A.length() - 1] == 0 || B[B.length() - 1] == 0)
        Error("Vanishing leading coefficient for sylvester_lzz_pX");
    dAy = A.length() - 1;
    dBy = B.length() - 1;
    dAx = deg(A);
    dBx = deg(B);
        
    a = A;
    b = B;
    revA.SetLength(a.length());
    for (long i = 0; i < a.length(); i++)
        revA[i] = a[a.length() - 1 - i];
    revB.SetLength(b.length());
    for (long i = 0; i < b.length(); i++)
        revB[i] = b[b.length() - 1 - i];

    n = dAy + dBy;
}

/*------------------------------------------------------------*/
/* getters                                                    */
/*------------------------------------------------------------*/
long sylvester_lzz_pX::NumCols() const
{
    return n;
}

long sylvester_lzz_pX::NumRows() const
{
    return n;
}

/*------------------------------------------------------------*/
/* right multiplication                                       */
/*------------------------------------------------------------*/
void sylvester_lzz_pX::mul_right(Vec<zz_pX>& out, const Vec<zz_pX>& in) const
{
    if (in.length() != n)
        Error("Bad input size for sylvester_lzz_pX mul right");

    zz_pX kro_A, kro_B, kro_F, kro_G;
    Vec<zz_pX> in_F, in_G;

    in_F.SetLength(dBy);
    in_G.SetLength(dAy);
    for (long i = 0; i < dBy; i++)
        in_F[i] = in[i];
    for (long i = 0; i < dAy; i++)
        in_G[i] = in[i + dBy];

    long dFx = deg(in_F);    
    long dGx = deg(in_G);

    to_kronecker(kro_A, a, dAx + dFx);  
    to_kronecker(kro_F, in_F, dAx + dFx); 
    to_kronecker(kro_B, b, dBx + dGx); 
    to_kronecker(kro_G, in_G, dBx + dGx); 

    if ((dAx + dFx) == (dBx + dGx))
    {
        zz_pX kro_res;
        kro_res = kro_A * kro_F + kro_B * kro_G;
        from_kronecker(out, kro_res, dAx + dFx);
        long ell = out.length();
        out.SetLength(n);
        for (long i = ell; i < n; i++)
            out[i] = 0;
    }
    else
    {
        zz_pX kro_res1, kro_res2;
        Vec<zz_pX> out1, out2;
        kro_res1 = kro_A * kro_F;
        from_kronecker(out1, kro_res1, dAx + dFx);
        kro_res2 = kro_B * kro_G;
        from_kronecker(out2, kro_res2, dBx + dGx);
        out.SetLength(n);
        for (long i = 0; i < n; i++)
        {
            out[i] = 0;
            if (i < out1.length())
                out[i] = out1[i];
            if (i < out2.length())
                out[i] = out[i] + out2[i];
        }
    }
}

/*------------------------------------------------------------*/
/* right multiplication, matrix version                       */
/*------------------------------------------------------------*/
void sylvester_lzz_pX::mul_right(Mat<zz_pX>& out, const Mat<zz_pX>& in) const
{
    if (&out == &in)
    {
        Mat<zz_pX> out2 = mul_right(in);
        out = out2;
        return;
    }

    if (in.NumRows() != n)
        Error("Bad dimensions in sylvester_lzz_pX matrix mul right");

    long p = in.NumCols();
    Vec<zz_pX> vec_in, vec_out;
    vec_in.SetLength(n);
    out.SetDims(n, p);

    for (long j = 0; j < p; j++)
    {
        for (long i = 0; i < n; i++)
            vec_in[i] = in[i][j];
        mul_right(vec_out, vec_in);
        for (long i = 0; i < n; i++)
            out[i][j] = vec_out[i];
    }
}

/*------------------------------------------------------------*/
/* right multiplication truncated mod x^s                     */
/*------------------------------------------------------------*/
void sylvester_lzz_pX::mul_right_trunc(Vec<zz_pX>& out, const Vec<zz_pX>& in, long s) const
{
    if (in.length() != n)
        Error("Bad input size for sylvester_lzz_pX mul right");

    zz_pX kro_A, kro_B, kro_F, kro_G;
    Vec<zz_pX> in_F, in_G;

    in_F.SetLength(dBy);
    in_G.SetLength(dAy);
    for (long i = 0; i < dBy; i++)
        trunc(in_F[i], in[i], s);
    for (long i = 0; i < dAy; i++)
        trunc(in_G[i], in[i + dBy], s);
    long dFx = deg(in_F);    
    long dGx = deg(in_G);

    Vec<zz_pX> truncA, truncB;
    truncA.SetLength(a.length());
    truncB.SetLength(b.length());
    for (long i = 0; i < a.length(); i++)
        trunc(truncA[i], a[i], s);
    for (long i = 0; i < b.length(); i++)
        trunc(truncB[i], b[i], s);
    long dTAx = deg(truncA);
    long dTBx = deg(truncB);

    to_kronecker(kro_A, a, dTAx + dFx);  
    to_kronecker(kro_F, in_F, dTAx + dFx); 
    to_kronecker(kro_B, b, dTBx + dGx); 
    to_kronecker(kro_G, in_G, dTBx + dGx); 

    if ((dTAx + dFx) == (dTBx + dGx))
    {
        zz_pX kro_res;
        kro_res = kro_A * kro_F + kro_B * kro_G;
        from_kronecker(out, kro_res, dTAx + dFx);
        long ell = out.length();
        out.SetLength(n);
        for (long i = 0; i < ell; i++)
            trunc(out[i], out[i], s);
        for (long i = ell; i < n; i++)
            out[i] = 0;
    }
    else
    {
        zz_pX kro_res1, kro_res2;
        Vec<zz_pX> out1, out2;
        kro_res1 = kro_A * kro_F;
        from_kronecker(out1, kro_res1, dTAx + dFx);
        kro_res2 = kro_B * kro_G;
        from_kronecker(out2, kro_res2, dTBx + dGx);
        out.SetLength(n);
        for (long i = 0; i < n; i++)
        {
            out[i] = 0;
            if (i < out1.length())
                trunc(out[i], out1[i], s);
            if (i < out2.length())
                out[i] = out[i] + trunc(out2[i], s);
        }
    }
}

/*------------------------------------------------------------*/
/* right multiplication mod x^s, matrix version               */
/*------------------------------------------------------------*/
void sylvester_lzz_pX::mul_right_trunc(Mat<zz_pX>& out, const Mat<zz_pX>& in, long s) const
{
    if (&out == &in)
    {
        Mat<zz_pX> out2 = mul_right_trunc(in, s);
        out = out2;
        return;
    }

    if (in.NumRows() != n)
        Error("Bad dimensions in sylvester_lzz_pX matrix mul right");

    long p = in.NumCols();
    Vec<zz_pX> vec_in, vec_out;
    vec_in.SetLength(n);
    out.SetDims(n, p);

    for (long j = 0; j < p; j++)
    {
        for (long i = 0; i < n; i++)
            vec_in[i] = in[i][j];
        mul_right_trunc(vec_out, vec_in, s);
        for (long i = 0; i < n; i++)
            out[i][j] = vec_out[i];
    }
}

/*------------------------------------------------------------*/
/* left multiplication                                        */
/*------------------------------------------------------------*/
void sylvester_lzz_pX::mul_left(Vec<zz_pX>& out, const Vec<zz_pX>& in) const
{
    if (in.length() != n)
        Error("Bad input size for sylvester_lzz_pX mul left");

    if (&out == &in)
    {
        Vec<zz_pX> out2 = mul_left(in);
        out = out2;
        return;
    }

    long d = deg(in);
    Vec<zz_pX> outA, outB;
    zz_pX kro_revA, kro_revB;
    to_kronecker(kro_revA, revA, dAx + d);  
    to_kronecker(kro_revB, revB, dBx + d);  
    
    if (dAx == dBx)
    {
        zz_pX kro_in, kro_outA, kro_outB;
        to_kronecker(kro_in, in, dAx + d);
        middle_product(kro_outA, kro_revA, kro_in << dAx, dAy * (dAx + d + 1) + dAx, dBy * (dAx + d + 1) - 1);
        middle_product(kro_outB, kro_revB, kro_in << dBx, dBy * (dAx + d + 1) + dBx, dAy * (dAx + d + 1) - 1);
        from_kronecker(outA, kro_outA, dAx + d);
        from_kronecker(outB, kro_outB, dAx + d);
    }
    else
    {
        zz_pX kro_inA, kro_inB, kro_outA, kro_outB;
        to_kronecker(kro_inA, in, dAx + d);
        to_kronecker(kro_inB, in, dBx + d);
        middle_product(kro_outA, kro_revA, kro_inA << dAx, dAy * (dAx + d + 1) + dAx, dBy * (dAx + d + 1) - 1);
        middle_product(kro_outB, kro_revB, kro_inB << dBx, dBy * (dBx + d + 1) + dBx, dAy * (dBx + d + 1) - 1);
        from_kronecker(outA, kro_outA, dAx + d);
        from_kronecker(outB, kro_outB, dBx + d);
    }

    out.SetLength(n);
    for (long i = 0; i < outA.length(); i++)
        out[i] = outA[i];
    for (long i = outA.length(); i < dBy; i++)
        out[i] = 0;
    for (long i = 0; i < outB.length(); i++)
        out[i + dBy] = outB[i];
    for (long i = outB.length(); i < dAy; i++)
        out[i + dBy] = 0;
}

/*------------------------------------------------------------*/
/* left multiplication, matrix version                        */
/*------------------------------------------------------------*/
void sylvester_lzz_pX::mul_left(Mat<zz_pX>& out, const Mat<zz_pX>& in) const
{
    if (&out == &in)
    {
        Mat<zz_pX> out2 = mul_left(in);
        out = out2;
        return;
    }

    if (in.NumCols() != n)
        Error("Bad dimensions in sylvester_lzz_pX matrix mul left");

    Vec<zz_pX> vec_in, vec_out;
    vec_in.SetLength(n);
    long p = in.NumRows();
    out.SetDims(p, n);

    for (long j = 0; j < p; j++)
    {
        for (long i = 0; i < n; i++)
            vec_in[i] = in[j][i];
        mul_left(vec_out, vec_in);
        for (long i = 0; i < n; i++)
            out[j][i] = vec_out[i];
    }
}

/*------------------------------------------------------------*/
/* turns M into a dense matrix                                */
/*------------------------------------------------------------*/
void sylvester_lzz_pX::to_dense(Mat<zz_pX>& Mdense) const
{
    Mdense.SetDims(n, n);
    for (long i = 0; i < n; i++)
        for (long j = 0; j < n; j++)
            Mdense[i][j] = 0;

    for (long i = 0; i < dBy; i++)
        for (long j = 0; j <= dAy; j++)
            Mdense[i+j][i] = a[j];
    for (long i = 0; i < dAy; i++)
        for (long j = 0; j <= dBy; j++)
            Mdense[i+j][i + dBy] = b[j];
}

/*------------------------------------------------------------*/
/* G, H such that Z1 M - M Z0 = G H^t                         */
/*------------------------------------------------------------*/
void sylvester_lzz_pX::phi_plus_generators(Mat<zz_pX>& G, Mat<zz_pX>& H) const
{
    G.SetDims(n, 2);
    H.SetDims(n, 2);

    for (long i = 0; i < n; i++)
    {
        H[i][0] = 0;
        H[i][1] = 0;
        G[i][0] = 0;
        G[i][1] = 0;
    }
    H[dBy - 1][0] = 1;
    H[dAy + dBy - 1][1] = 1;

    for (long i = 0; i < dAy; i++)
        G[i + dBy][0] += a[i];
    G[0][0] += a[dAy];
    for (long i = 0; i < dBy; i++)
    {
        G[i][0] -= b[i];
        G[i + dAy][1] += b[i];
    } 
    
    G[dBy % (dAy + dBy)][0] -= b[dBy];   // in case dAy = 0
    G[0][1] += b[dBy];
} 

/*------------------------------------------------------------*/
/* finds a sequence of degrees n0, n1, .. nk                  */
/* n0 <= 1, ni = {2*n{i-1}, 2*{n-1}-2}, nk >= n               */
/*------------------------------------------------------------*/
static vector<long> degrees(long n)
{
    vector<long> all_deg;

    while(n > 1)
    {
        if (n & 1)
            n++;
        all_deg.insert(all_deg.begin(), n);
        n >>= 1;
    }
    all_deg.insert(all_deg.begin(), n);
    return all_deg;
}

/*------------------------------------------------------------*/
/* Newton iteration for inverse                               */
/* assumes M(0) invertible; error otherwise                   */
/* return M^-1 mod x^m as a toeplitz_like_minus matrix        */
/*------------------------------------------------------------*/
void sylvester_lzz_pX::newton_inv_trunc(toeplitz_like_minus_lzz_pX& iM, long m) const
{
    zz_pX a0, b0;
    for (long i = 0; i < a.length(); i++)
        SetCoeff(a0, i, coeff(a[i], 0));
    for (long i = 0; i < b.length(); i++)
        SetCoeff(b0, i, coeff(b[i], 0));

    sylvester_lzz_p S(a0, b0);
    toeplitz_like_minus_lzz_p iS;
    long r = S.inv(iS);
    
    if (r == 0)
        Error("Sylvester matrix not invertible at zero");

    Mat<zz_pX> G, H; // generators of iM;
    G = conv(iS.G);
    H = conv(iS.H);

    vector<long> all_deg=degrees(m);
    long k = all_deg[0];
    long idx = 1;

    while (k < m) 
    {
        // Mat<zz_pX> y;
        // Mat<zz_pX> tr = trunc(a, 2*k);
        // middle_product(y, x, tr, k, k-1);
        // mul_trunc(y, y, x, k);
        // y <<= k;
        // x = x - y;
        k = 2 * k;
        if (k != all_deg[idx])
        {
            k = all_deg[idx];
            // trunc(x, x, k);
        }            
        idx++;
    }
    // trunc(x, x, m);
}



// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
