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


    cout << a << "     " << revA << endl;
    cout << b << "     " << revB << endl;

    long d = deg(in);
    Vec<zz_pX> outA, outB;
    zz_pX kro_revA, kro_revB;
    to_kronecker(kro_revA, revA, dAx + d);  
    to_kronecker(kro_revB, revB, dBx + d);  
    
    if (dAx == dBx)
    {
        zz_pX kro_in, kro_outA, kro_outB;
        to_kronecker(kro_in, in, dAx + d);
        cout << "kro_in" << kro_in << endl;
        cout << "kro_revA" << kro_revA << endl;
        cout << "kro_revB" << kro_revB << endl;
        middle_product(kro_outA, kro_revA, kro_in, dAy * (dAx + d), dBy * (dAx + d) - 1);
        middle_product(kro_outB, kro_revB, kro_in, dBy * (dAx + d), dAy * (dAx + d) - 1);
        cout << "kro_outA" << kro_outA << endl;
        cout << "kro_outB" << kro_outB << endl;
        from_kronecker(outA, kro_outA, dAx + d);
        from_kronecker(outB, kro_outB, dAx + d);
        cout << "outA" << outA << endl;
        cout << "outB" << outB << endl;
    }
    else
    {
        zz_pX kro_inA, kro_inB, kro_outA, kro_outB;
        to_kronecker(kro_inA, in, dAx + d);
        to_kronecker(kro_inB, in, dBx + d);
        middle_product(kro_outA, kro_revA, kro_inA, dAy * (dAx + d), dBy * (dAx + d) - 1);
        middle_product(kro_outB, kro_revB, kro_inB, dBy * (dBx + d), dAy * (dBx + d) - 1);
        from_kronecker(outA, kro_outA, dAx + d);
        from_kronecker(outB, kro_outB, dBx + d);
    }

    out.SetLength(n);
    for (long i = 0; i < outA.length(); i++)
        out[i] = outA[i];
    for (long i = outA.length(); i < dBy; i++)
        out[i] = 0;
    for (long i = 0; i < outB.length(); i++)
        out[i + dAy] = outB[i];
    for (long i = outB.length(); i < dBy; i++)
        out[i + dAy] = 0;
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

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
