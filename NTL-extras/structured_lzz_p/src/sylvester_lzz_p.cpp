#include <NTL/lzz_pX.h>
#include <NTL/mat_lzz_p.h>

#include "lzz_p_extra.h"
#include "lzz_pX_extra.h"
#include "mat_lzz_p_extra.h"
#include "lzz_pX_middle_product.h"
#include "structured_lzz_p.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/*                   SYLVESTER MATRICES                       */
/* = (square) matrix of (F,G) \mapsto AF + BG                 */
/* particular case of mosaic toeplitz                         */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* sets dimensions to 0                                       */
/*------------------------------------------------------------*/
sylvester_lzz_p::sylvester_lzz_p()
{
    n = 0;
}

/*------------------------------------------------------------*/
/* dimension = deg(A) + deg(B), A and B as above              */
/*------------------------------------------------------------*/
sylvester_lzz_p::sylvester_lzz_p(const zz_pX& A, const zz_pX& B)
{
    a = A;
    b = B;
    n = deg(A) + deg(B);
    revA.rep.SetLength(deg(A) + 1);
    for (long i = deg(A); i >= 0; i--)
        SetCoeff(revA, i, coeff(A, deg(A) - i));
    revB.rep.SetLength(deg(B) + 1);
    for (long i = deg(B); i >= 0; i--)
        SetCoeff(revB, i, coeff(B, deg(B) - i));

    invA = 0;
    invB = 0;
    singular = -1; // unknown
}

/*------------------------------------------------------------*/
/* getters                                                    */
/*------------------------------------------------------------*/
long sylvester_lzz_p::NumCols() const
{
    return n;
}

long sylvester_lzz_p::NumRows() const
{
    return n;
}


/*------------------------------------------------------------*/
/* right multiplication                                       */
/*------------------------------------------------------------*/
void sylvester_lzz_p::mul_right(Vec<zz_p>& out, const Vec<zz_p>& in) const
{
    if (&out == &in)
    {
        Vec<zz_p> out2 = mul_right(in);
        out = out2;
        return;
    }

    if (in.length() != n)
        Error("Bad input length in sylvester_lzz_p mul right");
    
    zz_pX F, G, H;
    for (long i = deg(b) - 1; i >= 0; i--)
        SetCoeff(F, i, in[i]);
    for (long i = deg(a) - 1; i >= 0; i--)
        SetCoeff(G, i, in[i + deg(b)]);
    H = F*a + G*b;
    out.SetLength(n);
    for (long i = 0; i < n; i++)
        out[i] = coeff(H, i);
}
 
/*------------------------------------------------------------*/
/* right matrix multiplication                                */
/*------------------------------------------------------------*/
void sylvester_lzz_p::mul_right(Mat<zz_p>& out, const Mat<zz_p>& in) const
{
    if (&out == &in)
    {
        Mat<zz_p> out2 = mul_right(in);
        out = out2;
        return;
    }

    if (in.NumRows() != n)
        Error("Bad dimensions in sylvester_lzz_p matrix mul right");

    Vec<zz_p> vec_in, vec_out;
    vec_in.SetLength(n);
    long p = in.NumCols();
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
void sylvester_lzz_p::mul_left(Vec<zz_p>& out, const Vec<zz_p>& in) const
{
    if (&out == &in)
    {
        Vec<zz_p> out2 = mul_left(in);
        out = out2;
        return;
    }

    if (in.length() != n)
        Error("Bad input length in sylvester_lzz_p mul right");

    out.SetLength(n);
    zz_pX F, poly_in;
    poly_in.SetLength(n);

    for (long i = n-1; i >= 0; i--)
        SetCoeff(poly_in, i, in[i]);
    middle_product(F, revA, poly_in, deg(a), deg(b)-1);
    for (long i = 0; i < deg(b); i++)
        out[i] = coeff(F, i);
    middle_product(F, revB, poly_in, deg(b), deg(a)-1);
    for (long i = 0; i < deg(a); i++)
        out[i + deg(b)] = coeff(F, i);
}
 
/*------------------------------------------------------------*/
/* left matrix multiplication                                 */
/*------------------------------------------------------------*/
void sylvester_lzz_p::mul_left(Mat<zz_p>& out, const Mat<zz_p>& in) const
{
    if (&out == &in)
    {
        Mat<zz_p> out2 = mul_left(in);
        out = out2;
        return;
    }

    if (in.NumCols() != n)
        Error("Bad dimensions in sylvester_lzz_p matrix mul left");

    Vec<zz_p> vec_in, vec_out;
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
void sylvester_lzz_p::to_dense(Mat<zz_p>& Mdense) const
{
    Mdense.SetDims(n, n);
    for (long i = 0; i < n; i++)
        for (long j = 0; j < n; j++)
            Mdense[i][j] = 0;

    for (long i = 0; i < deg(b); i++)
        for (long j = 0; j <= deg(a); j++)
            Mdense[i+j][i] = coeff(a, j);
    for (long i = 0; i < deg(a); i++)
        for (long j = 0; j <= deg(b); j++)
            Mdense[i+j][i + deg(b)] = coeff(b, j);
}

/*------------------------------------------------------------*/
/* precomputation for system solving                          */
/*------------------------------------------------------------*/
void sylvester_lzz_p::prepare_inverse()
{
    zz_pX d, s, t;
    XGCD(d, s, t, a, b);
    if (deg(d) == 0)
    {
        zz_p id = 1 / coeff(d, 0);
        invA = id * s;
        invB = id * t;
        build(modA, a);
        build(modB, b);
        build(multIA, invA, modB);
        build(multIB, invB, modA);

        invSerA = InvTrunc(revA, n);
        invSerB = InvTrunc(revB, n);

        singular = 0;
    }
    else
    {
        singular = 1;
    }
}

/*------------------------------------------------------------*/
/* right system solving                                       */
/* returns 0 if the matrix is non-invertible                  */
/* (even though there may be solutions)                       */
/* return 1 otherwise                                         */
/*------------------------------------------------------------*/
long sylvester_lzz_p::solve(Vec<zz_p>& out, const Vec<zz_p>& in) const
{
    if (&out == &in)
    {
        Vec<zz_p> out2;
        long r = solve(out2, in);
        out = out2;
        return r;
    }
    
    zz_pX F, FA, FB;
    for (long i = n-1; i >= 0; i--)
        SetCoeff(F, i, in[i]);
    FA = F % a;
    FB = F % b;
    
    if (singular == 1)
    {
        out.SetLength(0);
        return 0;
    }

    if (singular == 0) // inverse known
    {
        MulMod(FA, FA, multIB, modA);
        MulMod(FB, FB, multIA, modB);
    }
    else // not known 
    {
        zz_pX d, s, t;
        XGCD(d, s, t, a, b);
        if (deg(d) != 0)
        {
            out.SetLength(0);
            return 0;
        }
        zz_p id = 1 / coeff(d, 0);
        zz_pX iA = id * s;
        zz_pX iB = id * t;
        FA = MulMod(FA, iB, a);
        FB = MulMod(FB, iA, b);
    }
    
    out.SetLength(n);
    for (long i = 0; i < deg(b); i++)
        out[i] = coeff(FB, i);
    for (long i = 0; i < deg(a); i++)
        out[i + deg(b)] = coeff(FA, i);
    return 1;
}

/*------------------------------------------------------------*/
/* right system solving, matrix version                       */
/*------------------------------------------------------------*/
long sylvester_lzz_p::solve(Mat<zz_p>& out, const Mat<zz_p>& in) const
{
    if (&out == &in)
    {
        Mat<zz_p> out2;
        long r = solve(out2, in);
        out = out2;
        return r;
    }

    Vec<zz_p> vec_in, vec_out;
    vec_in.SetLength(n);
    long p = in.NumCols();
    out.SetDims(n, p);

    for (long j = 0; j < p; j++)
    {
        for (long i = 0; i < n; i++)
            vec_in[i] = in[i][j];
        long r = solve(vec_out, vec_in);
        if (r == 0)
        {
            out.SetDims(0, 0);
            return 0;
        }
        for (long i = 0; i < n; i++)
            out[i][j] = vec_out[i];
    }
    return 1;
}

/*------------------------------------------------------------*/
/* left system solving                                        */
/* returns 0 if the matrix is non-invertible                  */
/* (even though there may be solutions)                       */
/* return 1 otherwise                                         */
/*------------------------------------------------------------*/
long sylvester_lzz_p::solve_transpose(Vec<zz_p>& out, const Vec<zz_p>& in) const
{
    if (&out == &in)
    {
        Vec<zz_p> out2;
        long r = solve_transpose(out2, in);
        out = out2;
        return r;
    }

    if (singular == 1)
    {
        out.SetLength(0);
        return 0;
    }
    
    Vec<zz_p> vb, va;
    vb.SetLength(deg(b));
    va.SetLength(deg(a));

    for (long i = deg(b)-1; i >= 0; i--)
        vb[i] = in[i];

    for (long i = deg(a)-1; i >= 0; i--)
        va[i] = in[i + deg(b)];


    if (singular == 0) // inverse known
    {
        UpdateMap(va, va, multIB, modA);
        UpdateMap(vb, vb, multIA, modB);
    }
    else // not known 
    {
        zz_pX d, s, t;
        XGCD(d, s, t, a, b);
        if (deg(d) != 0)
        {
            out.SetLength(0);
            return 0;
        }
        zz_p id = 1 / coeff(d, 0);
        zz_pX iA = id * s;
        zz_pX iB = id * t;
        zz_pXModulus mA(a), mB(b);
        UpdateMap(va, va, zz_pXMultiplier(iB, mA), mA);
        UpdateMap(vb, vb, zz_pXMultiplier(iA, mB), mB);
    }

    zz_pX FA, FB; 
    for (long i = va.length()-1; i >= 0; i--)
        SetCoeff(FA, i, va[i]);
    for (long i = vb.length()-1; i >= 0; i--)
        SetCoeff(FB, i, vb[i]);

    zz_pX numerA = MulTrunc(FA, revA, deg(a));
    zz_pX numerB = MulTrunc(FB, revB, deg(b));

    if (singular == 0)
    {
        numerA = MulTrunc(numerA, invSerA, n);
        numerB = MulTrunc(numerB, invSerB, n);
    }
    else
    {
        zz_pX iSerA = InvTrunc(revA, n);
        zz_pX iSerB = InvTrunc(revB, n);
        numerA = MulTrunc(numerA, iSerA, n);
        numerB = MulTrunc(numerB, iSerB, n);
    }

    
    out.SetLength(n);
    for (long i = 0; i < n; i++)
        out[i] = coeff(numerA, i) + coeff(numerB, i);

    return 1;
}

/*------------------------------------------------------------*/
/* left system solving, matrix version                        */
/*------------------------------------------------------------*/
long sylvester_lzz_p::solve_transpose(Mat<zz_p>& out, const Mat<zz_p>& in) const
{
    if (&out == &in)
    {
        Mat<zz_p> out2;
        long r = solve_transpose(out2, in);
        out = out2;
        return r;
    }

    if (in.NumCols() != n)
        Error("Bad dimensions for sylvester_lzz_p solve transpose");

    Vec<zz_p> vec_in, vec_out;
    vec_in.SetLength(n);
    long p = in.NumRows();
    out.SetDims(p, n);

    for (long j = 0; j < p; j++)
    {
        for (long i = 0; i < n; i++)
            vec_in[i] = in[j][i];
        long r = solve_transpose(vec_out, vec_in);
        if (r == 0)
        {
            out.SetDims(0, 0);
            return 0;
        }
        for (long i = 0; i < n; i++)
            out[j][i] = vec_out[i];
    }
    return 1;
}

/*------------------------------------------------------------*/
/* G, H such that Z1 M - M Z0 = G H^t                         */
/*------------------------------------------------------------*/
void sylvester_lzz_p::phi_plus_generators(Mat<zz_p>& G, Mat<zz_p>& H) const
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
    H[deg(b) - 1][0] = 1;
    H[deg(a) + deg(b) - 1][1] = 1;

    for (long i = 0; i < deg(a); i++)
        G[i + deg(b)][0] += coeff(a, i);
    G[0][0] += coeff(a, deg(a));
    for (long i = 0; i < deg(b); i++)
    {
        G[i][0] -= coeff(b, i);
        G[i + deg(a)][1] += coeff(b, i);
    }
    
    G[deg(b) % (deg(a) + deg(b))][0] -= coeff(b, deg(b));   // in case deg(a) = 0
    G[0][1] += coeff(b, deg(b));
}

/*------------------------------------------------------------*/
/* finds the inverse of this as a topelitz_like_minus matrix  */
/* if return value is 0, singular matrix                      */
/*------------------------------------------------------------*/
long sylvester_lzz_p::inv(toeplitz_like_minus_lzz_p& inv)
{
    prepare_inverse();
    if (singular == 1)
        return 0;
    Mat<zz_p> U, V, iU, iV;
    phi_plus_generators(U, V);
    solve(iU, U);
    solve_transpose(iV, transpose(V));
    inv = toeplitz_like_minus_lzz_p(-iU, transpose(iV));
    return 1;
}

/*------------------------------------------------------------*/
/* compute the size-m top-right block of M^(-1)               */
/* r = 0 if matrix non invertible, r = 1 otherwise            */
/*------------------------------------------------------------*/
long sylvester_lzz_p::top_right_block_inverse(Mat<zz_p>& block, long m)
{
    if (m > NumCols())
        Error("Bad size for sylvester top right block inverse");

    zz_pX d, s, t;
    XGCD(d, s, t, a, b); 
    if (deg(d) != 0)
        return 0;

    zz_pX invA, invB;
    zz_p id = 1 / coeff(d, 0);

    if (id != 1)
    {
        s *= id;
        t *= id;
    }

    long dA = deg(a), dB = deg(b);
    Vec<zz_p> V0(INIT_SIZE, dA+dB);
    Vec<zz_p> V1(INIT_SIZE, dA+dB);
    Vec<zz_p> U0(INIT_SIZE, dA+dB);
    Vec<zz_p> U1(INIT_SIZE, dA+dB);

    zz_pX rS;
    rS.rep.SetLength(dB);
    for (long i = 0; i < dB; i++)
        rS[i] = s[dB - 1 - i];
    rS.normalize();

    zz_pX rT;
    rT.rep.SetLength(dA);
    for (long i = 0; i < dA; i++)
        rT[i] = t[dA - 1 - i];
    rT.normalize();

    zz_pX rS_over_rB, R_rS_rB, rT_over_rA, R_rT_rA;
    zz_pX top, bottom;  

    InvTruncMul(rS_over_rB, rS, revB, n);   
    for (long i = 0; i <= deg(rS_over_rB); i++)
        V0[i] = rS_over_rB[i];
    for (long i = deg(rS_over_rB)+1; i < n; i++)
        V0[i] = 0;
    
    R_rS_rB.rep.SetLength(n);
    for (long i = 0; i < max(m, dB); i++)
        R_rS_rB[i] = coeff(rS_over_rB, n - 1 - i);
    R_rS_rB.normalize();                    
    MulTrunc(top, R_rS_rB, trunc(b, max(m, dB)), max(m, dB)); 

    InvTruncMul(rT_over_rA, rT, revA, n);
    for (long i = 0; i <= deg(rT_over_rA); i++)
        V1[i] = rT_over_rA[i];
    for (long i = deg(rT_over_rA)+1; i < n; i++)
        V1[i] = 0;

    if (m > dB)  // if m <= dB, don't need bottom
    {
        R_rT_rA.rep.SetLength(n);
        for (long i = 0; i < dA; i++)
            R_rT_rA[i] = coeff(rT_over_rA, n - 1 - i);
        R_rT_rA.normalize();
        MulTrunc(bottom, R_rT_rA, a, dA+1);
    }

    zz_p lC = coeff(a, dA) * coeff(b, dB);
    for (long i = 0; i < max(m, dB); i++)
    {
        U0[i] = U1[i] = lC * (coeff(s, i) + coeff(top, i));
        U0[i] -= coeff(b, i);
    }
    if (m > dB)
        for (long i = 0; i < dA; i++)
        {
            U0[i + dB] = U1[i + dB] = lC * (coeff(t, i) + coeff(bottom, i));
            U1[i + dB] -= coeff(a, i);
        }
    U0[dB] -= coeff(b, dB);

    Mat<zz_p> L0, R0, L1, R1;
    L0.SetDims(m, m);
    R0.SetDims(m, m);
    L1.SetDims(m, m);
    R1.SetDims(m, m);
    for (long j = 0; j < m; j++)
        for (long i = j; i < m; i++)
        {
            L0[i][j] = U0[i-j];
            L1[i][j] = U1[i-j];
        }
    
    if (n - 2*m >= 0)
        for (long i = 0; i < m; i++)
            for (long j = 0, idx = n-m-1-i; j < m; j++, idx++)
            {
                R0[i][j] = V0[idx];
                R1[i][j] = V1[idx];
            }
    else
        for (long i = 0; i < m; i++)
            for (long j = 0; j < m; j++)
            {
                long idx = n-m-1+j-i;
                while (idx < 0)
                    idx += n;
                R0[i][j] = V0[idx];
                R1[i][j] = V1[idx];
            }
        

    block = L0*R0 + L1*R1;
    return 1;
}


// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
