#include <NTL/lzz_pX.h>
#include <NTL/matrix.h>
#include <NTL/vector.h>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <numeric>
#include <random>
#include <NTL/BasicThreadPool.h>
#include <cmath>

//#define SAFETY_CHECKS

#include "util.h"
#include "mat_lzz_pX_extra.h"
#include "mat_lzz_pX_approximant.h"
#include "sage_output.h"

NTL_CLIENT

std::ostream &operator<<(std::ostream &out, const VecLong &s)
{
    out << "[ ";
    for (auto &i: s)
        out << i << " ";
    return out << "]";
}

typedef Vec<Vec<Vec<zz_p>>> Coeffs;

static inline void MulMod_simple(zz_pX& upper, const zz_pX& a, const zz_pXModulus& F)
{
    long n = F.n;
    zz_pX P1(INIT_SIZE, n);
    fftRep R1(INIT_SIZE, F.l);
    TofftRep_trunc(R1, a, F.l, max(1L << F.k, 2*n-2));
    mul(R1, R1, F.HRep);
    FromfftRep(upper, R1, n-1, 2*n-3);
}


#if (1)
void MulMod_local(zz_pX& x, zz_pX& upper, const zz_pX& a, const zz_pXMultiplier& B, const zz_pXModulus& F)
{
    
    long n = F.n;
    long da;
    
    da = deg(a);
    
    if (da >= n)
        LogicError(" bad args to MulMod_local(zz_pX, zz_pX, zz_pXMultiplier, zz_pXModulus)");
    
    if (da < 0)
    {
        clear(x);
        return;
    }
    
    if (!B.UseFFT || !F.UseFFT || da <= NTL_zz_pX_MOD_CROSSOVER)
    {
        zz_pX P1;
        mul(P1, a, B.b);
        rem(x, P1, F);
        // todo: do this case
        upper = 0;
        return;
    }
    
    zz_pX P1(INIT_SIZE, n), P2(INIT_SIZE, n);
    fftRep R1(INIT_SIZE, F.l), R2(INIT_SIZE, F.l);
    
    long len;
    if (zz_p::IsFFTPrime())
        len = n;
    else
        len = 1L << F.k;
    
    TofftRep_trunc(R1, a, F.l, max(1L << F.k, 2*n-2));
    mul(R2, R1, B.B1);
    FromfftRep(P1, R2, n-1, 2*n-3);
    
    mul(R2, R1, F.HRep);
    FromfftRep(upper, R2, n-1, 2*n-3);
    
    reduce(R1, R1, F.k);
    mul(R1, R1, B.B2);
    TofftRep_trunc(R2, P1, F.k, len);
    mul(R2, R2, F.FRep);
    sub(R1, R1, R2);
    
    FromfftRep(x, R1, 0, n-1);
}

#else

void MulMod_local(zz_pX& x, zz_pX & upper, const zz_pX& a, const zz_pXMultiplier& B, const zz_pXModulus& F)
{
    
    long n = F.n;
    long da;
    
    da = deg(a);
    
    if (da >= n)
        LogicError(" bad args to MulMod_local(zz_pX, zz_pX, zz_pXMultiplier, zz_pXModulus)");
    
    if (da < 0)
    {
        clear(x);
        return;
    }
    
    if (!B.UseFFT || !F.UseFFT || da <= NTL_zz_pX_MOD_CROSSOVER)
    {
        zz_pX P1;
        mul(P1, a, B.b);
        rem(x, P1, F);
        // todo: do this case
        upper = 0;
        return;
    }
    zz_pX P1(INIT_SIZE, n), P2(INIT_SIZE, n);
    fftRep R1(INIT_SIZE, F.l), R2(INIT_SIZE, F.l);
    
    
    TofftRep(R1, a, F.l); // l = next power of 2n-3, n=deg(modulus)
    mul(R2, R1, B.B1);    // B.B1 =  b *rev(1/rev(f))
    FromfftRep(P1, R2, n-1, 2*n-3);
    
    mul(R2, R1, F.HRep);
    FromfftRep(upper, R2, n-1, 2*n-3);
    
    reduce(R1, R1, F.k);
    mul(R1, R1, B.B2);
    TofftRep(R2, P1, F.k);
    mul(R2, R2, F.FRep);
    sub(R1, R1, R2);
    
    FromfftRep(x, R1, 0, n-1);
}

#endif

/*------------------------------------------------------------*/
/* generates (t*a^(s*i) mod g) for 0 <= i < l                 */
/* if need_upper_product = 1, also returns                    */
/*   S * rev(t * a^(s*i) mod g) mod x^(n-1), 0 <= i < l       */
/*------------------------------------------------------------*/
void gen_pows (Vec<zz_pX> &pow, Vec<zz_pX>&upper,
               const zz_pX &t,
               const zz_pX &a,
               const long s,
               const zz_pX &g,
               const long l, const long need_upper_products = 0)
{
    zz_pXModulus g_mod;
    build(g_mod, g);
    
    zz_pX a_pow;
    PowerMod(a_pow, a, s, g_mod); // a_pow = a^s
    
    zz_pXMultiplier a_pow_mul;
    build(a_pow_mul, a_pow, g_mod);
    
    pow.SetLength(l);
    pow[0] = t;
    
    
    if (need_upper_products)
    {
        upper.SetLength(l);
        long n = deg(g);
        for (long i = 1; i < l; i++)
        {
            MulMod_local(pow[i], upper[i - 1], pow[i - 1], a_pow_mul, g_mod);
            reverse(upper[i - 1], upper[i - 1], n - 2);
        }
        zz_pX tmp;
        MulMod_local(tmp, upper[l - 1], pow[l - 1], a_pow_mul, g_mod);
        reverse(upper[l - 1], upper[l - 1], n - 2);
    }
    else
    {
        for (long i = 1; i < l; i++)
        {
            MulMod(pow[i], pow[i - 1], a_pow_mul, g_mod);
        }
    }
    
}

/*------------------------------------------------------------*/
/* matrix multiplication using the algorithm of Giorgi et al. */
/* uses matrix multiplication for evaluation and interpolation*/
/*------------------------------------------------------------*/
void mul(Mat<zz_pX>& c, const Mat<zz_pX>& a, const Mat<zz_pX>& b)
{
    long dA = deg(a);
    long dB = deg(b);
    long sA = dA + 1;
    long sB = dB + 1;
    long m = a.NumRows();
    long n = a.NumCols();
    long p = b.NumCols();
    long ell;
    
    Mat<zz_p> tmp_mat(INIT_SIZE, sA, m * n);
    Mat<zz_p> valA, valB, valC;
    Mat<zz_p> vA, vB, iv;
    Mat<zz_p> valAp, valBp, valCp;
    
    vandermonde(vA, vB, iv, dA, dB);
    long nb_points = vA.NumRows();
    
    ell = 0;
    for (long i = 0; i < m; i++)
        for (long j = 0; j < n; j++)
        {
            long d = deg(a[i][j]);
            if (d >= 0)
            {
                const zz_p * cAij = a[i][j].rep.elts();
                long k;
                for (k = 0; k <= d; k++)  // k <= d-2 so k+1 <= d-1
                {
                    tmp_mat[k][ell] = cAij[k];
                }
            }
            for (long k = d+1; k <= dA; k++)
            {
                tmp_mat[k][ell] = 0;
            }
            
            ell++;
        }
    valA = vA * tmp_mat;
    
    tmp_mat.SetDims(sB, n * p);
    ell = 0;
    for (long i = 0; i < n; i++)
        for (long j = 0; j < p; j++)
        {
            long d = deg(b[i][j]);
            if (d >= 0)
            {
                const zz_p * cBij = b[i][j].rep.elts();
                long k;
                for (k = 0; k <= d; k++)  // k <= d-2 so k+1 <= d-1
                {
                    tmp_mat[k][ell] = cBij[k];
                }
            }
            for (long k = d+1; k <= dA; k++)
            {
                tmp_mat[k][ell] = 0;
            }
            ell++;
        }
    valB = vB * tmp_mat;
    
    valAp.SetDims(m, n);
    valBp.SetDims(n, p);
    valC.SetDims(nb_points, m * p);
    for (long i = 0; i < nb_points; i++)
    {
        long ell;
        ell = 0;
        for (long u = 0; u < m; u++)
            for (long v = 0; v < n; v++)
            {
                valAp[u][v] = valA[i][ell++];
            }
        
        ell = 0;
        for (long u = 0; u < n; u++)
            for (long v = 0; v < p; v++)
            {
                valBp[u][v] = valB[i][ell++];
            }
        
        valCp = valAp * valBp;
        
        ell = 0;
        for (long u = 0; u < m; u++)
            for (long v = 0; v < p; v++)
            {
                valC[i][ell++] = valCp[u][v];
            }
    }
    tmp_mat = iv*valC;
    
    c.SetDims(m, p);
    ell = 0;
    for (long u = 0; u < m; u++)
        for (long v = 0; v < p; v++)
        {
            c[u][v].rep.SetLength(nb_points);
            zz_p * cc = c[u][v].rep.elts();
            for (long i = 0; i < nb_points; i++)
                cc[i] = tmp_mat[i][ell];
            c[u][v].normalize();
            ell++;
        }
}


/*------------------------------------------------------------*/
/* matrix multiplication using the algorithm of Giorgi et al. */
/* uses matrix multiplication for evaluation and interpolation*/
/*                                                            */
/* variant does one product and a row-shifted product at once */
/*------------------------------------------------------------*/
void mul_special(Mat<zz_pX>& c, Mat<zz_pX>& cs, const Mat<zz_pX>& a, const Mat<zz_pX>& b)
{
    long dA = deg(a);
    long dB = deg(b);
    long sA = dA + 1;
    long sB = dB + 1;
    long m = a.NumRows();
    long n = a.NumCols();
    long p = b.NumCols();
    long ell;
    
    Mat<zz_p> tmp_mat(INIT_SIZE, sA, m * n);
    Mat<zz_p> valA, valB, valC, valCs;
    Mat<zz_p> vA, vB, iv;
    Mat<zz_p> valAp, valBp, valBsp, valCp, valCsp;
    
    vandermonde(vA, vB, iv, dA, dB);
    long nb_points = vA.NumRows();
    
    ell = 0;
    for (long i = 0; i < m; i++)
        for (long j = 0; j < n; j++)
        {
            long d = deg(a[i][j]);
            if (d >= 0)
            {
                const zz_p * cAij = a[i][j].rep.elts();
                for (long k = 0; k <= d; k++)  // k <= d-2 so k+1 <= d-1
                {
                    tmp_mat[k][ell] = cAij[k];
                }
            }
            for (long k = d+1; k <= dA; k++)
            {
                tmp_mat[k][ell] = 0;
            }
            ell++;
        }
    valA = vA * tmp_mat;
    
    tmp_mat.SetDims(sB, n * p);
    ell = 0;
    for (long i = 0; i < n; i++)
        for (long j = 0; j < p; j++)
        {
            long d = deg(b[i][j]);
            if (d >= 0)
            {
                const zz_p * cBij = b[i][j].rep.elts();
                for (long k = 0; k <= d; k++)  // k <= d-2 so k+1 <= d-1
                {
                    tmp_mat[k][ell] = cBij[k];
                }
            }
            for (long k = d+1; k <= dB; k++)
            {
                tmp_mat[k][ell] = 0;
            }
            ell++;
        }
    valB = vB * tmp_mat;
    
    valAp.SetDims(m, n);
    valBp.SetDims(n, p);
    valBsp.SetDims(n, p);
    valC.SetDims(nb_points, m * p);
    valCs.SetDims(nb_points, m * p);
    for (long i = 0; i < nb_points; i++)
    {
        long ell;
        ell = 0;
        for (long u = 0; u < m; u++)
            for (long v = 0; v < n; v++)
            {
                valAp[u][v] = valA[i][ell++];
            }
        
        ell = 0;
        for (long u = 0; u < n-1; u++)
            for (long v = 0; v < p; v++)
            {
                valBsp[u][v] = valB[i][ell + p];
                valBp[u][v] = valB[i][ell++];
            }
        for (long v = 0; v < p; v++)
        {
            valBsp[n-1][v] = 0;
            valBp[n-1][v] = valB[i][ell++];
        }
        
        valCp = valAp * valBp;
        ell = 0;
        for (long u = 0; u < m; u++)
            for (long v = 0; v < p; v++)
            {
                valC[i][ell++] = valCp[u][v];
            }
        
        valCsp = valAp * valBsp;
        ell = 0;
        for (long u = 0; u < m; u++)
            for (long v = 0; v < p; v++)
            {
                valCs[i][ell++] = valCsp[u][v];
            }
    }
    tmp_mat = iv*valC;
    c.SetDims(m, p);
    ell = 0;
    for (long u = 0; u < m; u++)
        for (long v = 0; v < p; v++)
        {
            c[u][v].rep.SetLength(nb_points);
            zz_p * cc = c[u][v].rep.elts();
            for (long i = 0; i < nb_points; i++)
                cc[i] = tmp_mat[i][ell];
            c[u][v].normalize();
            ell++;
        }
    
    tmp_mat = iv*valCs;
    cs.SetDims(m, p);
    ell = 0;
    for (long u = 0; u < m; u++)
        for (long v = 0; v < p; v++)
        {
            cs[u][v].rep.SetLength(nb_points);
            zz_p * cc = cs[u][v].rep.elts();
            for (long i = 0; i < nb_points; i++)
                cc[i] = tmp_mat[i][ell];
            cs[u][v].normalize();
            ell++;
        }
}

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
void get_quos (Mat<zz_pX> &quos,
               const Vec<zz_pX> &alphas,
               const Vec<zz_pX> &As,
               const Vec<zz_pX> &upper,
               const zz_pX &g,
               const long m)
{
    cout << "-enter quos\n";
    double t;
    
    t = GetWallTime();
    const long n = deg(g);
    const long pad = 2*m - ((n-2*m-1) % (2*m));
    const long len = ceil((n*1.0) / (2*m));
    Mat<zz_pX> mat1, mat2;
    mat1.SetDims(alphas.length(), len);
    mat2.SetDims(len, alphas.length());
    cout << "setup: " << GetWallTime()-t << endl;
    
    t = GetWallTime();
    for (long i = 0; i < alphas.length(); i++)
    {
        zz_pX alpha_rev = reverse(alphas[i], n-1);
        long at;
        
        at = -pad;
        for (long j = 0; j < len; j++)
        {
            mat1[i][j].rep.SetLength(2*m);
            zz_p * c = mat1[i][j].rep.elts();
            for (long k = 0; k <= 2*m-1; k++)
                c[k] = coeff(alpha_rev, at++);
            mat1[i][j].normalize();
        }
        
        // not sure we need to make this so complicated.
        at = 0;
        const zz_p * s = upper[i].rep.elts();
        for (long j = 0; j < len; j++)
        {
            mat2[len - j - 1][i].rep.SetLength(2*m);
            zz_p * c = mat2[len - j - 1][i].rep.elts();
            long dg = min(2*m - 1, deg(upper[i]) - at);
            for (long k = 0; k <= dg; k++)
                c[k] = s[k];
            s += dg + 1;
            at += dg + 1;
            mat2[len - j - 1][i].normalize();
        }
    }
    cout << "build matrices: " << GetWallTime()-t << endl;
    
    t = GetWallTime();
    Mat<zz_pX> res1, res2;
    mul_special(res1, res2, mat1, mat2);
    cout << "mat mul: " << GetWallTime()-t << endl;
    
    t = GetWallTime();
    trunc(res1, res1, 2*m);
    RightShift(res2, res2, 2*m);
    quos = res1 + res2;
    reverse(quos, quos, 2 * m - 1);
    cout << "clean up: " << GetWallTime()-t << endl;
}

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
void SetDims (Coeffs &res, const long r, const long c, const long d)
{
    res.SetLength(r);
    for (long i = 0; i < r; i++)
    {
        res[i].SetLength(c);
        for (long j = 0; j < c; j++)
            res[i][j].SetLength(d);
    }
}

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
void format (Vec<Mat<zz_p>> &res, const Vec<Coeffs> &coeffs, const long d, const long m)
{
    res.SetLength(d);
    const long sqrt_d = ceil(sqrt(d));
    for (long i = 0; i < d; i++)
        res[i].SetDims(m, m);
    for (long i = 0; i < m; i++)
    {
        for (long u = 0; u < sqrt_d; u++)
        {
            for (long v = 0; v < sqrt_d; v++)
            {
                long pow = sqrt_d * v + u;
                for (long j = 0; j < m; j++)
                {
                    if (pow < d)
                        res[d-pow-1][m-i-1][j] = coeffs[i][u][v][j];
                }
            }
        }
    }
}

/*
void get_first_row (Coeffs &res,
                    const zz_pX &a,
                    const zz_pX &g,
                    const long m)
{
    const long n = deg(g);
    const long d = 2*ceil((n*1.0)/m)+1;
    const long sqrt_d = ceil(sqrt(d));
    
    Vec<zz_pX> alphas, As;
    zz_pX t;
    SetCoeff(t,m-1,1);
    
    double t0, t1;
    
    t1 = GetWallTime();
    t0 = GetWallTime();
    Vec<zz_pX> upper;
    cout << "-enter pows\n";
    gen_pows(alphas, upper, t, a, 1, g, sqrt_d, 0);
    cout << "pow1: " << GetWallTime() - t1 << endl;
    t1 = GetWallTime();
    gen_pows(As, upper, zz_pX(1), a, sqrt_d, g, sqrt_d, 1);
    cout << "pow2: " << GetWallTime() - t1 << endl;
    cout << "-total pow: " << GetWallTime() - t0 << endl;
    
    t1 = GetWallTime();
    Mat<zz_pX> quos;
    get_quos(quos, alphas, As, upper, g, m);
    cout << "-total quos: " << GetWallTime() - t1 << endl;
    
    t1 = GetWallTime();
    SetDims(res, sqrt_d, sqrt_d, 2*m);
    zz_pX g_trunc;
    trunc(g_trunc,g,2*m);
    for (long u = 0; u < sqrt_d; u++)
        for (long v = 0; v < sqrt_d; v++)
        {
            zz_pX rem1, rem2;
            auto alpha_trunc = trunc(alphas[u],2*m);
            auto A_trunc = trunc(As[v],2*m);
            auto quo_trunc = trunc(quos[u][v],2*m);
            MulTrunc(rem1, alpha_trunc, A_trunc, 2*m);
            MulTrunc(rem2, quo_trunc, g_trunc, 2*m);
            rem1 = rem1 - rem2;
            for (long s = 0; s < 2*m; s++)
                res[u][v][s] = coeff(rem1,s);
        }
    cout << "-rest from row1: " << GetWallTime() - t1 << endl;
}
*/

void get_first_row (Coeffs &res,
                    const zz_pX &t,
                    const zz_pX &a,
                    const zz_pX &g,
                    const long m)
{
    const long n = deg(g);
    const long d = 2*ceil((n*1.0)/m)+1;
    const long sqrt_d = ceil(sqrt(d));
    
    Vec<zz_pX> alphas, As;
    
    double t0, t1;
    
    t1 = GetWallTime();
    t0 = GetWallTime();
    Vec<zz_pX> upper;
    cout << "-enter pows\n";
    gen_pows(alphas, upper, t, a, 1, g, sqrt_d, 0);
    cout << "pow1: " << GetWallTime() - t1 << endl;
    t1 = GetWallTime();
    gen_pows(As, upper, zz_pX(1), a, sqrt_d, g, sqrt_d, 1);
    cout << "pow2: " << GetWallTime() - t1 << endl;
    cout << "-total pow: " << GetWallTime() - t0 << endl;

    t1 = GetWallTime();
    Mat<zz_pX> quos;
    get_quos(quos, alphas, As, upper, g, m);
    cout << "-total quos: " << GetWallTime() - t1 << endl;

    t1 = GetWallTime();
    SetDims(res, sqrt_d, sqrt_d, 2*m);
    zz_pX g_trunc;
    trunc(g_trunc,g,2*m);
    for (long u = 0; u < sqrt_d; u++)
        for (long v = 0; v < sqrt_d; v++)
        {
            zz_pX rem1, rem2;
            auto alpha_trunc = trunc(alphas[u],2*m);
            auto A_trunc = trunc(As[v],2*m);
            auto quo_trunc = trunc(quos[u][v],2*m);
            MulTrunc(rem1, alpha_trunc, A_trunc, 2*m);
            MulTrunc(rem2, quo_trunc, g_trunc, 2*m);
            rem1 = rem1 - rem2;
            for (long s = 0; s < 2*m; s++)
                res[u][v][s] = coeff(rem1,s);
        }
    cout << "-rest from row1: " << GetWallTime() - t1 << endl;
}

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
void get_coeffs (Vec<Mat<zz_p>> &mats, Vec<Coeffs> &res,
                 const zz_pX &a,
                 const zz_pX &g,
                 const long m)
{
    const long n = deg(g);
    const long d = 2*ceil((n*1.0)/m)+1;
    const long sqrt_d = ceil(sqrt(d));
    
    double t;
    
    t = GetWallTime();
    res.SetLength(m);
		zz_pX x_poly;
		SetCoeff(x_poly, m-1, zz_p(1));
    get_first_row(res[0],x_poly,a,g,m);
    cout << "--all first row: " << GetWallTime() - t << endl;
    
    t = GetWallTime();
    const zz_p inv_c0 = -(zz_p(1)/ConstTerm(g));
    for (long i = 1; i < m; i++)
    {
        SetDims(res[i], sqrt_d, sqrt_d, 2 * m - i);
        for (long u = 0; u < sqrt_d; u++)
            for (long v = 0; v < sqrt_d; v++)
            {
                zz_p l_n_1 =  inv_c0 * res[i-1][u][v][0];
                for (long j = 0; j < 2*m-i; j++)
                {
                    res[i][u][v][j] = res[i-1][u][v][j+1] + coeff(g, j+1) * l_n_1;
                }
            }
    }
    cout << "fill all: " << GetWallTime() - t << endl;
    
    // format
    t = GetWallTime();
    for (long i = 0; i < m; i++)
        SetDims(res[i], sqrt_d, sqrt_d, m);
    format(mats, res, d, m);
    cout << "format: " << GetWallTime() - t << endl;
}

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
void get_coeffs_naive (Vec<Coeffs> &res,
                       const zz_pX &a,
                       const zz_pX &g,
                       const long m)
{
    const long n = deg(g);
    const long d = 2*ceil((n*1.0)/m)+1;
    const long sqrt_d = ceil(sqrt(d));
    
    res.SetLength(m);
    
    zz_pXModulus g_mod;
    build(g_mod, g);
    for (long i = 0; i < m; i++)
    {
        SetDims(res[i], sqrt_d, sqrt_d, m);
        for (long u = 0; u < sqrt_d; u++)
            for (long v = 0; v < sqrt_d; v++)
            {
                zz_pX x;
                SetCoeff(x,m-i-1,1);
                zz_pX a_temp;
                PowerMod(a_temp, a, u+v*sqrt_d, g_mod);
                MulMod(a_temp, a_temp, x, g_mod);
                for (long t = 0; t < m; t++)
                    res[i][u][v][t] = coeff(a_temp,t);
            }
    }
}

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
void print(const Vec<Coeffs> &c)
{
    for (long i = 0; i < c.length(); i++)
        cout << c[i] << endl;
}

void print(const zz_pX &p)
{
    for (long i = 0; i <= deg(p); i++)
    {
        cout << coeff(p,i) << "*y^" <<i;
        if (i != deg(p)) cout << "+";
        else cout << endl;
    }
}


int main(int argc, char *argv[])
{
    SetSeed(ZZ(0));
    if (argc!=5)
        throw std::invalid_argument("Usage: ./test_charpoly_amodg degree expnt nbits nthreads");

    // main parameter: degree
    long n = atoi(argv[1]);
    // parameter: exponent for m = n^exponent; default = 0.333 ~ 1/omega for omega=3
    double exponent = (atof(argv[2]) <= 0.) ? 0.333 : (atof(argv[2]));

    // size of prime field
    long nbits = atoi(argv[3]);
    // whether we should verify the result
    SetNumThreads(atoi(argv[4]));

    // used for timings
    double t1,t2;
    double t_charpoly=0.0;

    t1 = GetWallTime();
    if (nbits==0)
        zz_p::FFTInit(0);
    else
        zz_p::init(NTL::GenPrime_long(nbits));

    {
        t_charpoly=0.0;
        // modulus
        zz_pX g;
        while (deg(g)<n)
            random(g, n+1);
        MakeMonic(g);

        // evaluation point
        zz_pX a;
        random(a, n); // deg < n: assumes a reduced modulo g
        
        zz_pX b;
        random(b, n);

        // parameters
        long m = pow(n,exponent);
        long d = ceil( (double)n / m );
        t2 = GetWallTime();

        std::cout << "Testing charpoly with random input polynomials" << std::endl;
        std::cout << "-- n = " << n << std::endl;
        std::cout << "-- m = " << m << std::endl;
        std::cout << "-- d = " << d << std::endl;
        std::cout << "--prime = " << zz_p::modulus() << std::endl;
        std::cout << "--nthreads = " << AvailableThreads() << std::endl;
        if ((nbits==0 && n<30) || (nbits != 0 && n*nbits < 1000))
        {
            std::cout << "--modulus g = " << g << std::endl;
            std::cout << "--point a = " << a << std::endl;
            std::cout << "--point h = " << b << std::endl;
        }
        else
        {
            std::cout << "--modulus g = <degree " << n << " polynomial>" << std::endl;
            std::cout << "--point a = <degree " << n << " polynomial>" << std::endl;
            std::cout << "--point h = <degree " << n << " polynomial>" << std::endl;
        }
        
        

#ifdef SAFETY_CHECKS
        std::cout << "###~~~Warning~~~### SAFETY_CHECKS is on; use with non-small degrees may be very slow." << std::endl;
#endif // SAFETY_CHECKS

        std::cout << "TIME ~~ build polynomials and compute parameters: " << (t2-t1) << std::endl;

        t1 = GetWallTime();
        zz_pXModulus G(g); // precompute things to work mod g
        zz_pXMultiplier A(a, G); // precompute things to multiply by a mod g
        t2 = GetWallTime();
        std::cout << "TIME ~~ pre-computations to work with a mod g: " << (t2-t1) << std::endl;
        t_charpoly += t2-t1;
        t1 = GetWallTime();

        // will hold powers of a, initially it is a
        zz_pX pow_a = a; 

#ifdef SAFETY_CHECKS
        // compute the block-Wiedemann sequence (naive)
        Vec<Mat<zz_p>> seq_naive;
        seq_naive.SetLength(2*d+1);
        ident(seq_naive[2*d], m); // pmat[2*d] = identity m x m
        for (long k = 2*d-1; k >= 0; --k)
        {
            seq_naive[k].SetDims(m, m);
            zz_pX f = pow_a;  // a^{2d-k} mod g
            for (long i = 0; i < m; ++i)
            {
                // here f = x^i * f mod g = x^i a^{2d-k} mod g
                for (long j = 0; j < m; ++j)
                    seq_naive[k][i][j] = coeff(f, j);  
                // note: here we have tranposed, because pmbasis works row-wise
                // f = x * f mod g  (recall g is monic)
                f <<= 1;
                f = f - coeff(f,n) * g;
            }
            MulMod(pow_a, pow_a, A, G); // a^{2d-(k+1)} mod g
        }
        t2 = GetWallTime();
        std::cout << "TIME ~~ compute sequence for block Wiedemann (naive): " << (t2-t1) << std::endl;
        pow_a = a; // set back to its state before this #ifdef
        t1 = GetWallTime();
#endif // SAFETY_CHECKS

        Vec<Mat<zz_p>> seq;
        Vec<Coeffs> res;
        get_coeffs(seq,res,a,g,m);
        t2 = GetWallTime();
       
#ifdef SAFETY_CHECKS
        bool correct=true;
        for (long k = 0; k < 2*d+1; ++k)
            if (seq[k] != seq_naive[k])
                correct=false;
#endif // SAFETY_CHECKS
        std::cout << "TIME ~~ compute sequence for block Wiedemann: " << (t2-t1);
#ifdef SAFETY_CHECKS
        std::cout << (correct ? " (correct)" : " (wrong)");
#endif // SAFETY_CHECKS
        std::cout << std::endl;
        t_charpoly += t2-t1;
        
        t1 = GetWallTime();
        // make Sb
        Vec<Mat<zz_p>> Sb;
        Sb.SetLength(2*d + 1);
        
				Coeffs coeffs;
        get_first_row(coeffs,b,a,g,m);

				long row = 0;
				long col = 0;
        for (long i = 0; i < 2*d+1; i++)
        {
            auto &temp = coeffs[row][col];
            Sb[Sb.length() - 1 - i].SetDims(m,1);
            for (long j = 0; j < m; j++)
                Sb[Sb.length() - 1 - i][j][0] = temp[j];
						row++;
						if (row == coeffs.length())
						{
							row = 0;
							col++;
						}
        }
				t2 = GetWallTime();
        std::cout << "TIME ~~ compute rhs sequence: " << (t2-t1) << endl;
        t_charpoly += t2-t1;
        
				
        t1 = GetWallTime();

        // convert the sequence into a polynomial matrix
        Mat<zz_pX> pmat;
        conv(pmat, seq, 2*d+1);
        
        Mat<zz_pX> pmat2;
        conv(pmat2, Sb, 2*d+1);
        
        // Matrix fraction reconstruction: add below
        // 1) polynomial matrix for Sb
        // 2) identity
        Mat<zz_pX> P;
        P.SetDims(2*m+1, m);
        for (long i = 0; i < m; i++)
            P[0][i] = -pmat2[i][0];
        for (long i = 0; i < m; i++)
            for (long j = 0; j < m; j++)
                P[i+1][j] = pmat[i][j];
        for (long i = 0; i < m; ++i)
            P[i+m+1][i] = zz_p(-1);
         
        //cout << "P:  " << P << endl;

        t2 = GetWallTime();
        std::cout << "TIME ~~ convert sequence to polynomial matrix: " << (t2-t1) << std::endl;
        t_charpoly += t2-t1;
        t1 = GetWallTime();

        // reconstruct fraction
        Mat<zz_pX> appbas;
        VecLong shift(2*m+1, 0);
        shift[0] = 2*d+1;
        VecLong pivdeg = pmbasis(appbas, P, 2*d+1, shift);

        t2 = GetWallTime();
        std::cout << "TIME ~~ matrix fraction reconstruction: " << (t2-t1) << std::endl;
        t_charpoly += t2-t1;
        
        t1 = GetWallTime();
        Mat<zz_pX> sysmat;
        sysmat.SetDims(m+1, m-1);
        
        for (long i = 0; i < m-1; i++)
            sysmat[0][i] = appbas[0][i+2];
        for (long i = 0; i < m; i++)
            for (long j = 0; j < m-1; j++)
                sysmat[i+1][j] = appbas[i+1][j+2];
                
        //cout << "sysmat: " << degree_matrix(sysmat) << endl;
        
        /*
				Mat<zz_pX> kerbas1;
        shift = VecLong(m+1, d);
        shift[0] += m*d;
        kernel_basis_zls_via_interpolation(kerbas1, sysmat, shift);
        cout << "shift: " << shift << endl;
        cout << "kerbas: " << degree_matrix(kerbas1) << endl;
        Mat<zz_pX> test_temp;
        multiply(test_temp, kerbas1,sysmat);
        cout << "test_temp: " << degree_matrix(test_temp) << endl;
        */

        Mat<zz_pX> kerbas;
        shift = VecLong(m+1, 0);
        shift[0] += m*d;
        pmbasis(kerbas, sysmat, m*d, shift);
       // cout << "kernel: " << degree_matrix(kerbas) << endl;
        
        zz_pX h = appbas[0][1];
        for (long i = 0; i < m; i++)
            h += kerbas[0][i+1] * appbas[i+1][1];
        
        t2 = GetWallTime();
        std::cout << "TIME ~~ compute kerbas and output: " << (t2-t1)<<endl;
        t_charpoly += t2-t1;
        
        zz_pXModulus g_mod;
        build(g_mod,g);
        zz_pX modcomp;
        CompMod(modcomp, h, a, g_mod);
        std::cout << "TIME ~~ compute h st h(a) mod g = b: " << t_charpoly;
        if (modcomp == b) cout << " (correct)" << endl;
        else cout << " (wrong)" << endl;
    }
    return 0;
}


// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
