#include "mat_lzz_pX_utils.h"
#include "mat_lzz_pX_arith.h"
#include "mat_lzz_pX_approximant.h"
#include "mat_lzz_pX_interpolant.h"
#include "mat_lzz_pX_sequence.h"

/*------------------------------------------------------------*/
/* TODO: what is this?                                        */
/*------------------------------------------------------------*/
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


/*------------------------------------------------------------*/
/* TODO: what is this?                                        */
/*------------------------------------------------------------*/
void gen_pows (
               Vec<zz_pX> &pow,
               Vec<zz_pX>&upper,
               const zz_pX &t,
               const zz_pX &a,
               const long s,
               const zz_pX &g,
               const long l,
               const bool need_upper_products
              )
{
    // build modulus for g
    zz_pXModulus g_mod(g);

    // build multiplier for a^s
    zz_pX a_pow;
    PowerMod(a_pow, a, s, g_mod); // a_pow = a^s
    zz_pXMultiplier a_pow_mul(a_pow, g_mod);

    // initialize list of powers with t * a**0 == t
    pow.SetLength(l);
    pow[0] = t;

    if (need_upper_products)
    {
        upper.SetLength(l);
        long n = deg(g);
        for (long i = 1; i < l; i++)
        {
            MulMod_local(pow[i], upper[i-1], pow[i-1], a_pow_mul, g_mod);
            reverse(upper[i-1], upper[i-1], n-2);
        }
        zz_pX tmp;
        MulMod_local(tmp, upper[l-1], pow[l-1], a_pow_mul, g_mod);
        reverse(upper[l-1], upper[l-1], n-2);
    }
    else
        for (long i = 1; i < l; ++i)
            MulMod(pow[i], pow[i-1], a_pow_mul, g_mod);

}

/*------------------------------------------------------------*/
/* TODO: what is this?                                        */
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
/* TODO: what is this?                                        */
/*------------------------------------------------------------*/
static void get_quos (Mat<zz_pX> &quos,
                      const Vec<zz_pX> &alphas,
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
    cout << "mat mul ("  << mat1.NumRows() << " x " << mat2.NumRows() << " x " << mat2.NumCols() << ", deg " << max(deg(mat1),deg(mat2)) << "): " << GetWallTime()-t << endl;

    t = GetWallTime();
    trunc(res1, res1, 2*m);
    RightShift(res2, res2, 2*m);
    quos = res1 + res2;
    reverse(quos, quos, 2 * m - 1);
    cout << "clean up: " << GetWallTime()-t << endl;
}

/*------------------------------------------------------------*/
/* TODO: what is this?                                        */
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
/* TODO: what is this?                                        */
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

/*------------------------------------------------------------*/
/* TODO: what is this?                                        */
/*------------------------------------------------------------*/
void gen_sequence (Coeffs &res,
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
    get_quos(quos, alphas, upper, g, m);
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
/* TODO: what is this?                                        */
/*------------------------------------------------------------*/
void gen_sequence (Vec<Mat<zz_p>> &mats,
                   const zz_pX &a,
                   const zz_pX &g,
                   const long m)
{
    const long n = deg(g);
    const long d = 2*ceil((n*1.0)/m)+1;
    const long sqrt_d = ceil(sqrt(d));

    double t;

    t = GetWallTime();
	Vec<Coeffs> res;
    res.SetLength(m);
    zz_pX x_poly;
    SetCoeff(x_poly, m-1, zz_p(1));
    gen_sequence(res[0],x_poly,a,g,m);
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
/* TODO: what is this?                                        */
/*------------------------------------------------------------*/
void gen_sequence_naive (Vec<Coeffs> &res,
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
/* TODO: what is this?                                        */
/*------------------------------------------------------------*/
void matrix_recon_approximation(Mat<zz_pX> &basis,
                                const Vec<Mat<zz_p>> &seq)
{
    Mat<zz_pX> pmat;
    conv(pmat, seq, seq.length());

    // reconstruct fraction
    matrix_pade_generic(basis,pmat,seq.length());
}

/*------------------------------------------------------------*/
/* TODO: what is this?                                        */
/*------------------------------------------------------------*/
void matrix_recon_interpolation(Mat<zz_pX> & basis,
                                const Vec<zz_p> & pts,
                                Vec<Mat<zz_p>> & seq)
{
    // length of sequence
    const long len = seq.length();
    if (len == 0)
        Error("empty sequence for matrix reconstruction");

    // dimension
    const long m = seq[0].NumRows();
    // uniform shift (0,..,0)
    VecLong shift(2*m);

    // add identity at the bottom of each matrix in seq
    for (long j = 0; j < len; ++j)
    {
        seq[j].SetDims(2*m, m);
        for (long i = 0; i < m; ++i)
            set(seq[j][m+i][i]);
    }

    // call pmbasis
    Mat<zz_pX> intbas;
    pmbasis(intbas, seq, pts, shift, 0, pts.length());

    basis.SetDims(m, m);
    for (long i = 0; i < m; ++i)
        for (long j = 0; j < m; ++j)
            basis[i][j].swap(intbas[i][j]);
}

/*------------------------------------------------------------*/
/* matrix reconstruction by geometric interpolants            */
/* TODO: why do we need pts and r?                            */
/*------------------------------------------------------------*/
void matrix_recon_interpolation_geometric(Mat<zz_pX> & basis,
                                          const Vec<zz_p> & pts,
                                          const zz_p & r,
                                          Vec<Mat<zz_p>> & seq)
{
    // length of sequence
    const long len = seq.length();
    if (len == 0)
        Error("empty sequence for matrix reconstruction");

    // dimension
    long m = seq[0].NumRows();
    // uniform shift (0,..,0)
    VecLong shift(2*m);

    // add identity at the bottom of each matrix in seq
    for (long j = 0; j < len; j++)
    {
        seq[j].SetDims(2*m, m);
        for (long i = 0; i < m; i++)
            set(seq[j][m+i][i]);
    }

    // call pmbasis
    Mat<zz_pX> intbas;
    pmbasis_geometric(intbas, seq, pts, r, shift, 0, len);

    basis.SetDims(m, m);
    for (long i = 0; i < m; ++i)
        for (long j = 0; j < m; ++j)
            basis[i][j].swap(intbas[i][j]);
}


// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
