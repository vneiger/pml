#include <NTL/vec_lzz_p.h>
#include <assert.h>

#include "vec_lzz_p_extra.h"
#include "lzz_pX_CRT.h"
#include "structured_lzz_p.h"

PML_CLIENT

/* BEGIN TAKEN FROM src/mosaic_hankel_lzz_p.cpp */
/*------------------------------------------------------------------*/
/* finds c such that                                                */
/* - c != 0                                                         */
/* - a^i - c a^j != 0 for 0 <= i < n and 0 <= j < m                 */
/*------------------------------------------------------------------*/
static void find_c(zz_p& c, const zz_p& a, long n, long m)
{
    Vec<zz_p> pow_a, pow_inva;
    pow_a.SetLength(n);
    pow_inva.SetLength(m);

    pow_a[0] = to_zz_p(1);
    pow_inva[0] = to_zz_p(1);

    zz_p inva = 1/a;
    for (long i = 1; i < n; i++)
        pow_a[i] = a*pow_a[i-1];
    for (long i = 1; i < m; i++)
        pow_inva[i] = inva*pow_inva[i-1];

    bool done;
    do
    {
        c = random_zz_p();
        done = true;
        if (c == 0)
            done = false;
        for (long i = 0; i < n; i++)
            if (c == pow_a[i])
                done = false;
        for (long i = 0; i < m; i++)
            if (c == pow_inva[i])
                done = false;
    } 
    while (done != true);
}

/*------------------------------------------------------------------*/
/* preconditions M                                                  */
/* builds the matrix CL = (D_e X_int) M (D_f Y_int)^t, where:       */
/* - X_int, Y_int are geometric interpolation                       */
/* - D_e, D_f are diagonal matrix built on vectors e and f          */
/* - CL is cauchy-like geometric                                    */
/* - CL is expected to have generic rank profile                    */
/* return value is 0 if field to small for preconditioning          */
/*------------------------------------------------------------------*/
static long to_cauchy_grp(cauchy_like_geometric_lzz_p& CL, 
                          zz_pX_Multipoint_Geometric& X_int, zz_pX_Multipoint_Geometric& Y_int,
                          Vec<zz_p> &e, Vec<zz_p> &f,
                          const mosaic_hankel_lzz_p& M)
{
    Mat<zz_p> X, Y;
    Mat<zz_p> G, H;
    M.phi_plus_generators(G, H);
    zz_p a, b, c;
    long n = M.NumRows();
    long m = M.NumCols();

    if (max(m, n) > zz_p::modulus()/ 10)
    {
        return 0;
    }

    element_of_order(a, 2 * max(m, n));
    b = a*a;
    find_c(c, b, n, m);
    X_int = zz_pX_Multipoint_Geometric(a, to_zz_p(1), n);
    Y_int = zz_pX_Multipoint_Geometric(a, c, m);

    long alpha = G.NumCols();
    X.SetDims(n, alpha+2);
    Y.SetDims(m, alpha+2);

    e.SetLength(n);
    for (long i = 0; i < n; i++)
        e[i] = random_zz_p();

    f.SetLength(m);
    for (long i = 0; i < m; i++)
        f[i] = random_zz_p();

    Vec<zz_p> tmp_v;
    for (long j = 0; j < alpha; j++)
    {
        zz_pX tmp_p;
        tmp_p.rep.SetLength(n);
        zz_p* coef_p = tmp_p.rep.elts();
        for (long i = 0; i < n; i++)
        {
            coef_p[i] = G[i][j];
        }
        tmp_p.normalize();
        X_int.evaluate(tmp_v, tmp_p);
        for (long i = 0; i < n; i++)
        {
            X[i][j] = tmp_v[i] * e[i];
        }
    }

    zz_p tmp_z = to_zz_p(1);
    for (long i = 0; i < n; i++)
    {
        X[i][alpha] = e[i]*(power(tmp_z, n)-1);
        tmp_z = tmp_z * b;
    }

    M.last_column_of_block(tmp_v, M.NumBlockCols()-1);
    zz_pX tmp_p;
    tmp_p.rep.SetLength(n);
    zz_p* coef_p = tmp_p.rep.elts();
    for (long i = 0; i < n; i++)
    {
        coef_p[i] = tmp_v[i];
    }
    tmp_p.normalize();
    X_int.evaluate(tmp_v, tmp_p);
    for (long i = 0; i < n; i++)
    {
        X[i][alpha+1] = tmp_v[i] * e[i];
    }

    Vec<zz_p> tmp_w;
    for (long j = 0; j < alpha; j++)
    {
        zz_pX tmp_q;
        tmp_q.rep.SetLength(m);
        zz_p* coef_q = tmp_q.rep.elts();
        for (long i = 0; i < m; i++)
        {
            coef_q[i] = H[i][j];
        }
        tmp_q.normalize();
        Y_int.evaluate(tmp_w, tmp_q);
        for (long i = 0; i < m; i++)
        {
            Y[i][j] = tmp_w[i] * f[i];
        }
    }

    M.last_row_of_block(tmp_w, M.NumBlockRows()-1);
    zz_pX tmp_q;
    tmp_q.rep.SetLength(m);
    zz_p* coef_q = tmp_q.rep.elts();
    for (long i = 0; i < m; i++)
    {
        coef_q[i] = tmp_w[i];
    }
    tmp_q.normalize();
    Y_int.evaluate(tmp_w, tmp_q);
    for (long i = 0; i < m; i++)
    {
        Y[i][alpha] = tmp_w[i] * f[i];
    }

    tmp_z = c;
    for (long i = 0; i < m; i++)
    {
        Y[i][alpha+1] = -f[i]*power(tmp_z, m);
        tmp_z = tmp_z * b;
    }

    CL = cauchy_like_geometric_lzz_p(X, Y, to_zz_p(1), c, b);
    return 1;
}
/* END TAKEN FROM src/mosaic_hankel_lzz_p.cpp */


/*------------------------------------------------------------*/
/* returns the diagonal matrix for e                          */
/*------------------------------------------------------------*/
void to_dense(Mat<zz_p>& M, const Vec<zz_p>& e)
{
    long n = e.length();
    M.SetDims(n, n);
    for (long i = 0; i < n; i++)
        for (long j = 0; j < n; j++)
            M[i][j] = 0;
    for (long i = 0; i < n; i++)
        M[i][i] = e[i];
}

/*------------------------------------------------------------*/
/* returns the companion matrix of x^n-c                      */
/*------------------------------------------------------------*/
void do_Z(Mat<zz_p>& M, long n, const zz_p& c)
{
    M.SetDims(n, n);
    for (long i = 0; i < n; i++)
        for (long j = 0; j < n; j++)
            M[i][j] = 0;
    for (long i = 0; i < n-1; i++)
        M[i+1][i] = 1;
    M[0][n-1] = c;
}

/*------------------------------------------------------------*/
/* creates mosaic hankel matrices, transforms to cauchy grp   */
/*------------------------------------------------------------*/
void check(long p)
{
    if (p == 0)
        zz_p::FFTInit(0);
    else
        zz_p::init(p);

    for (long i = 2; i < 110; i += 1)
    {

        Vec<zz_p> dat00, dat01, dat02, dat10, dat11, dat12;

        dat00.SetLength(2 + i - 1);
        for (long k = 0; k < 2 + i -1; k++)
            dat00[k] = 0;
        random(dat01, 2 + 2 - 1);
        random(dat02, 2 + i - 1);
        random(dat10, i-1 + i - 1);
        random(dat11, i-1 + 2 - 1);
        random(dat12, i-1 + i - 1);

        hankel_lzz_p h00(dat00, 2, i), h01(dat01, 2, 2), h02(dat02, 2, i), h10(dat10, i-1, i), h11(dat11, i-1, 2), h12(dat12, i-1, i);
        Vec<hankel_lzz_p> row0, row1;

        row0.SetLength(3);
        row0[0] = h00;
        row0[1] = h01;
        row0[2] = h02;
        row1.SetLength(3);
        row1[0] = h10;
        row1[1] = h11;
        row1[2] = h12;
        Vec< Vec<hankel_lzz_p> > H;
        H.SetLength(2);
        H[0] = row0;
        H[1] = row1;

        mosaic_hankel_lzz_p MH;
        MH = mosaic_hankel_lzz_p(H);

        cauchy_like_geometric_lzz_p CL;
        Mat<zz_p> X, Y;
        Vec<zz_p> e, f;
        zz_pX_Multipoint_Geometric X_int, Y_int;
        
        to_cauchy_grp(CL, X_int, Y_int, e, f, MH);
        
        Mat<zz_p> MM, MC, MX_i, MY_i, Me, Mf;
        MM = MH.to_dense();
        MC = CL.to_dense();
        X_int.to_dense(MX_i);
        Y_int.to_dense(MY_i);
        to_dense(Me, e);
        to_dense(Mf, f);

        Mat<zz_p> T = Me*MX_i*MM*transpose(Mf*MY_i);
        assert (T == MC);
    }
}


/*------------------------------------------------------------*/
/* main just calls check()                                    */
/*------------------------------------------------------------*/
int main(int argc, char** argv)
{
    check(0);
    check(786433);
    check(288230376151711813);
    return 0;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
