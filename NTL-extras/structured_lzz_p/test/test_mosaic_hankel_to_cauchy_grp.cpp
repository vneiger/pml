#include <NTL/vec_lzz_p.h>
#include <assert.h>

#include "vec_lzz_p_extra.h"
#include "lzz_pX_CRT.h"
#include "structured_lzz_p.h"

NTL_CLIENT


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
        random_vec_zz_p(dat01, 2 + 2 - 1);
        random_vec_zz_p(dat02, 2 + i - 1);
        random_vec_zz_p(dat10, i-1 + i - 1);
        random_vec_zz_p(dat11, i-1 + 2 - 1);
        random_vec_zz_p(dat12, i-1 + i - 1);

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
