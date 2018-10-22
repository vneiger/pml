#include <NTL/vec_lzz_p.h>
#include <assert.h>

#include "lzz_p_extra.h"
#include "vec_lzz_p_extra.h"
#include "structured_lzz_p.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* ensures that G and H have generic rank profiles            */
/*------------------------------------------------------------*/
void generate_generic_rank_profile(Mat<zz_p>& G, Mat<zz_p>& H, long rk)
{
    long n = G.NumRows();
    long m = H.NumRows();
    long alpha = G.NumCols();

    for (long i = 0; i < rk; i++)
        for (long j = 0; j < alpha; j++)
            G[i][j] = random_zz_p();
    for (long i = rk; i < n; i++)
        for (long j = 0; j < alpha; j++)
            G[i][j] = 0;

    for (long i = 0; i < rk; i++)
        for (long j = 0; j < alpha; j++)
            H[i][j] = random_zz_p();
    for (long i = rk; i < m; i++)
        for (long j = 0; j < alpha; j++)
            H[i][j] = 0;
}

/*------------------------------------------------------------*/
/* ensures that G and H have not generic rank profile         */
/* puts some rows fulls of zeros, then a non-zero one         */
/*------------------------------------------------------------*/
void generate_non_generic_rank_profile(Mat<zz_p>& G, Mat<zz_p>& H, long rk, long extra)
{
    long alpha = G.NumCols();

    if (extra <= rk)
    {
        Error("Matrix will have generic rank profile");
    }

    generate_generic_rank_profile(G, H, rk);

    for (long j = 0; j < alpha; j++){
        G[extra][j] = random_zz_p();
        H[extra][j] = random_zz_p();
    }
}


/*------------------------------------------------------------*/
/* test in size (i,j), with rank rk                           */
/* tests a generic rank profile matrix                        */
/* tests a non-generic rank profile one                       */
/*------------------------------------------------------------*/
void test(long i, long j, long alpha, long rk)
{
    if (rk < 0 || rk > i || rk > j)
    {
        return;
    }

    cauchy_like_geometric_lzz_p M;
    Mat<zz_p> A, B;
    long r;
    zz_p a;

    a = element_of_order(2*max(i, j));
    A.SetDims(i, alpha);
    B.SetDims(j, alpha);

    // 1. generic rank profile
    {
        generate_generic_rank_profile(A, B, rk);
        M = cauchy_like_geometric_lzz_p(A, B, to_zz_p(1), power(a, i), a);
        Vec<zz_p> u, v, w;
        u = random_vec_zz_p(j);
        v = M.mul_right(u);
        
        r = M.solve_grp(w, v, 50, 8);
        assert (r == 1);
        assert (M.mul_right(w) == v);

        r = M.solve_grp(w, v, 10, 15);
        assert (r == 1);
        assert (M.mul_right(w) == v);
        
        r = M.solve_grp(w, v, 10000, 15);
        assert (r == 1);
        assert (M.mul_right(w) == v);
    }

    // 2. non-generic rank profile
    if (rk < min(i, j)-1) 
    {
        if ((rk+5 < i) && (rk+5 < j))
            generate_non_generic_rank_profile(A, B, rk, rk + 5);
        else
            generate_non_generic_rank_profile(A, B, rk, min(i-1, j-1));

        M = cauchy_like_geometric_lzz_p(A, B, to_zz_p(1), power(a, i), a);
        Vec<zz_p> u, v, w;
        u = random_vec_zz_p(j);
        v = M.mul_right(u);

        r = M.solve_grp(w, v, 50, 8);
        assert (r == 0);

        r = M.solve_grp(w, v, 10, 15);
        assert (r == 0);

        r = M.solve_grp(w, v, 10000, 15);
        assert (r == 0);
    }
}


/*------------------------------------------------------------*/
/* does some right matrix / vector multiplications            */
/*------------------------------------------------------------*/
void check(long p)
{
    if (p == 0)
        zz_p::FFTInit(0);
    else
        zz_p::init(p);

    for (long i = 1; i < 1000; i += (i < 100 ? 1 : 50))
    {
        long j;
        j = max(1, i-2);
        test(i, j, 5, max(1, min(i-2, 19)));
        test(i, j, min(30, max(1, min(i,j)/2)), max(1, min(i-2, 17)));
        test(i, j, min(30, max(1, min(i,j)-7)), max(1, min(i-2, 22)));
        j = i+2;
        test(i, j, 5, max(1, min(i-2, 19)));
        test(i, j, min(30, max(1, min(i,j)/2)), max(1, min(i-2, 17)));
        test(i, j, min(30, max(1, min(i,j)-7)), max(1, min(i-2, 22)));
        j = max(1, i/2);
        test(i, j, 5, max(1, min(i-2, 19)));
        test(i, j, min(30, max(1, min(i,j)/2)), max(1, min(i-2, 17)));
        test(i, j, min(30, max(1, min(i,j)-7)), max(1, min(i-2, 22)));
        j = i*2;
        test(i, j, 5, max(1, min(i-2, 19)));
        test(i, j, min(30, max(1, min(i,j)/2)), max(1, min(i-2, 17)));
        test(i, j, min(30, max(1, min(i,j)-7)), max(1, min(i-2, 22)));
        j = i;
        test(i, j, 5, max(1, min(i-2, 19)));
        test(i, j, min(30, max(1, min(i,j)/2)), max(1, min(i-2, 17)));
        test(i, j, min(30, max(1, min(i,j)-7)), max(1, min(i-2, 22)));
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
