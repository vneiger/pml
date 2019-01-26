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

#include "util.h"
#include "mat_lzz_pX_extra.h"
#include "mat_lzz_pX_sequence.h"
#include "structured_lzz_p.h"

NTL_CLIENT

int main(int argc, char *argv[])
{
    if (argc != 3) throw std::invalid_argument("usage: ./time_compare m d");

    zz_p::FFTInit(0);

    long m = atoi(argv[1]);
    long d = atoi(argv[2]);

    std::cout << "Comparing:" << std::endl;
    std::cout << "-- computing denominator P such that deg(PF) < d, where F random m x m matrix of degree < 2d (done with matrix Pade generic code: timings include checking the result is correct)" << std::endl;
    std::cout << "-- computing approximant basis of [[F], [-1]] at order 2d, with F random m x m matrix of degree < 2d" << std::endl;
    std::cout << "-- solving a block Hankel system with m x m random Hankel blocks, each d x d" << std::endl;
    std::cout << "-- m = " << m << ", d = " << d << std::endl;

    double t, tt;
    long nb_iter;
    bool matpade = true;

    // generic matrix Pade
    nb_iter = 0; t = 0.0;
    while (t<0.5 && matpade)
    {
        // set pmat random m x m of degree < 2d
        Mat<zz_pX> pmat;
        random(pmat, m, m, 2*d);

        // compute fraction denominator
        tt = GetWallTime();
        Mat<zz_pX> denom;
        matrix_pade_generic(denom, pmat, 2*d);
        Mat<zz_pX> prod;
        mul_trunc(prod, denom, pmat, 2*d);
        if ((deg(denom) > d) || deg(prod) >= d)
            matpade = false;
        t += GetWallTime() - tt;
        ++nb_iter;
    }
    cout << "matrix Pade:\t" << t/nb_iter << (matpade ? " (correct)" : " (error)" ) << endl;

    // pmbasis
    nb_iter = 0; t = 0.0;
    while (t<0.5)
    {
        // set pmat to [[F], [-1]] with F random m x m of degree < 2d
        Mat<zz_pX> pmat;
        random(pmat, m, m, 2*d);
        pmat.SetDims(2*m, m);
        for (long i = 0; i < m; ++i)
            pmat[i+m][i] = -1;

        // compute approx basis
        tt = GetWallTime();
        VecLong shift(2*m,0);
        Mat<zz_pX> appbas;
        pmbasis(appbas, pmat, 2*d, shift);
        t += GetWallTime() - tt;
        ++nb_iter;
    }
    cout << "pmbasis:\t" << t/nb_iter << endl;

    // structured matrix
    nb_iter = 0; t = 0.0;
    while (t<0.5)
    {
        Vec<Vec<hankel_lzz_p>> h_vec; //set to m blocks by m blocks
        h_vec.SetLength(m);
        for (long i = 0; i < m ; i++)
            h_vec[i].SetLength(m);

        for (long i = 0; i < m; i++)
            for (long j = 0; j < m; j++)
                h_vec[i][j] = hankel_lzz_p(random_vec_zz_p(2*d-1),d,d);

        mosaic_hankel_lzz_p mh{h_vec};

        Vec<zz_p> b;
        random(b, m*d);

        tt = GetWallTime();
        Vec<zz_p> x;
        mh.solve(x,b);
        t += GetWallTime() - tt;

        ++nb_iter;
    }

    cout << "block Hankel:\t" << t/nb_iter << endl;

    std::cout << std::endl << "Comparing:" << std::endl;
    std::cout << "-- computing approximant basis of [[F], [-1]] at order (m+1)d, with F random 1 x m matrix of degree < (m+1)d" << std::endl;
    std::cout << "-- solving a block Hankel system with 1 x m random Hankel blocks, each md x d" << std::endl;
    std::cout << "-- m = " << m << ", d = " << d << std::endl;


    // pmbasis
    nb_iter = 0; t = 0.0;
    while (t<0.5)
    {
        // set pmat to [[F], [-1]] with F random 1 x m of degree < 2md
        Mat<zz_pX> pmat;
        random(pmat, 1, m, (m+1)*d);
        pmat.SetDims(m+1, m);
        for (long i = 0; i < m; ++i)
            pmat[i+1][i] = -1;

        // compute approx basis
        tt = GetWallTime();
        VecLong shift(m+1,0);
        Mat<zz_pX> appbas;
        pmbasis(appbas, pmat, (m+1)*d, shift);
        t += GetWallTime() - tt;
        ++nb_iter;
    }
    cout << "pmbasis:\t" << t/nb_iter << endl;

    // structured matrix
    nb_iter = 0; t = 0.0;
    while (t<0.5)
    {
        Vec<Vec<hankel_lzz_p>> h_vec; //set to m blocks by 1 blocks
        h_vec.SetLength(1);
        //for (long i = 0; i < 1 ; i++)
            h_vec[0].SetLength(m);

        //for (long i = 0; i < 1; i++)
            for (long j = 0; j < m; j++)
                h_vec[0][j] = hankel_lzz_p(random_vec_zz_p(m*d+d-1),m*d,d);

        mosaic_hankel_lzz_p mh{h_vec};

        Vec<zz_p> b;
        random(b, m*d);

        tt = GetWallTime();
        Vec<zz_p> x;
        mh.solve(x,b);
        t += GetWallTime() - tt;

        ++nb_iter;
    }

    cout << "block Hankel:\t" << t/nb_iter << endl;
}


// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
