#include <NTL/lzz_pX.h>
#include <NTL/vector.h>

#include "lzz_pX_CRT.h"
#include "vec_lzz_p_extra.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* creates, t-evaluates, t-interpolates a matrix of polys.    */
/*------------------------------------------------------------*/
void check(){

    long p = 1125899906842679;
    zz_p::init(p);

    for (long j = 1; j < 500; j+=1)
    {
        Vec<Mat<zz_p>> u, valF, v;
        Mat<zz_pX> f, tvalU;
        zz_pX_Multipoint_Geometric points;
        long rdim, cdim;

        rdim = 10;
        cdim = 3;

        u.SetLength(j);
        f.SetDims(rdim, cdim);

        for (long i = 0; i < j; i++)
            u[i].SetDims(rdim, cdim);

        for (long r = 0; r < rdim; r++)
        {
            for (long s = 0; s < cdim; s++)
            {
                for (long i = 0; i < j; i++)
                    u[i][r][s] = random_zz_p();
                f[r][s] = random_zz_pX(j);
            }
        }

        points = get_geometric_points(j);

        points.t_evaluate_matrix(tvalU, u);
        points.evaluate_matrix(valF, f);

        for (long r = 0; r < rdim; r++)
        {
            for (long s = 0; s < cdim; s++)
            {
                zz_p res1, res2;
                res1 = 0;
                for (long i = 0; i < j; i++)
                    res1 += u[i][r][s] * valF[i][r][s];

                res2 = 0;
                for (long i = 0; i < j; i++)
                    res2 += coeff(f[r][s], i) * coeff(tvalU[r][s], i);

                if (res1 != res2) 
                {
                    cerr << "error matrix t_evaluate for j=" << j << endl;
                    exit (-1);
                }
            }
        }

        points.t_interpolate_matrix(v, tvalU);

        if (u != v)
        {
            cerr << "error with matrix t_interpolate at j=" << j << endl;
        }
    }
}

/*------------------------------------------------------------*/
/* main calls check                                           */
/*------------------------------------------------------------*/
int main(int argc, char ** argv){
    check();
    return 0;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
