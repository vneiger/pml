#include <NTL/lzz_pX.h>
#include <NTL/vector.h>
#include <NTL/matrix.h>

#include "lzz_pX_CRT.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* creates, evaluates and interpolates a matrix of polynomials*/
/*------------------------------------------------------------*/
void check(){

    long p = 1125899906842679;
    zz_p::init(p);

    for (long j = 1; j < 500; j+=1)
    {
        Mat<zz_pX> u;
        long rdim = 10;
        long cdim = 3;
        u.SetDims(rdim, cdim);

        for (long r = 0; r < rdim; r++)
        {
            for (long c = 0; c < cdim; c++)
            {
                u[r][c] = random_zz_pX(j);
            }
        }

        Vec<Mat<zz_p>> valU;

        zz_pX_Multipoint_Geometric points = get_geometric_points(j);
        points.evaluate_matrix(valU, u);

        for (long r = 0; r < rdim; r++)
        {
            for (long c = 0; c < cdim; c++)
            {
                zz_p pt;
                for (long s = 0; s < points.length(); s++)
                {
                    points.get_point(pt, s);
                    if (eval(u[r][c], pt) != valU[s][r][c])
                    {
                        cerr << "error with matrix evaluate at j=" << j << endl;
                        exit(-1);
                    }
                }
            }
        }       

        Mat<zz_pX> v;
        points.interpolate_matrix(v, valU);

        if (u != v)
        {
            cerr << "error with matrix interpolate at j=" << j << endl;
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
