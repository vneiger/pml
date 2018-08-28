#include <NTL/lzz_pX.h>
#include <NTL/vector.h>

#include "lzz_pX_CRT.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* creates, evaluates and interpolates a vector of polynomials*/
/*------------------------------------------------------------*/
void check(){

    long p = 1125899906842679;
    zz_p::init(p);

    for (long j = 1; j < 500; j+=1)
    {
        Vec<zz_pX> u;
        long len = 10;
        u.SetLength(len);

        for (long r = 0; r < len; r++)
        {
            u[r] = random_zz_pX(j);
        }

        Vec<Vec<zz_p>> valU;

        zz_pX_Multipoint_Geometric points = get_geometric_points(j);
        points.evaluate_vector(valU, u);

        for (long r = 0; r < len; r++)
        {
            zz_p pt;
            for (long s = 0; s < points.length(); s++)
            {
                points.get_point(pt, s);
                if (eval(u[r], pt) != valU[s][r])
                {
                    cerr << "error with vector evaluate at j=" << j << endl;
                    exit(-1);
                }
            }
        }

        Vec<zz_pX> v;
        points.interpolate_vector(v, valU);

        if (u != v)
        {
            cerr << "error with vector interpolate at j=" << j << endl;
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
