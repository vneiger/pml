#include <NTL/lzz_pX.h>
#include <NTL/vector.h>

#include "lzz_pX_CRT.h"
#include "vec_lzz_p_extra.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* creates, t-evaluates, t-interpolates a vector of polys.    */
/*------------------------------------------------------------*/
void check(){

    long p = 1125899906842679;
    zz_p::init(p);

    for (long j = 1; j < 500; j+=1)
    {
        Vec<Vec<zz_p>> u, valF, v;
        Vec<zz_pX> f, tvalU;
        zz_pX_Multipoint_Geometric points;
        long len;

        len = 10;
        u.SetLength(j);
        f.SetLength(len);

        for (long i = 0; i < j; i++)
            u[i].SetLength(len);

        for (long r = 0; r < len; r++)
        {
            for (long i = 0; i < j; i++)
                u[i][r] = random_zz_p();
            f[r] = random_zz_pX(j);
        }

        points = get_geometric_points(j);

        points.t_evaluate_vector(tvalU, u);
        points.evaluate_vector(valF, f);


        for (long r = 0; r < len; r++)
        {
            zz_p res1, res2;
            res1 = 0;
            for (long i = 0; i < j; i++)
                res1 += u[i][r] * valF[i][r];

            res2 = 0;
            for (long i = 0; i < j; i++)
                res2 += coeff(f[r], i) * coeff(tvalU[r], i);

            if (res1 != res2) 
            {
                cerr << "error vector t_evaluate for j=" << j << endl;
                exit (-1);
            }
        }

        points.t_interpolate_vector(v, tvalU);

        if (u != v)
        {
            cerr << "error with vector t_interpolate at j=" << j << endl;
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
