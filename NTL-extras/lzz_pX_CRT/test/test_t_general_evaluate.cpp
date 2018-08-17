#include <NTL/lzz_pX.h>
#include <NTL/vector.h>
#include <iomanip>

#include "vec_lzz_p_extra.h"
#include "lzz_pX_CRT.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* computes rnd . eval(f) using direct / tranposed eval       */
/*------------------------------------------------------------*/
void check()
{
    long p = 1125899906842679;
    zz_p::init(p);

    for (long j = 1; j < 500; j++)
    {
        zz_p res1, res2;
        zz_pX f, g;
        zz_pX_Multipoint_General ev;
        Vec<zz_p> valF, rnd;

        ev = get_general_points(j);
        f = random_zz_pX(j);
        random_vec_zz_p(rnd, j);

        ev.evaluate(valF, f);
        ev.t_evaluate(g, rnd);

        res1 = 0;
        for (long i = 0; i < j; i++)
            res1 += rnd[i] * valF[i];

        res2 = 0;
        for (long i = 0; i < j; i++)
            res2 += coeff(f, i) * coeff(g, i);

        if (res1 != res2) 
        {
            cerr << "error for j=" << j << endl;
            exit (-1);
        }
    }
}  

/*------------------------------------------------------------*/
/* main calls check                                           */
/*------------------------------------------------------------*/
int main(int argc, char ** argv)
{
    check();
    return 0;
}
