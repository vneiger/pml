#include <NTL/lzz_pX.h>
#include <NTL/vector.h>
#include <iomanip>

#include "lzz_pX_CRT.h"
#include "vec_lzz_p_extra.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* does a direct and inverse transposed geometric evaluation  */
/*------------------------------------------------------------*/
void check(long p)
{
    if (p == 0)
        zz_p::FFTInit(0);
    else
        zz_p::init(p);

    for (long j = 1; j < 500; j++)
    {
        zz_p a, b;
        zz_pX g;
        zz_pX_Multipoint_Geometric ev;
        Vec<zz_p> rnd1, rnd2;

        a = random_zz_p();
        b = random_zz_p();

        ev = zz_pX_Multipoint_Geometric(a, j);
        random_vec_zz_p(rnd1, j);
        ev.t_evaluate(g, rnd1);
        ev.t_interpolate(rnd2, g);
        if (rnd1 != rnd2) 
        {
            cerr << "error for j=" << j << endl;
            exit (-1);
        }

        ev = zz_pX_Multipoint_Geometric(a, b, j);
        random_vec_zz_p(rnd1, j);
        ev.t_evaluate(g, rnd1);
        ev.t_interpolate(rnd2, g);
        if (rnd1 != rnd2) 
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
    check(0);
    check(1125899906842679);
    return 0;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
