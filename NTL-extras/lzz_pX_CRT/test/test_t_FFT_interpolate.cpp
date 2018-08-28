#include <NTL/lzz_pX.h>
#include <NTL/vector.h>
#include <iomanip>

#include "util.h"
#include "lzz_pX_CRT.h"
#include "vec_lzz_p_extra.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* does a direct and inverse transposed FFT evaluation        */
/*------------------------------------------------------------*/
void check()
{
    zz_p::FFTInit(0);

    for (long j = 1; j < 500; j++)
    {
        zz_pX g;
        zz_pX_Multipoint_FFT ev;
        Vec<zz_p> rnd1, rnd2;

        ev = get_FFT_points(j);
        random_vec_zz_p(rnd1, ev.length());

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
