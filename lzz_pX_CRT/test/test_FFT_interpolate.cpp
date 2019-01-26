#include <NTL/lzz_pX.h>
#include <NTL/vector.h>
#include <NTL/version.h>

#include "lzz_pX_CRT.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* tests FFT evaluation, then interpolation                   */
/*------------------------------------------------------------*/
void check()
{
    zz_p::FFTInit(0);

    for (long j = 1; j < 500; j++)
    {
        zz_pX f, g;
        zz_pX_Multipoint_FFT ev;
        Vec<zz_p> val;

        f = random_zz_pX(j);

        ev = get_FFT_points(j);
        ev.evaluate(val, f);
        ev.interpolate(g, val);

        if (g != f) 
        {
            cerr << "error for FFT interpolate at j=" << j << endl;
            exit(-1);
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
