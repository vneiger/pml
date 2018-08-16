#include <NTL/lzz_pX.h>
#include <NTL/vector.h>

#include "lzz_pX_CRT.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* tests geometric evaluation, then interpolation             */
/*------------------------------------------------------------*/
void check()
{
    long p = 1125899906842679;
    zz_p::init(p);

    for (long j = 1; j < 500; j++)
    {
        zz_p a;
        zz_pX f, g;
        zz_pX_Multipoint_Geometric ev;
        Vec<zz_p> val;

        a = random_zz_p();
        f = random_zz_pX(j);

        ev = zz_pX_Multipoint_Geometric(a, j);
        ev.evaluate(val, f);
        ev.interpolate(g, val);

        if (g != f) 
        {
            cerr << "error for geometric interpolate at j=" << j << endl;
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
