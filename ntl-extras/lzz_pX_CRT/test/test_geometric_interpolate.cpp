#include <NTL/lzz_pX.h>
#include <NTL/vector.h>

#include "lzz_pX_CRT.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* tests geometric evaluation, then interpolation             */
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
        zz_pX f, g;
        zz_pX_Multipoint_Geometric ev;
        Vec<zz_p> val;

        a = random_zz_p();
        b = random_zz_p();
        f = random_zz_pX(j);

        ev = zz_pX_Multipoint_Geometric(a, j);
        ev.set_FFT_evaluate();
        ev.set_FFT_interpolate();
        ev.evaluate(val, f);
        ev.interpolate(g, val);

        if (g != f) 
        {
            cerr << "error for geometric interpolate at j=" << j << endl;
            exit(-1);
        }

        ev = zz_pX_Multipoint_Geometric(a, j);
        ev.unset_FFT_evaluate();
        ev.unset_FFT_interpolate();
        ev.evaluate(val, f);
        ev.interpolate(g, val);

        if (g != f) 
        {
            cerr << "error for geometric interpolate at j=" << j << endl;
            exit(-1);
        }

        ev = zz_pX_Multipoint_Geometric(a, b, j);
        ev.set_FFT_evaluate();
        ev.set_FFT_interpolate();
        ev.evaluate(val, f);
        ev.interpolate(g, val);

        if (g != f) 
        {
            cerr << "error for geometric interpolate at j=" << j << endl;
            exit(-1);
        }

        ev = zz_pX_Multipoint_Geometric(a, b, j);
        ev.unset_FFT_evaluate();
        ev.unset_FFT_interpolate();
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
    check(0);
    check(23068673);
    check(288230376151711813);
    return 0;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
