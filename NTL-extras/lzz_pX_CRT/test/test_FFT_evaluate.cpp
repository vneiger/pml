#include <NTL/lzz_pX.h>
#include <NTL/vector.h>

#include "lzz_pX_CRT.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* compares FFT evaluation to subproduct tree one             */
/*------------------------------------------------------------*/
void check()
{
    zz_p::FFTInit(0);

    for (long j = 1; j < 500; j++)
    {
        zz_pX f;
        zz_pX_Multipoint_General evG;
        zz_pX_Multipoint_FFT ev;
        Vec<zz_p> q, val, valG;

        ev = get_FFT_points(j);
        long n = ev.length();

        q.SetLength(n);
        for (long i = 0; i < n; i++)
            ev.get_point(q[i], i);
        evG = zz_pX_Multipoint_General(q);

        f = random_zz_pX(j);
        evG.evaluate(valG, f);
        ev.evaluate(val, f);

        if (valG != val) 
        {
            cerr << "error for FFT evaluate at j=" << j << endl;
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
