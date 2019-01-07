#include <NTL/lzz_pX.h>
#include <NTL/vector.h>
#include <memory>

#include "vec_lzz_p_extra.h"
#include "lzz_pX_CRT.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* does a multipoint evaluation                               */
/*------------------------------------------------------------*/
void check(){

    long p = 1125899906842679;
    zz_p::init(p);

    for (long j = 1; j < 500; j++)
    {
        zz_pX f;
        zz_pX_Multipoint_General ev;
        Vec<zz_p> q, val;

        random(q, j);
        ev = zz_pX_Multipoint_General(q);
        f = random_zz_pX(2*j);

        ev.evaluate(val, f);
        for (long i = 0; i < j; i++)
        {
            if (val[i] != eval(f, q[i]))
            {
                cerr << "error Multipoint_General evaluate with j=" << j << endl;
                exit(-1);
            }
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
