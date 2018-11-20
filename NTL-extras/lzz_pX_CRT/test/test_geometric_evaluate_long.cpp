#include <NTL/lzz_pX.h>
#include <NTL/vector.h>

#include "lzz_pX_CRT.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* checks geometric evaluation at many points                 */
/*------------------------------------------------------------*/
void check(long p)
{
    if (p == 0)
        zz_p::FFTInit(0);
    else
        zz_p::init(p);

    for (long j = 1; j < 100; j+= 5)
    {
        Vec<zz_p> valF;
        zz_pX f;
        zz_pX_Multipoint_Geometric ev;
        zz_p a, b;
        
        f = random_zz_pX(j);
        a = random_zz_p();
        b = random_zz_p();
        long maxd = 5*j;
        valF.SetLength(maxd);
        for (long k = 0; k < maxd; k++)
            valF[k] = eval(f, b * power(a, 2*k));
        ev = zz_pX_Multipoint_Geometric(a, b, j);

        for (long t = 1; t < maxd; t+= 7)
        {
            Vec<zz_p> valF2;
            ev.evaluate(valF2, f, t);
            for (long k = 0; k < t; k++)
                if (valF2[k] != valF[k])
                    Error("bad value in geometric evaluate long");
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
