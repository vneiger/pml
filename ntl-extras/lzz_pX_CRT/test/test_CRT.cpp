#include <NTL/lzz_pX.h>
#include <NTL/vector.h>

#include "lzz_pX_CRT.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* does a CRT + multimod                                      */
/*------------------------------------------------------------*/
void check()
{

    long p = 1125899906842679;
    zz_p::init(p);

    for (long j = 1; j < 50; j++)
    {
        Vec<zz_pX> val, q, mmod;
        q.SetLength(j);
        val.SetLength(j);
        for (long i = 0; i < j; i++)
        {
            q[i] = random_zz_pX(i*3+5);
            val[i] = random_zz_pX(i*3+5) % q[i];
        }
        zz_pX_CRT crt = zz_pX_CRT(q);
        zz_pX f;
        crt.combine(f, val);
        crt.multimod(mmod, f);
        for (long i = 0; i < j; i++)
        {
            if (val[i] != (f % q[i]))
            {
                cerr << "error for CRT combine with j=" << j << endl;
                exit(-1);
            }
            if (val[i] != mmod[i])
            {
                cerr << "error for CRT multimod with j=" << j << endl;
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
