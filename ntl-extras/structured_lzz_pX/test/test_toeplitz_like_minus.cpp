#include <NTL/vec_lzz_p.h>
#include <assert.h>

#include "structured_lzz_pX.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* creates toeplitz-like matrices                             */
/*------------------------------------------------------------*/
void check(long p)
{
    if (p == 0)
        zz_p::FFTInit(0);
    else
        zz_p::init(p);

    for (long i = 2; i < 300; i += (i < 10 ? 1 : 19))
        for (long k = 1; k < 300; k += (k < 10 ? 1 : 21))
        {
            Vec<zz_pX> dat;
            dat.SetLength(i);
            for (long d = 2; d < 300 &&  d*d*i*k < 4000000; d += (d < 100 ? 11 : 19))
            {
                Mat<zz_pX> G = random_mat_zz_pX(i, 2, d);
                Mat<zz_pX> H = random_mat_zz_pX(k, 2, d);
                toeplitz_like_minus_lzz_pX T(G, H);
            }
        }
}

/*------------------------------------------------------------*/
/* main just calls check()                                    */
/*------------------------------------------------------------*/
int main(int argc, char** argv)
{
    check(0);
    check(786433);
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
