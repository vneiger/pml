#include <NTL/vec_lzz_p.h>
#include <assert.h>

#include "vec_lzz_p_extra.h"
#include "structured_lzz_p.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* creates sylvester matrices                                 */
/*------------------------------------------------------------*/
void check(long p)
{
    if (p == 0)
        zz_p::FFTInit(0);
    else
        zz_p::init(p);

    for (long i = 2; i < 100; i += 3)
    {
        for (long j = 3; j < 100; j += 3)
        {
            zz_pX F, G;
            do
                F = random_zz_pX(i);
            while (deg(F) < 1);
            
            do
                G = random_zz_pX(j);
            while (deg(G) < 1);
            sylvester_lzz_p S(F, G);
            Mat<zz_p> M = S.to_dense();
            
            Mat<zz_p> U, V;
            S.phi_plus_generators(U, V);
            assert (toeplitz_lzz_p_phi_plus(M) == (U * transpose(V)));
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
