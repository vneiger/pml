#include <NTL/vec_lzz_p.h>
#include <assert.h>

#include "structured_lzz_pX.h"

PML_CLIENT

/*------------------------------------------------------------*/
/* creates toeplitz-like matrices                             */
/* does left truncated products                               */
/*------------------------------------------------------------*/
void check(long p)
{
    if (p == 0)
        zz_p::FFTInit(0);
    else
        zz_p::init(p);

    for (long i = 2; i < 100; i += (i < 10 ? 5 : 19))
        for (long k = 1; k < 300; k += (k < 10 ? 7 : 21))
        {
            Vec<zz_pX> dat;
            dat.SetLength(i);
            for (long d = 2; d < 300 &&  (d+10)*(i+k+5) < 7000; d += (d < 100 ? 11 : 19))
            {
                Mat<zz_pX> G = random_mat_zz_pX(i, 2, d);
                Mat<zz_pX> H = random_mat_zz_pX(k, 2, d);
                toeplitz_like_minus_lzz_pX T(G, H);
                Mat<zz_pX> M = T.to_dense();
                trunc(M, M, d/2);
                Mat<zz_pX> in, out, out2;

                in = random_mat_zz_pX(2, i, d);
                out = T.mul_left_trunc(in, d/2);
                
                out2 = in * M ;
                assert (out == trunc(out2, d/2));
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
