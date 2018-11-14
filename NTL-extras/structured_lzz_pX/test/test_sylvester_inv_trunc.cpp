#include <NTL/vec_lzz_p.h>
#include <assert.h>

#include "util.h"
#include "vec_lzz_p_extra.h"
#include "structured_lzz_p.h"
#include "mat_lzz_pX_extra.h"
#include "structured_lzz_pX.h"

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

    for (long i = 2; i < 100; i += 7)
    {
        for (long j = 2; j < 100; j += 9)
        {
            for (long d = 7; (d < 1000) && (d*i*j < 70000); d += 15)
            {
                Vec<zz_pX> F, G;
                do
                    F = random_vec_zz_pX(i, d);
                while (F[i - 1] == 0);
                do
                    G = random_vec_zz_pX(j, d + 1);
                while (G[j - 1] == 0);

                sylvester_lzz_pX S(F, G);
                toeplitz_like_minus_lzz_pX iS;
                S.newton_inv_trunc(iS, 19);
                S.newton_inv_trunc(iS, 18);
                exit(0);
                // Mat<zz_pX> M = S.to_dense();

                // Mat<zz_pX> in, out, out2;

                // in = random_mat_zz_pX(S.NumCols(), 2, d + 5);
                // for (long k = 0; k < S.NumCols(); k++)
                //     in[k][1] = random_zz_pX(d);
                // out = S.mul_right(in);
                // multiply(out2, M, in, 0);
                // assert (out == out2);

                // in = random_mat_zz_pX(S.NumCols(), 2, d);
                // for (long k = 0; k < S.NumCols(); k++)
                //     in[k][1] = random_zz_pX(2);
                // out = S.mul_right(in);
                // multiply(out2, M, in, 0);
                // assert (out == out2);
            }
        }
    }
}

/*------------------------------------------------------------*/
/* main just calls check()                                    */
/*------------------------------------------------------------*/
int main(int argc, char** argv)
{
    check(786433);
    check(0);
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
