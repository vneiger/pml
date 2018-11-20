#include <NTL/vec_lzz_p.h>
#include <assert.h>

#include "util.h"
#include "vec_lzz_p_extra.h"
#include "lzz_pX_extra.h"
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

    for (long i = 500; i < 600; i += 700)
    {
        long d = i;
        long j = i;
        Vec<zz_pX> F, G;
        do
            F = random_vec_zz_pX(i, d);
        while (F[i - 1] == 0);
        do
            G = random_vec_zz_pX(j, d);
        while (G[j - 1] == 0);
        
        sylvester_lzz_pX S(F, G);
        toeplitz_like_minus_lzz_pX iS;
        double t;
        cout << i << " " << j << " " << d << " ";

        long factor = 2;

        zz_pX a0, b0;
        for (long i = 0; i < F.length(); i++)
            SetCoeff(a0, i, coeff(F[i], 0));
        for (long i = 0; i < G.length(); i++)
            SetCoeff(b0, i, coeff(G[i], 0));


        t = get_time();
        for (long nb = 0; nb < factor * d; nb++)
        {
            Mat<zz_p> block;
            sylvester_lzz_p S(a0, b0);
            long r = S.top_right_block_inverse(block, (long) cbrt(d));
        }
        cout << get_time()-t << " ";

        t = get_time();
        for (long nb = 0; nb < factor * d; nb++)
            zz_p r = resultant(a0, b0);
        cout << get_time()-t << " ";

        // t = get_time();
        // S.high_precision_inv_trunc(iS, factor * d);
        // cout << get_time()-t << endl;

        cout << endl;
        if (i < 100)
        {
            Mat<zz_pX> M = trunc(S.to_dense(), 10 * d);
            Mat<zz_pX> iM = trunc(iS.to_dense(), 10 * d);
            Mat<zz_pX> Id = trunc(M * iM, 10 * d);
            assert (IsIdent(Id, M.NumRows()));
        }
    }
}

/*------------------------------------------------------------*/
/* main just calls check()                                    */
/*------------------------------------------------------------*/
int main(int argc, char** argv)
{
    // warmup();
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
