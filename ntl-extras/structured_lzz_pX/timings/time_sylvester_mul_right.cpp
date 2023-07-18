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

    for (long i = 20; i < 100; i += 10)
    {
        long j = i;
        long d = i;

        Vec<zz_pX> F, G;
        do
            F = random_vec_zz_pX(i, d);
        while (F[i - 1] == 0);
        do
            G = random_vec_zz_pX(j, d);
        while (G[j - 1] == 0);
        
        sylvester_lzz_pX S(F, G);
        Mat<zz_pX> M = S.to_dense();
        
        double u1, u2;
        
        Mat<zz_pX> in, out, out2;
        
        in = random_mat_zz_pX(S.NumCols(), 2, d);
        u1 = get_time();
        out = S.mul_right(in);
        u1 = get_time() - u1;
        u2 = get_time();
        out2 = M * in;
        u2 = get_time() - u2;
        assert (out == out2);
        cout << i << " " << u1 << " " << u2;
        cout << endl;
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
