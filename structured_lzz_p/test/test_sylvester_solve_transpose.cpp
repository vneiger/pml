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
        for (long j = 2; j < 100; j += 3)
        {
            zz_pX F, G;

            do{
                do
                    F = random_zz_pX(i);
                while (deg(F) < 1);
                
                do
                    G = random_zz_pX(j);
                while (deg(G) < 1);
            }
            while (deg(GCD(F, G)) != 0);

            sylvester_lzz_p S;
            long n = deg(F) + deg(G);
            S = sylvester_lzz_p(F, G);

            Mat<zz_p> in = random_mat_zz_p(3, n), in2 = random_mat_zz_p(3, n);
            Mat<zz_p> out;
            out = S.mul_left(in);
            long r;

            r = S.solve_transpose(in2, out);
            assert (r == 1);
            assert (in == in2);

            S.prepare_inverse();

            r = S.solve_transpose(in2, out);
            assert (r == 1);
            assert (in == in2);
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
