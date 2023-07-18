#include <NTL/vec_lzz_p.h>
#include <assert.h>

#include "vec_lzz_p_extra.h"
#include "structured_lzz_p.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* creates lower triangular toeplitz matrices                 */
/*------------------------------------------------------------*/
void check(long p)
{
    if (p == 0)
        zz_p::FFTInit(0);
    else
        zz_p::init(p);

    for (long i = 1; i < 300; i += 1)
    {
        Vec<zz_p> dat;
        dat = random_vec_zz_p(i);
        lower_triangular_toeplitz_lzz_p h = lower_triangular_toeplitz_lzz_p(dat);
        Vec<zz_p> in, out1, out2;
        in = random_vec_zz_p(i);
        out1 = h.mul_left(in);
        out2 = in * h.to_dense();
        assert (out1 == out2);
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
