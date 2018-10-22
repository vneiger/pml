#include <NTL/lzz_pX.h>
#include <NTL/vector.h>
#include <assert.h>

#include "lzz_pX_CRT.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* tests geometric evaluation, then interpolation             */
/*------------------------------------------------------------*/
void check(long p)
{
    if (p == 0)
        zz_p::FFTInit(0);
    else
        zz_p::init(p);

    for (long j = 1; j < 200; j++)
    {
        zz_p a, b;
        zz_pX_Multipoint_Geometric ev;
        Vec<zz_p> in, out, out2;
        Mat<zz_p> M;

        a = random_zz_p();
        b = random_zz_p();
        ev = zz_pX_Multipoint_Geometric(a, b, j);
        M = ev.to_dense();
        in = random_vec_zz_p(j);
        
        ev.mul_right(out, in);
        out2 = M * in;
        assert (out2 == out);

        ev.inv_mul_right(out2, out);
        assert (out2 == in);

        ev.mul_left(out, in);
        out2 = in * M;
        assert (out2 == out);

        ev.inv_mul_left(out2, out);
        assert (out2 == in);
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
