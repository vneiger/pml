#include <NTL/vec_lzz_p.h>
#include <assert.h>

#include "structured_lzz_pX.h"

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


    for (long i = 1; i < 300; i += (i < 100 ? 1 : 17))
    {
        Vec<zz_pX> dat;
        dat.SetLength(i);
        for (long d = 0; d < 300 &&  d*i < 2000; d += (d < 100 ? 1 : 19))
        {
            lower_triangular_toeplitz_lzz_pX h;
            Mat<zz_pX> in, out, out2;

            for (long j = 0; j < i; j++)
                dat[j] = random_zz_pX(d);

            in = random_mat_zz_pX(2, i, d);
            h = lower_triangular_toeplitz_lzz_pX(dat);
            out = h.mul_left_trunc(in, d/2);
            out2 = in * h.to_dense();
            out2 = trunc(out2, d/2);
            assert (out == out2);

            for (long j = 0; j < i/3; j++)
                dat[j] = 0;
            for (long j = 0; j < i/2; j++)
                in[0][j] = 0;

            h = lower_triangular_toeplitz_lzz_pX(dat);
            out = h.mul_left_trunc(in, d/2);
            out2 = in * h.to_dense();
            out2 = trunc(out2, d/2);
            assert (out == out2);
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
