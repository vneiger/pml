#include <NTL/vec_lzz_p.h>
#include <assert.h>

#include "structured_lzz_pX.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* creates circulant matrices and right multiplies            */
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
                circulant_row_lzz_pX c;
                Mat<zz_pX> in, out, out2;

                for (long j = 0; j < i; j++)
                    dat[j] = random_zz_pX(d);
                in = random_mat_zz_pX(2, k, d);
                c = circulant_row_lzz_pX(dat, k);
                out = c.mul_left_trunc(in, d/2);
                out2 = in * c.to_dense();
                assert (out == trunc(out2, d/2));

                for (long j = 0; j < i/2; j++)
                    dat[j] = 0;
                for (long j = i/3; j < k; j++)
                    in[0][j] = 0;

                c = circulant_row_lzz_pX(dat, k);
                out = c.mul_left_trunc(in, d/2);
                out2 = in * c.to_dense();
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
