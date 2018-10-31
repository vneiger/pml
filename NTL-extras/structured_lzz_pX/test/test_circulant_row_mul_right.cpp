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

    for (long i = 1; i < 300; i += (i < 10 ? 1 : 19))
        for (long k = 1; k < 300; k += (k < 10 ? 1 : 21))
        {
            Vec<zz_pX> dat;
            dat.SetLength(i);
            for (long d = 1; d < 300 &&  d*d*i*k < 4000000; d += (d < 100 ? 1 : 19))
            {
                circulant_row_lzz_pX c;
                Vec<zz_pX> in, out, out2;
                Mat<zz_pX> inM, outM;
                in.SetLength(i);

                for (long j = 0; j < i; j++)
                    dat[j] = random_zz_pX(d);
                for (long j = 0; j < i; j++)
                    in[j] = random_zz_pX(d);

                c = circulant_row_lzz_pX(dat, k);

                out = c.mul_right(in);

                inM.SetDims(i, 1);
                for (long j = 0; j < i; j++)
                    inM[j][0] = in[j];
                outM = c.to_dense() * inM;
                out2.SetLength(k);
                for (long j = 0; j < k; j++)
                    out2[j] = outM[j][0];
                assert (out == out2);

                c.prepare_degree(max(1, d/2));
                out2 = c.mul_right(in);
                assert (out == out2);

                c.prepare_degree(d-1);
                out2 = c.mul_right(in);
                assert (out == out2);

                for (long j = 0; j < i/2; j++)
                    dat[j] = 0;
                
                for (long j = i/3; j < i; j++)
                    in[j] = 0;

                c = circulant_row_lzz_pX(dat, k);

                out = c.mul_right(in);

                inM.SetDims(i, 1);
                for (long j = 0; j < i; j++)
                    inM[j][0] = in[j];
                outM = c.to_dense() * inM;
                out2.SetLength(k);
                for (long j = 0; j < k; j++)
                    out2[j] = outM[j][0];
                assert (out == out2);

                c.prepare_degree(max(1, d/2));
                out2 = c.mul_right(in);
                assert (out == out2);

                c.prepare_degree(d-1);
                out2 = c.mul_right(in);
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
