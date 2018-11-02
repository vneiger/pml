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
            Vec<zz_pX> in, out, out2;
            Mat<zz_pX> inM, outM;

            in.SetLength(i);
            out.SetLength(i);

            for (long j = 0; j < i; j++)
                dat[j] = random_zz_pX(d);
            for (long j = 0; j < i; j++)
                in[j] = random_zz_pX(d);


            h = lower_triangular_toeplitz_lzz_pX(dat);

            out = h.mul_right(in);
            
            inM.SetDims(i, 1);
            for (long j = 0; j < i; j++)
                inM[j][0] = in[j];
            outM = h.to_dense() * inM;
            out2.SetLength(i);
            for (long j = 0; j < i; j++)
                out2[j] = outM[j][0];
            assert (out == out2);

            h.prepare_degree(max(1, d/2));
            out2 = h.mul_right(in);
            assert (out == out2);

            h.prepare_degree(d-1);
            out2 = h.mul_right(in);
            assert (out == out2);

            for (long j = 0; j < i/3; j++)
                dat[j] = 0;
            for (long j = 0; j < i/2; j++)
                in[j] = 0;

            h = lower_triangular_toeplitz_lzz_pX(dat);

            out = h.mul_right(in);
            
            inM.SetDims(i, 1);
            for (long j = 0; j < i; j++)
                inM[j][0] = in[j];
            outM = h.to_dense() * inM;
            out2.SetLength(i);
            for (long j = 0; j < i; j++)
                out2[j] = outM[j][0];
            assert (out == out2);

            h.prepare_degree(max(1, d/2));
            out2 = h.mul_right(in);
            assert (out == out2);

            h.prepare_degree(d-1);
            out2 = h.mul_right(in);
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
