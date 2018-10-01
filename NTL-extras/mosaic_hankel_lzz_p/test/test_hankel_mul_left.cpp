#include <NTL/vec_lzz_p.h>
#include <assert.h>

#include "vec_lzz_p_extra.h"
#include "mosaic_hankel_lzz_p.h"

NTL_CLIENT


/*------------------------------------------------------------*/
/* left multiplies by a vector and 2-column matrix            */
/*------------------------------------------------------------*/
void check(long p)
{
    if (p == 0)
        zz_p::FFTInit(0);
    else
        zz_p::init(p);

    for (long i1 = 1; i1 < 200; i1 += 1)
        for (long i2 = 1; i2 < 200; i2 += 1)
        {
            Vec<zz_p> dat, input, output, output2;
            hankel_lzz_p h;
            Mat<zz_p> M;

            dat = random_vec_zz_p(i1+i2-1);
            h = hankel_lzz_p(dat, i1, i2);
            M = h.to_dense();

            input = random_vec_zz_p(i1);
            output = random_vec_zz_p(i2);
            output = h.mul_left(input);
            output2 = input*M;
            
            assert (output == output2);
            cout << i1 << " " << i2 << endl;
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
