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

    for (long i = 1; i < 40; i += 1)
        for (long k = 1; k < 40; k += 1)
        {
            Vec<zz_pX> dat;
            dat.SetLength(i);
            for (long d = 1; d < 100; d += (d < 20 ? 1 : 17))
            {
                for (long j = 0; j < i; j++)
                    dat[j] = random_zz_pX(d);
                circulant_row_lzz_pX c(dat, k);
                if (p != 0 && p < (1L << 23) && (i * (d+1)) < 15 && i * k < 13 )
                {
                    cout << "(" << i << " x " << k << ") circulant:\n" << c.to_dense() << endl;
                    cout << dat << endl;
                }
                
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
