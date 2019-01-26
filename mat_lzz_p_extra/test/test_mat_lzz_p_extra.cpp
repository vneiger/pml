#include <NTL/vec_lzz_p.h>
#include <assert.h>

#include "util.h"
#include "vec_lzz_p_extra.h"
#include "mat_lzz_p_extra.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* creates hankel matrices                                    */
/*------------------------------------------------------------*/
void check(long p)
{
    if (p == 0)
        zz_p::FFTInit(0);
    else
        zz_p::init(p);
    
    for (long i = 1; i < 100; i += 10)
    {
        Mat<zz_p> z = Z_lzz_p(i, random_zz_p());
        if (i <= 20 && p < (1L << 30) && p != 0)
            cout << "(" << i << " x " << i << ") Z:\n" << z << endl;

        Mat<zz_p> j = J_lzz_p(i);
        if (i <= 20 && p < (1L << 30) && p != 0)
            cout << "(" << i << " x " << i << ") J:\n" << j << endl;
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
