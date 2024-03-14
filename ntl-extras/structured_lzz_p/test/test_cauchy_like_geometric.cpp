#include <NTL/vec_lzz_p.h>
#include <assert.h>

#include "lzz_p_extra.h"
#include "vec_lzz_p_extra.h"
#include "structured_lzz_p.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* creates cauchy-like matrices                               */
/*------------------------------------------------------------*/
void check(long p)
{
    if (p == 0)
        zz_p::FFTInit(0);
    else
        zz_p::init(p);

    for (long i = 1; i < 1000; i += 1)
    {
        zz_p a = element_of_order(2*i);
        long j = i+10;
        long alpha = 4;
        Mat<zz_p> A, B;
        random(A, i, alpha);
        random(B, j, alpha);
        cauchy_like_geometric_lzz_p M(A, B, to_zz_p(1), power(a, i), a);
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
