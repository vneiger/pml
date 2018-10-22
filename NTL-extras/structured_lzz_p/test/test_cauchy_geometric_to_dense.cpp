#include <NTL/vec_lzz_p.h>
#include <assert.h>

#include "structured_lzz_p.h"

NTL_CLIENT


/*------------------------------------------------------------*/
/* creates cauchy matrices, makes them dense                  */
/*------------------------------------------------------------*/
void check(long p)
{
    if (p == 0)
        zz_p::FFTInit(0);
    else
        zz_p::init(p);

    for (long i = 1; i < 100; i += 1)
    {
        zz_p a = random_zz_p();
        cauchy_geometric_lzz_p C;
        Mat<zz_p> M;
        C = cauchy_geometric_lzz_p(to_zz_p(1), power(a, i), a, i, i+19);
        M = C.to_dense();
        C = cauchy_geometric_lzz_p(to_zz_p(1), power(a, i+19), a, i+19, i);
        M = C.to_dense();
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
