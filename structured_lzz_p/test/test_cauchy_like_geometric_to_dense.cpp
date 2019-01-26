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

    for (long i = 1; i < 100; i += 1)
    {
        zz_p a = element_of_order(2*i);
        long j = i+10;
        long alpha = 4;
        Mat<zz_p> A, B;
        random(A, i, alpha);
        random(B, j, alpha);
        cauchy_like_geometric_lzz_p M(A, B, to_zz_p(1), power(a, i), a);

        Mat<zz_p> MD = M.to_dense();
        mat_zz_p D1, D2;
        D1.SetDims(i, i);
        for (long k = 0; k < i; k++)
        {
            D1[k][k] = power(a, k);
        }
        D2.SetDims(j, j);
        for (long k = 0; k < j; k++)
        {
            D2[k][k] = power(a, i)*power(a, k);
        }
        assert ((D1*MD - MD*D2) == (A*transpose(B)));
        
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
