#include <NTL/vec_lzz_p.h>
#include <assert.h>

#include "vec_lzz_p_extra.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* runs a check up to size 1000                               */
/*------------------------------------------------------------*/
void check()
{
    zz_p::FFTInit(0);

    for (long i = 0; i < 1000; i += 1)
    {
        Vec<zz_p> A, oldA, invA2;
        random(A, i);
        oldA = A;

        inv_naive(invA2, A);
        inv(A, A);
        for (long j = 0; j < i; j++)
        {
            assert (A[j] == invA2[j]);
            assert (A[j] == 1/oldA[j]);
        }
    }
}

/*------------------------------------------------------------*/
/* main just calls check()                                    */
/*------------------------------------------------------------*/
int main(int argc, char** argv)
{
    check();
    return 0;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
