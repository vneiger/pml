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
        Vec<zz_p> A, invA1, invA2;
        random_vec_zz_p(A, i);
     
        inv(invA1, A);
        inv_naive(invA2, A);
        for (long j = 0; j < i; j++)
        {
            assert (invA1[j] == invA2[j]);
            assert (invA1[j] == 1/A[j]);
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
