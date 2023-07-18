#include <NTL/lzz_pX.h>
#include <NTL/vector.h>

#include "lzz_p_extra.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* tries to find an element of order p-1                      */
/*------------------------------------------------------------*/
void check(long p)
{
    long ord;
    zz_p a;
    zz_p::init(p);

    ord = p - 1;
    element_of_order(a, ord);
    if (a != 0)
        for (long i = 1; i < ord; i++)
            1/(power(a, i)-1);
}  

/*------------------------------------------------------------*/
/* main just calls check()                                    */
/*------------------------------------------------------------*/
int main(int argc, char ** argv)
{
    check(1024);
    check(9001);
    return 0;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
