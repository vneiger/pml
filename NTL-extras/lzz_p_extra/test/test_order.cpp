#include <NTL/lzz_pX.h>
#include <NTL/vector.h>

#include "lzz_p_extra.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* computes a few orders                                      */
/*------------------------------------------------------------*/
void check()
{
    long p = 65537;
    zz_p::init(p);
    zz_p a;

    a = to_zz_p(-1);
    order(a);

    a = random_zz_p();
    order(a);
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
