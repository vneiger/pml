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
