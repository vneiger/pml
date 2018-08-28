#include <NTL/lzz_pX.h>
#include <NTL/vector.h>

#include "sage_output.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* prints a vector and a polynomial                           */
/*------------------------------------------------------------*/
void check()
{
    long p = 1125899906842679;
    zz_p::init(p);

    sage_init();
    sage_init_X();

    Vec<zz_p> v;
    v.SetLength(10);
    for (long i = 0; i < v.length(); i++)
        v[i] = random_zz_p();
    sage_assign(v, "v");

    zz_pX f = random_zz_pX(10);
    sage_assign(f, "f");
    sage_output(f, "x");
}  

/*------------------------------------------------------------*/
/* main calls check                                           */
/*------------------------------------------------------------*/
int main(int argc, char ** argv)
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
