#include <NTL/lzz_pX.h>
#include <NTL/matrix.h>
#include <NTL/vector.h>

#include "mat_lzz_pX_extra.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* creates random polynomial matrix                           */
/* check takes an extra argument, not used here               */
/*------------------------------------------------------------*/
void check()
{
    zz_p::init(1125899906842679);
    Mat<zz_pX> a;
    random(a, 3, 4, 2);
    cout << a << endl;

}  

int main()
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
