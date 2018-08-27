#include <NTL/lzz_pX.h>
#include <NTL/vector.h>

#include "lzz_pX_CRT.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* displays (to cerr) thresholds for geometric direct/FFT     */
/*------------------------------------------------------------*/
void check(long p)
{
    if (p == 0)
        zz_p::FFTInit(0);
    else
        zz_p::init(p);

    cerr << p << " " << min_geometric_FFT() << " " << max_geometric_FFT() << endl;
}  

/*------------------------------------------------------------*/
/* main calls check                                           */
/*------------------------------------------------------------*/
int main(int argc, char ** argv)
{
    check(0);
    check(23068673);
    check(288230376151711813);
    return 0;
}
