#include <NTL/lzz_pX.h>
#include <NTL/matrix.h>
#include <NTL/vector.h>
#include <iomanip>

#include "util.h"
#include "mat_lzz_pX_extra.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
void one_check(long sz, long deg)
{
    Mat<zz_pX> a, x, residue;
    Mat<zz_p> a0;

    do
    {
        random_mat_zz_pX(a, sz, sz, deg);
        GetCoeff(a0, a, 0);
    }
    while (determinant(a0) == 0);
    plain_inv_trunc(x, a, 2*deg);
    mul_trunc(residue, x, a, 2*deg);
}

/*------------------------------------------------------------*/
/* for a give prime, checks some (size, degree)               */
/*------------------------------------------------------------*/
void all_checks()
{
    std::vector<long> szs =
    {
        1, 2, 3, 5, 10, 20, 30
    };

    std::vector<long> degs =
    {
        20, 500
    };

    for (size_t si = 0; si < szs.size(); si++)
        for (size_t di = 0; di < degs.size(); di++)
            one_check(szs[si], degs[di]);

}



/*------------------------------------------------------------*/
/* checks some primes                                         */
/*------------------------------------------------------------*/
void check()
{
    zz_p::FFTInit(0);
    all_checks();
    zz_p::UserFFTInit(786433);
    all_checks();
    zz_p::init(288230376151711813);
    all_checks();
    zz_p::init(786433);
    all_checks();
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
