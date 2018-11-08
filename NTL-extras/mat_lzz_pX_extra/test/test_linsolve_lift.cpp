#include <NTL/lzz_pX.h>
#include <NTL/matrix.h>
#include <NTL/vector.h>
#include <iomanip>

#include "util.h"
#include "mat_lzz_pX_extra.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* checks an (sz,sz) matrix in degree < deg                   */
/*------------------------------------------------------------*/
void one_check(long sz, long deg)
{
    Mat<zz_pX> A;
    Vec<zz_pX> b, b2, u, res;
    zz_pX den;

    random(A, sz, sz, deg);
    random(b, sz, deg);

    linsolve_via_series(u, den, A, b);
    b2.SetLength(sz);
    for (long i = 0; i < sz; i++)
        b2[i] = b[i] * den;
    multiply(res, A, u);
    res = res - b2;
    if (!IsZero(res))
        Error("non-zero residue in linsolve via series");
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
        1, 2, 3, 4, 5, 10, 15, 20, 25, 50, 60, 70, 100, 150, 200, 250, 300, 400
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
