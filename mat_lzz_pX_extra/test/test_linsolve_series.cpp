#include <NTL/lzz_pX.h>
#include <NTL/matrix.h>
#include <NTL/vector.h>
#include <iomanip>

#include "util.h"
#include "mat_lzz_pX_utils.h"
#include "mat_lzz_pX_multiply.h"
#include "mat_lzz_pX_linsolve.h"

NTL_CLIENT


/*------------------------------------------------------------*/
/* checks an (sz,sz) matrix in degree < deg                   */
/*------------------------------------------------------------*/
void one_check(long sz, long len)
{
    Mat<zz_pX> A, b, u, res;

    long deg = len-1;
    long hlen = 5*len;

    do
        random(A, sz, sz, len);
    while (determinant(coeff(A,0)) == 0);

    random(b, sz, 1, len);
    Mat<zz_pX> truncB = trunc(b,deg); // for high-order lifting

    solve_series_low_precision(u, A, b, len);
    multiply(res, A, u);
    res = trunc(res - b, len);
    if (!IsZero(res))
        LogicError("Bad output for low precision series solve");

    solve_series_high_precision(u, A, b, len);
    multiply(res, A, u);
    res = trunc(res - b, len);
    if (!IsZero(res))
        LogicError("Bad output for low precision series solve");

    solve_series_high_order_lifting(u, A, truncB, len);
    multiply(res, A, u);
    res = trunc(res - truncB, len);
    if (!IsZero(res))
        LogicError("Bad output for high-order-lifting series solve");

    random(b, sz, 1, hlen);
    trunc(truncB, b, deg); // for high-order lifting

    solve_series_low_precision(u, A, b, hlen);
    multiply(res, A, u);
    res = trunc(res - b, hlen);
    if (!IsZero(res))
        LogicError("Bad output for high precision series solve");

    solve_series_high_precision(u, A, b, hlen);
    multiply(res, A, u);
    res = trunc(res - b, hlen);
    if (!IsZero(res))
        LogicError("Bad output for high precision series solve");

    solve_series_high_order_lifting(u, A, truncB, hlen);
    multiply(res, A, u);
    res = trunc(res - truncB, hlen);
    if (!IsZero(res))
        LogicError("Bad output for high-order-lifting series solve");
}

/*------------------------------------------------------------*/
/* for a give prime, checks some (size, degree)               */
/*------------------------------------------------------------*/
void all_checks()
{
    VecLong szs =
    {
        1, 2, 3, 5, 10, 20, 30
    };

    VecLong degs =
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
