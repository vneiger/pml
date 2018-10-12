#include <NTL/lzz_pX.h>
#include <NTL/matrix.h>
#include <NTL/vector.h>
#include <iomanip>

#include "util.h"
#include "mat_lzz_pX_extra.h"

NTL_CLIENT

static Mat<zz_pX> check(const Mat<zz_pX>& a, long d, long i)
{
    Mat<zz_pX> s = inv_trunc(a, i);
    return s >> (i - (2*d-1));
}

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
void check_inverse(long sz, long deg)
{
    Mat<zz_pX> a, slice;
    Mat<zz_p> a0;
    std::unique_ptr<mat_lzz_pX_lmultiplier> ma, minvA;

    do
    {
        random(a, sz, sz, deg+1);
        GetCoeff(a0, a, 0);
    }
    while (determinant(a0) == 0);

    ma = get_lmultiplier(a, deg-1);
    minvA = get_lmultiplier(inv_trunc(a, deg-1), deg-2);

    slice = inv_trunc(a, 2*deg) >> 1;

    high_order_lift_inverse_odd(slice, slice, ma, minvA, deg);
    if (slice != check(a, deg, 3*deg)) 
        LogicError("Error with high order lift--odd");

    high_order_lift_inverse_odd(slice, slice, ma, minvA, deg);
    if (slice != check(a, deg, 5*deg)) 
        LogicError("Error with high order lift--odd");

    high_order_lift_inverse_odd(slice, slice, ma, minvA, deg);
    if (slice != check(a, deg, 9*deg)) 
        LogicError("Error with high order lift--odd");

}

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
void check_solution(long sz, long deg)
{
    Mat<zz_pX> a, b, sol;
    Mat<zz_p> a0;

    do
    {
        random(a, sz, sz, deg+1);
        GetCoeff(a0, a, 0);
    }
    while (determinant(a0) == 0);

    random(b, sz, 1, deg);
    solve_series_high_order_lifting(sol, a, b, sz*deg);

    Mat<zz_pX> check = solve_series_high_precision(a, b, sz*deg);
    if (check != sol) 
    {
        cout << sz << endl;
        cout << deg << endl;
        LogicError("Error with high order lifting-solution");
    }

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
        20, 50, 75, 100, 150
    };

    for (size_t si = 0; si < szs.size(); si++)
        for (size_t di = 0; di < degs.size(); di++)
        {
            check_inverse(szs[si], degs[di]);
            check_solution(szs[si], degs[di]);
        }
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
