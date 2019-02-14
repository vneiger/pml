#include <NTL/lzz_pX.h>
#include <NTL/matrix.h>
#include <NTL/vector.h>
#include <iomanip>

#include "util.h"
#include "mat_lzz_pX_extra.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
void check_solution(long sz, long deg)
{
    double t;
    long nb;
    const double thresh = 0.2;
    Mat<zz_pX> a, b, sol;
    Mat<zz_p> a0;

    do
    {
        random(a, sz, sz, deg+1);
        GetCoeff(a0, a, 0);
    }
    while (determinant(a0) == 0);

    cout << sz << " " << deg << " ";

    random(b, sz, 1, deg);

    t = get_time();
    nb = 0;
    do
    {
        solve_series_high_order_lifting(sol, a, b, sz*deg);
        nb++;
    }
    while ((get_time()-t) <= thresh);
    t = (get_time()-t) / nb;
    cout << t << " ";

    Mat<zz_pX> check;

    t = get_time();
    nb = 0;
    do
    {
        check = solve_series_high_precision(a, b, sz*deg);
        nb++;
    }
    while ((get_time()-t) <= thresh);
    t = (get_time()-t) / nb;
    cout << t << " ";

    cout << endl;
}


/*------------------------------------------------------------*/
/* for a give prime, checks some (size, degree)               */
/*------------------------------------------------------------*/
void all_checks()
{

    VecLong szs =
    {
        32
    };

    VecLong degs =
    {
        1024
    };
    // VecLong szs =
    // {
    //     1, 2, 3, 5, 10, 20, 30, 50, 100
    // };

    // VecLong degs =
    // {
    //     20, 50, 75, 100, 150, 250, 350
    // };

    for (size_t si = 0; si < szs.size(); si++)
        for (size_t di = 0; di < degs.size(); di++)
            check_solution(szs[si], degs[di]);

}



/*------------------------------------------------------------*/
/* checks some primes                                         */
/*------------------------------------------------------------*/
void check()
{
    //zz_p::FFTInit(0);
    //all_checks();
    // zz_p::UserFFTInit(786433);
    // all_checks();
    zz_p::init(288230376151711813);
    all_checks();
    // zz_p::init(786433);
    // all_checks();
}  

/*------------------------------------------------------------*/
/* main calls check                                           */
/*------------------------------------------------------------*/
int main(int argc, char ** argv)
{
    warmup();
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
