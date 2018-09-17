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
    Mat<zz_pX> A, b, u, res;
    const double thresh = 0.01;
    random_mat_zz_pX(A, sz, sz, deg);
    random_mat_zz_pX(b, sz, 1, deg);

    double t;
    long nb;
    
    cout << sz << " " << deg << " ";

    t = get_time();
    nb = 0;
    do
    {
        solve_series_low_precision(u, A, b, deg);
        nb++;
    }
    while ((get_time()-t) <= thresh);
    t = (get_time()-t) / nb;
    cout << t << " ";

    t = get_time();
    nb = 0;
    do
    {
        solve_series_low_precision2(u, A, b, deg);
        nb++;
    }
    while ((get_time()-t) <= thresh);
    t = (get_time()-t) / nb;
    cout << t << " ";

    t = get_time();
    nb = 0;
    do
    {
        solve_series_low_precision3(u, A, b, deg);
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
    std::vector<long> szs =
    {
        30
        // 1, 2, 3, 5, 10, 20, 30
    };

    std::vector<long> degs =
    {
        200
        // 1, 2, 3, 4, 5, 10, 15, 20, 25, 50, 60, 70, 100, 150, 200, 250, 300, 400
    };

    for (size_t si = 0; si < szs.size(); si++)
        for (size_t di = 0; di < degs.size(); di++)
            one_check(szs[si], degs[di]);
}


/*------------------------------------------------------------*/
/* checks some primes                                         */
/*------------------------------------------------------------*/
void check(long p)
{
    if (p == 0)
        zz_p::FFTInit(0);
    else
        zz_p::init(p);
    all_checks();
}  

/*------------------------------------------------------------*/
/* main calls check                                           */
/*------------------------------------------------------------*/
int main(int argc, char ** argv)
{
    std::cout << std::fixed;
    std::cout << std::setprecision(8);
    warmup();
    check(0);
    check(23068673);
    check(288230376151711813);
    return 0;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
