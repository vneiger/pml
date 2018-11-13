#include <NTL/lzz_pX.h>
#include <NTL/matrix.h>
#include <NTL/vector.h>
#include <iomanip>

#include "util.h"
#include "lzz_p_extra.h"
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
    double t, tmin;
    long idx;

    random(A, sz, sz, deg);
    random(b, sz, deg);
    
    cout << sz << " " << deg << " ";
    tmin = 10000.0;
    idx = -1;
    for (long nb = 1; nb < 6; nb++)
    {
        t = get_time();
        linsolve_via_series(u, den, A, b, nb);
        t = get_time()-t;
        cout << t << " ";
        if (t < tmin)
        {
            tmin = t;
            idx = nb;
        }
    }
    cout << "   " << idx << endl;
}

/*------------------------------------------------------------*/
/* for a give prime, checks some (size, degree)               */
/*------------------------------------------------------------*/
void all_checks()
{
    VecLong szs =
    {
        // 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 175, 200, 250
        40
    };

    VecLong degs =
    {
        // 15, 20, 25, 50, 60, 70, 100, 150, 200, 250, 300, 400
        1000
    };

    cout << zz_p::modulus() << " " << is_FFT_prime() << endl;
    for (size_t si = 0; si < szs.size(); si++)
        for (size_t di = 0; di < degs.size(); di++)
            if (szs[si] * degs[di] < 150000)
                one_check(szs[si], degs[di]);

}


/*------------------------------------------------------------*/
/* checks some primes                                         */
/*------------------------------------------------------------*/
void check()
{
    zz_p::FFTInit(0);
    all_checks();
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
