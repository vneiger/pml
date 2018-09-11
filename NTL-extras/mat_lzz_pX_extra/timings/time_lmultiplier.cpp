#include <NTL/lzz_pX.h>
#include <NTL/matrix.h>
#include <NTL/vector.h>
#include <iomanip>

#include "util.h"
#include "mat_lzz_pX_extra.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* checks a product (s,s) x (s,1) in degree < deg             */
/* checks a product (s,s) x (s,s) in degree < deg             */
/*------------------------------------------------------------*/
void one_check(long sz, long deg)
{
    Mat<zz_pX> a, b, c1, c2;
    mat_lzz_pX_lmultiplier_FFT mula;
    double t_plain, t_lmul;
    long nb;
    const double thres = 0.01;

    cout << sz << " " << deg << " ";

    random_mat_zz_pX(a, sz, sz, deg);
    mula = mat_lzz_pX_lmultiplier_FFT(a, deg-1);

    random_mat_zz_pX(b, sz, 1, deg);
    
    t_plain = get_time();
    nb = 0;
    do
    {
        multiply(c1, a, b);
        nb++;
    }
    while ((get_time()-t_plain) <= thres);
    t_plain = (get_time()-t_plain) / nb;

    t_lmul = get_time();
    nb = 0;
    do
    {
        mula.multiply(c2, b);
        nb++;
    }
    while ((get_time()-t_lmul) <= thres);
    t_lmul = (get_time()-t_lmul) / nb;

    cout << t_plain << " " << t_lmul << " ";


    random_mat_zz_pX(b, sz, sz, deg);
    
    t_plain = get_time();
    nb = 0;
    do
    {
        multiply(c1, a, b);
        nb++;
    }
    while ((get_time()-t_plain) <= thres);
    t_plain = (get_time()-t_plain) / nb;

    t_lmul = get_time();
    nb = 0;
    do
    {
        mula.multiply(c2, b);
        nb++;
    }
    while ((get_time()-t_lmul) <= thres);
    t_lmul = (get_time()-t_lmul) / nb;

    cout << t_plain << " " << t_lmul << " ";

    cout << endl;
}

/*------------------------------------------------------------*/
/* for a give prime, checks some (size, degree)               */
/*------------------------------------------------------------*/
void all_checks()
{
    std::vector<long> szs =
    {
        1, 2, 3, 5, 10, 20, 30, 50, 100, 150, 200
    };

    std::vector<long> degs =
    {
        1, 2, 3, 4, 5, 10, 15, 20, 25, 50, 60, 70, 100, 150, 200
    };

    cout << "p=" << zz_p::modulus() << " FFT=" << is_FFT_prime() << endl;
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
    // zz_p::init(288230376151711813);
    // all_checks();
    // zz_p::init(786433);
    // all_checks();
}  

/*------------------------------------------------------------*/
/* main calls check                                           */
/*------------------------------------------------------------*/
int main(int argc, char ** argv)
{
    std::cout << std::fixed;
    std::cout << std::setprecision(8);
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
