#include <NTL/BasicThreadPool.h>
#include <NTL/lzz_pX.h>
#include <NTL/matrix.h>
#include <NTL/vector.h>
#include <iomanip>

#include "util.h"
#include "mat_lzz_pX_extra.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* checks a product (s,s) x (s,s) in degree < deg             */
/*------------------------------------------------------------*/
void one_check(long sz, long deg)
{
    const double thresh = 0.001;
    long nb;
    Mat<zz_pX> a, b, c1, c2;
    double t_FFT, t_waksman, t_transform, t_geometric, t_dense, t_3_primes, t_multiply;

    cout << " " << sz << " " << deg << " ";
    random(a, sz, sz, deg);
    random(b, sz, sz, deg);

    if (sz < 20)
    {
        t_waksman = get_time();
        nb = 0;
        do
        {
            multiply_waksman(c2, a, b);
            nb++;
        }
        while ((get_time()-t_waksman) <= thresh);
        t_waksman = (get_time()-t_waksman) / nb;
    }
    else
        t_waksman = 10000.0;
    cout << t_waksman << " ";

    if (deg < 8)
    {
        t_transform = get_time();
        nb = 0;
        do
        {
            multiply_transform(c2, a, b);
            nb++;
        }
        while ((get_time()-t_transform) <= thresh);
        t_transform = (get_time()-t_transform) / nb;
    }
    else
        t_transform = 10000.0;
    cout << t_transform << " ";

    t_geometric = get_time();
    nb = 0;
    do
    {
        multiply_evaluate_geometric(c2, a, b);
        nb++;
    }
    while ((get_time()-t_geometric) <= thresh);
    t_geometric = (get_time()-t_geometric) / nb;
    cout << t_geometric << " ";

    if (deg <= 300)
    {
        t_dense = get_time();
        nb = 0;
        do
        {
            multiply_evaluate_dense(c1, a, b);
            nb++;
        }
        while ((get_time()-t_dense) <= thresh);
        t_dense = (get_time()-t_dense) / nb;
    }
    else
        t_dense = 10000.0;
    cout << t_dense << "(" << (c2 == c1) << ") ";

    t_3_primes = get_time();
    nb = 0;
    do
    {
        multiply_3_primes(c2, a, b);
        nb++;
    }
    while ((get_time()-t_3_primes) <= thresh);
    t_3_primes = (get_time()-t_3_primes) / nb;
    cout << t_3_primes << " ";

    if (is_FFT_ready(NextPowerOfTwo(2*deg - 1)))
    {
        t_FFT = get_time();
        nb = 0;
        do
        {
            multiply_evaluate_FFT(c2, a, b);
            nb++;
        }
        while ((get_time()-t_FFT) <= thresh);
        t_FFT = (get_time()-t_FFT) / nb;
    }
    else
        t_FFT = 10000.0;
    cout << t_FFT << " ";

    t_multiply = get_time();
    nb = 0;
    do
    {
        multiply(c2, a, b);
        nb++;
    }
    while ((get_time()-t_multiply) <= thresh);
    t_multiply = (get_time()-t_multiply) / nb;
    cout << "   " << t_multiply << " ";

    cout << endl;
}


/*------------------------------------------------------------*/
/* for a give prime, checks some (size, degree)               */
/*------------------------------------------------------------*/
void all_checks()
{
    std::vector<long> szs =
    {
        20, 30, 50, 100, 150, 200, 300
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
    // over an FFT prime, FFT always wins
    zz_p::FFTInit(0);
    cout << 0 << endl;
    all_checks();

    zz_p::init(288230376151711813);
    cout << zz_p::modulus() << endl;
    all_checks();

    zz_p::init(23068673);
    cout << zz_p::modulus() << endl;
    all_checks();
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
