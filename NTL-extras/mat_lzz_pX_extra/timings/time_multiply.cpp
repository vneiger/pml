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
    const double thresh = 1;
    long nb;
    Mat<zz_pX> a, b, c1, c2;
    double t_FFT, t_waksman, t_transform, t_geometric, t_dense, t_dense2, t_3_primes, t_multiply;

    cout << sz << "\t" << deg << "\t";
    random(a, sz, sz, deg);
    random(b, sz, sz, deg);

    t_multiply = get_time();
    nb = 0;
    do
    {
        multiply(c2, a, b);
        nb++;
    }
    while ((get_time()-t_multiply) <= 2*thresh);
    t_multiply = (get_time()-t_multiply) / nb;

    cout << t_multiply << "\t";

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
        if (t_waksman/t_multiply<10.)
            cout << t_waksman/t_multiply << "\t";
        else
            std::cout << ">10" << "\t";
    }
    else
        cout << "inf\t";

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
        if (t_transform/t_multiply<10.)
            cout << t_transform/t_multiply << "\t";
        else
            std::cout << ">10" << "\t";
    }
    else
        cout << "inf\t";

    t_geometric = get_time();
    nb = 0;
    do
    {
        multiply_evaluate_geometric(c2, a, b);
        nb++;
    }
    while ((get_time()-t_geometric) <= thresh);
    t_geometric = (get_time()-t_geometric) / nb;
    if (t_geometric/t_multiply<10.)
        cout << t_geometric/t_multiply << "\t";
    else
        std::cout << ">10" << "\t";

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
        if (t_dense/t_multiply<10.)
            cout << t_dense/t_multiply << "\t";
        else
            std::cout << ">10" << "\t";
    }
    else
        cout << "inf\t";

    if (deg <= 300)
    {
        t_dense2 = get_time();
        nb = 0;
        do
        {
            multiply_evaluate_dense2(c1, a, b);
            nb++;
        }
        while ((get_time()-t_dense2) <= thresh);
        t_dense2 = (get_time()-t_dense2) / nb;
        if (t_dense2/t_multiply<10.)
            cout << t_dense2/t_multiply << "\t";
        else
            std::cout << ">10" << "\t";
    }
    else
        cout << "inf\t";

    t_3_primes = get_time();
    nb = 0;
    do
    {
        multiply_3_primes(c2, a, b);
        nb++;
    }
    while ((get_time()-t_3_primes) <= thresh);
    t_3_primes = (get_time()-t_3_primes) / nb;
    if (t_3_primes/t_multiply<10.)
        cout << t_3_primes/t_multiply << "\t";
    else
        std::cout << ">10" << "\t";

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
        if (t_FFT/t_multiply<10.)
            cout << t_FFT/t_multiply << "\t";
        else
            std::cout << ">10" << "\t";
    }
    else
        cout << "inf\t";

    cout << endl;
}


/*------------------------------------------------------------*/
/* for a give prime, checks some (size, degree)               */
/*------------------------------------------------------------*/
void all_checks()
{
    VecLong szs =
    {
        2, 4, 8, 14, 20, 30, 50, 100, 150, 200
    };

    VecLong degs =
    {
        1, 2, 3, 4, 5, 10, 15, 20, 25, 50, 60, 70, 100, 150, 200, 250, 300, 400, 1000, 2000
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
    zz_p::UserFFTInit(23068673);
    cout << endl << "20 bit FFT prime " << zz_p::modulus() << endl;
    cout << "sz\tdeg\ttmul\twak\ttrsf\tgeom\tdens\tdens2\t3prm\tfft" << endl;
    one_check(2,10000);

    // over an FFT prime, FFT always wins
    zz_p::FFTInit(0);
    cout << endl << "60 bit FFT prime" << endl;
    cout << "sz\tdeg\ttmul\twak\ttrsf\tgeom\tdens\tdens2\t3prm\tfft" << endl;
    //all_checks();
    one_check(2,1000000);

    zz_p::init(288230376151711813);
    cout << endl << "60 bit prime " << zz_p::modulus() << endl;
    cout << "sz\tdeg\ttmul\twak\ttrsf\tgeom\tdens\tdens2\t3prm\tfft" << endl;
    //all_checks();
    one_check(2,1000000);

    zz_p::init(23068673);
    cout << endl << "20 bit prime " << zz_p::modulus() << endl;
    cout << "sz\tdeg\ttmul\twak\ttrsf\tgeom\tdens\tdens2\t3prm\tfft" << endl;
    //all_checks();
    one_check(2,1000000);
}  

/*------------------------------------------------------------*/
/* main calls check                                           */
/*------------------------------------------------------------*/
int main(int argc, char ** argv)
{
    cout << std::fixed;
    cout << std::setprecision(5);
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
