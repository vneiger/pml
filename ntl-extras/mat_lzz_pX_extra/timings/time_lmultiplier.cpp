#include <NTL/lzz_pX.h>
#include <NTL/matrix.h>
#include <NTL/vector.h>
#include <iomanip>

#include "util.h"
#include "mat_lzz_pX_extra.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* times a product (s,s) x (s,1) in degree < deg              */
/*------------------------------------------------------------*/
void one_check(long sz, long deg)
{
    Mat<zz_pX> a, b, c1, c2;
    double s, t_plain;
    long nb;
    const double thres = 0.1;

    cout << sz << "\t" << deg << "\t";

    random(a, sz, sz, deg);
    random(b, sz, 1, deg);

    t_plain = get_time();
    nb = 0;
    do
    {
        multiply(c1, a, b);
        nb++;
    }
    while ((get_time()-t_plain) <= thres);
    t_plain = (get_time()-t_plain) / nb;
    cout << t_plain << "\t";

    if (is_FFT_prime())
    {
        double t_FD, t_FM;

        mat_lzz_pX_lmultiplier_FFT_direct mulFD(a, deg-1);
        t_FD = get_time();
        nb = 0;
        do
        {
            mulFD.multiply(c2, b);
            nb++;
        }
        while ((get_time()-t_FD) <= thres);
        t_FD = (get_time()-t_FD) / nb;
        cout << t_FD << "\t";

        mat_lzz_pX_lmultiplier_FFT_matmul mulFM(a, deg-1);
        t_FM = get_time();
        nb = 0;
        do
        {
            mulFM.multiply(c2, b);
            nb++;
        }
        while ((get_time()-t_FM) <= thres);
        t_FM = (get_time()-t_FM) / nb;
        cout << t_FM << "\t";
        
        s = min(t_FM, t_FD);
    }
    else
    {
        double t_g, t_3, t_d;

        mat_lzz_pX_lmultiplier_3_primes mul3(a, deg-1);
        t_3 = get_time();
        nb = 0;
        do
        {
            mul3.multiply(c2, b);
            nb++;
        }
        while ((get_time()-t_3) <= thres);
        t_3 = (get_time()-t_3) / nb;
        cout << t_3 << "\t";
        
        mat_lzz_pX_lmultiplier_geometric mulg(a, deg-1);
        t_g = get_time();
        nb = 0;
        do
        {
            mulg.multiply(c2, b);
            nb++;
        }
        while ((get_time()-t_g) <= thres);
        t_g = (get_time()-t_g) / nb;
        cout << t_g << "\t";

        mat_lzz_pX_lmultiplier_dense muld(a, deg-1);
        t_d = get_time();
        nb = 0;
        do
        {
            muld.multiply(c2, b);
            nb++;
        }
        while ((get_time()-t_d) <= thres);
        t_d = (get_time()-t_d) / nb;
        cout << t_d << "\t";
        
        s = min(t_d, min(t_g, t_3));
    }

    double t_get;
    std::unique_ptr<mat_lzz_pX_lmultiplier> mul = get_lmultiplier(a, deg-1);
    t_get = get_time();
    nb = 0;
    do
    {
        mul->multiply(c2, b);
        nb++;
    }
    while ((get_time()-t_get) <= thres);
    t_get = (get_time()-t_get) / nb;
    cout << t_get << "(" << t_get/s << ")";

    cout << endl;
}

/*------------------------------------------------------------*/
/* for a give prime, checks some (size, degree)               */
/*------------------------------------------------------------*/
void all_checks()
{

    VecLong szs =
    {
        5, 10, 20, 30, 40, 50, 60, 75, 100, 150, 200, 300
    };

    VecLong degs =
    {
        10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 150, 200, 250, 300, 400, 500, 1000,
    };

    cout << "p=" << zz_p::modulus() << " FFT=" << is_FFT_prime() << endl;
    if (is_FFT_prime())
        std::cout << "sz\tdeg\tplain\tf_dir\tf_mm\tmul" << std::endl;
    else
        std::cout << "sz\tdeg\tplain\t3pri\tgeom\tdens\tmul" << std::endl;

    for (size_t si = 0; si < szs.size(); si++)
        for (size_t di = 0; di < degs.size(); di++)
            one_check(szs[si], degs[di]);

}

/*------------------------------------------------------------*/
/* checks some primes                                         */
/*------------------------------------------------------------*/
void check()
{
    //zz_p::FFTInit(0);
    //all_checks();
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
    std::cout << std::fixed;
    std::cout << std::setprecision(5);
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
