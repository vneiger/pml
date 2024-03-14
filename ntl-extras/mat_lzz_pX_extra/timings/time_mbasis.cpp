#include <iomanip>
#include <NTL/BasicThreadPool.h>
#include <numeric> // for std::iota
#include <random> // for mt19937

#include "util.h"
#include "mat_lzz_pX_approximant.h"
#include "mat_lzz_pX_interpolant.h"

static std::ostream &operator<<(std::ostream &out, const VecLong &s)
{
    out << "[ ";
    for (auto &i: s)
        out << i << " ";
    return out << "]";
}

NTL_CLIENT

/*------------------------------------------------------------*/
/* run one bench for specified rdim,cdim,order                */
/*------------------------------------------------------------*/
void one_bench_mbasis(long rdim, long cdim, long degree, long order)
{
    VecLong shift(rdim,0);

    double t1,t2;

    long nb_iter=0;

    double t_mbasis_rescomp=0.0;
    while (t_mbasis_rescomp<0.1)
    {
        Mat<zz_pX> pmat;
        random(pmat, rdim, cdim, degree+1);

        t1 = GetWallTime();
        Mat<zz_pX> appbas;
        VecLong rdeg(shift);
        mbasis_rescomp(appbas,pmat,order,rdeg);
        t2 = GetWallTime();

        t_mbasis_rescomp += t2-t1;
        ++nb_iter;
    }
    t_mbasis_rescomp /= nb_iter;

    double t_mbasis_resupdate=0.0;
    nb_iter=0;
    while (t_mbasis_resupdate<0.1)
    {
        Mat<zz_pX> pmat;
        random(pmat, rdim, cdim, degree+1);
        VecLong pivdeg;

        t1 = GetWallTime();
        Mat<zz_pX> appbas;
        VecLong rdeg(shift);
        mbasis_resupdate(appbas,pmat,order,rdeg);
        t2 = GetWallTime();

        t_mbasis_resupdate += t2-t1;
        ++nb_iter;
    }
    t_mbasis_resupdate /= nb_iter;

    // TODO : when degree<order-1, should probably take a matrix of degree 'degree'
    // and evaluate at 'order' points
    double t_mbasis_int_rescomp=0.0;
    if (order<zz_p::modulus())
    {
        nb_iter=0;
        while (t_mbasis_int_rescomp<0.1)
        {
            Vec<zz_p> pts(INIT_SIZE, order);
            std::iota(pts.begin(), pts.end(), 0);
            std::shuffle(pts.begin(), pts.end(), std::mt19937{std::random_device{}()});
            Vec<Mat<zz_p>> evals(INIT_SIZE, order);
            for (long i = 0; i < order; ++i)
                random(evals[i], rdim, cdim);

            t1 = GetWallTime();
            Mat<zz_pX> intbas;
            mbasis_rescomp(intbas,evals,pts,shift,0,order);
            t2 = GetWallTime();

            t_mbasis_int_rescomp += t2-t1;
            ++nb_iter;
        }
        t_mbasis_int_rescomp /= nb_iter;
    }
    else // not enough distinct points can be found
        t_mbasis_int_rescomp = -1.0;

    // TODO : when degree<order-1, should probably take a matrix of degree 'degree'
    // and evaluate at 'order' points
    double t_mbasis_int_resupdate=0.0;
    if (order<zz_p::modulus())
    {
        nb_iter=0;
        while (t_mbasis_int_resupdate<0.1)
        {
            Vec<zz_p> pts(INIT_SIZE, order);
            std::iota(pts.begin(), pts.end(), 0);
            std::shuffle(pts.begin(), pts.end(), std::mt19937{std::random_device{}()});
            Vec<Mat<zz_p>> evals(INIT_SIZE, order);
            for (long i = 0; i < order; ++i)
                random(evals[i], rdim, cdim);

            t1 = GetWallTime();
            Mat<zz_pX> intbas;
            mbasis_resupdate(intbas,evals,pts,shift,0,order);
            t2 = GetWallTime();

            t_mbasis_int_resupdate += t2-t1;
            ++nb_iter;
        }
        t_mbasis_int_resupdate /= nb_iter;
    }
    else // not enough distinct points can be found
        t_mbasis_int_resupdate = -1.0;

    double t_app_2x1=-1.0;
    if (rdim==2 && cdim==1)
    {
        nb_iter=0;
        t_app_2x1=0.0;
        while (t_app_2x1<0.1)
        {
            zz_pX f0,f1;
            random(f0, degree+1);
            random(f1, degree+1);

            t1 = GetWallTime();
            zz_pX p00,p01,p10,p11;
            long s0 = shift[0];
            long s1 = shift[1];
            appbas_iterative_2x1(p00,p01,p10,p11,f0,f1,order,s0,s1);
            t2 = GetWallTime();

            t_app_2x1 += t2-t1;
            ++nb_iter;
        }
        t_app_2x1 /= nb_iter;
    }


    cout << rdim << "\t" << cdim << "\t" << degree << "\t" << order;
    cout << "\t" << t_mbasis_rescomp << "\t" << t_mbasis_resupdate;
    cout << "\t" << t_mbasis_int_rescomp << "\t" << t_mbasis_int_resupdate;
    cout << "\t" << t_app_2x1;
    cout << endl;
}

/*------------------------------------------------------------*/
/* run bench on variety of parameters                         */
/*------------------------------------------------------------*/
void run_bench(long nthreads, long nbits, bool fftprime, long rdim=-1, long cdim=-1, long degree=-1, long order=-1)
{
    SetNumThreads(nthreads);

    if (fftprime)
    {
        cout << "Bench mbasis, FFT prime p = ";
        if (nbits < 25)
        {
            zz_p::UserFFTInit(786433); // 20 bits
            cout << zz_p::modulus() << ", bit length = " << 20 << endl;
        }
        else if (nbits < 35)
        {
            zz_p::UserFFTInit(2013265921); // 31 bits
            cout << zz_p::modulus() << ", bit length = " << 31 << endl;
        }
        else if (nbits < 45)
        {
            zz_p::UserFFTInit(2748779069441); // 42 bits
            cout << zz_p::modulus() << ", bit length = " << 42 << endl;
        }
        else if (nbits < 61)
        {
            zz_p::FFTInit(0); // 60 bits
            cout << zz_p::modulus() << ", bit length = " << NumBits(zz_p::modulus()) << endl;
        }
        else
        {
            std::cout << "Asking for FFT prime with too large bitsize (> 60). Exiting." << std::endl;
            return;
        }
    }
    else
    {
        cout << "Bench mbasis, random prime p = ";
        zz_p::init(NTL::GenPrime_long(nbits));
        cout << zz_p::modulus() << ", bit length = " << nbits << endl;
    }

    std::cout << "Note: negative timings for interpolant variants indicate that not enough interpolation points could be found in the base field." << std::endl;
    cout << "rdim\tcdim\tdegree\torder\tapp-rcmp\tapp-rupd\tint-rcmp\tint-rupd" << endl;
    
    if (rdim==-1) // then cdim==-1 && order==-1, default case
    {
        VecLong szs = {2, 4, 8, 16, 32, 64, 128, 256, 512, 1024};

        for (size_t i=0; i<szs.size(); ++i)
        {
            VecLong cdims = {szs[i]/4, szs[i]/2, 3*szs[i]/4};
            for (long j : cdims)
                if (j > 0)
                {
                    long max_order=128;
                    if (szs[i]==512)
                        max_order=64;
                    else if (szs[i]==1024)
                        max_order=32;
                    for (long k=2; k<=max_order; k=2*k)
                    {
                        one_bench_mbasis(szs[i],j,k-1,k); // degree ~ order
                        //one_bench_mbasis(szs[i],j,k-1,2*k); // degree ~ order/2
                    }
                }
        }
        cout << endl;
    }
    else
        one_bench_mbasis(rdim,cdim,degree,order);
}

/*------------------------------------------------------------*/
/* main calls check                                           */
/*------------------------------------------------------------*/
int main(int argc, char ** argv)
{
    std::cout << std::fixed;
    std::cout << std::setprecision(8);

    if (argc!=3 and argc!=7)
        throw std::invalid_argument("Usage: ./time_mbasis nbits fftprime (rdim cdim degree order)");
    // assume rdim>0 , cdim>0, degree>0, order>0

    const long nbits = atoi(argv[1]);
    const bool fftprime = (atoi(argv[2])==1);

    if (argc==7)
    {
        const long rdim = atoi(argv[3]);
        const long cdim = atoi(argv[4]);
        const long degree = atoi(argv[5]);
        const long order = atoi(argv[6]);
        run_bench(1,nbits,fftprime,rdim,cdim,degree,order);
    }
    else
        run_bench(1,nbits,fftprime);

    return 0;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
