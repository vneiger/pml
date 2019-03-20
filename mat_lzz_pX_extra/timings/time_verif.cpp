#include <iomanip>
#include <NTL/BasicThreadPool.h>
#include <numeric>
#include <random>

#include "util.h"
#include "mat_lzz_pX_approximant.h"
#include "mat_lzz_pX_interpolant.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* run one bench for specified rdim,cdim,order                */
/*------------------------------------------------------------*/
void one_bench_pmbasis(long rdim, long cdim, long degree, long order)
{
    VecLong shift(rdim,0);

    double t1,t2;

    long nb_iter=0;

    double t_pmbasis_app=0.0;
    double t_verif_app=0.0;
    while (t_pmbasis_app<0.2)
    {
        Mat<zz_pX> pmat;
        random(pmat, rdim, cdim, degree+1);
        Mat<zz_pX> copy(pmat);

        t1 = GetWallTime();
        Mat<zz_pX> appbas;
        VecLong rdeg(shift);
        pmbasis(appbas,pmat,order,rdeg);
        t2 = GetWallTime();

        t_pmbasis_app += t2-t1;

        t1 = GetWallTime();
        is_approximant_basis(appbas, copy, order, shift, ORD_WEAK_POPOV, true);
        t2 = GetWallTime();

        t_verif_app += t2-t1;

        ++nb_iter;
    }
    t_pmbasis_app /= nb_iter;
    t_verif_app /= nb_iter;

    nb_iter=0;
    double t_pmbasis_intgeom=0.0;
    double t_verif_int=0.0;
    if (2*order+1 < zz_p::modulus())
    {
        nb_iter=0;
        while (t_pmbasis_intgeom<0.2)
        {
            Mat<zz_pX> pmat;
            random(pmat, rdim, cdim, degree+1);
            // geometric in degree 'order' (bound on degree of intbas) requires an
            // element order at least 2*deg+1
            zz_p r;
            element_of_order(r, 2*order+1); 
            if (IsZero(r))
            {
                t_pmbasis_intgeom=-1.0; nb_iter=1;
                break;
            }

            t1 = GetWallTime();
            Vec<zz_p> pts;
            Mat<zz_pX> intbas;
            VecLong rdeg(shift);
            pmbasis_geometric(intbas,pmat,r,order,rdeg,pts);
            t2 = GetWallTime();
        
            t_pmbasis_intgeom += t2-t1;

            t1 = GetWallTime();
            is_interpolant_basis_geometric(intbas, pmat, r, order, shift, ORD_WEAK_POPOV, false);
            t2 = GetWallTime();

            t_verif_int += t2-t1;

            ++nb_iter;
        }
        t_pmbasis_intgeom /= nb_iter;
        t_verif_int /= nb_iter;
    }
    else
        t_pmbasis_intgeom=-1.0;

    std::cout << rdim << "\t" << cdim << "\t" << order << "\t";
    std::cout << t_pmbasis_app << "\t" << t_verif_app << "\t";
    std::cout << t_pmbasis_intgeom << "\t" << t_verif_int << "\t" << std::endl;
}

/*------------------------------------------------------------*/
/* run bench on variety of parameters                         */
/*------------------------------------------------------------*/
void run_bench(long nthreads, long nbits, bool fftprime, long rdim=-1, long cdim=-1, long degree=-1, long order=-1)
{
    SetNumThreads(nthreads);

    if (fftprime)
    {
        cout << "Bench pmbasis, FFT prime p = ";
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
            zz_p::FFTInit(0);
            std::cout << zz_p::modulus() << ", bit length = " << NumBits(zz_p::modulus()) << std::endl;
        }
        else
        {
            std::cout << "Asking for FFT prime with too large bitsize (> 60). Exiting." << std::endl;
            return;
        }
    }
    else
    {
        cout << "Bench pmbasis, random prime p = ";
        zz_p::init(NTL::GenPrime_long(nbits));
        cout << zz_p::modulus() << ", bit length = " << nbits << endl;
    }

    std::cout << "Note: negative timings for interpolant variants indicate that not enough interpolation points could be found in the base field." << std::endl;
    cout << "rdim\tcdim\torder\tapp\t\tverif-app\tint\t\tverif-int" << endl;

    if (rdim==-1) // then cdim==-1 && order==-1, default case
    {
        VecLong szs = {2, 4, 8, 16, 32, 64, 128, 256, 512};

        for (size_t i=0; i<szs.size(); ++i)
        {
            //long interval = ceil( (double)szs[i] / 4);
            //for (long j=1; 2*j<3*szs[i]; j+=interval)
            long j = szs[i]/2;
            {
                long max_order= 65536/szs[i];
                for (long k=32; k<=max_order; k=2*k)
                {
                    one_bench_pmbasis(szs[i],j,k-1,k); // degree ~ order
                    //one_bench_pmbasis(szs[i],j,k-1,2*k); // degree ~ order/2
                }
            }
        }
        cout << endl;
    }
    else
        one_bench_pmbasis(rdim,cdim,degree,order);
}

/*------------------------------------------------------------*/
/* main calls check                                           */
/*------------------------------------------------------------*/
int main(int argc, char ** argv)
{
    std::cout << std::fixed;
    std::cout << std::setprecision(8);

    if (argc!=3 and argc!=7)
        throw std::invalid_argument("Usage: ./time_verif nbits fftprime (rdim cdim degree order)");
    // assume rdim>0 , cdim>0, order>0

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
