#include <iomanip>
#include <NTL/BasicThreadPool.h>
#include <numeric>
#include <random>

#include "util.h"
#include "mat_lzz_pX_extra.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* run one bench for specified rdim,cdim,degree               */
/*------------------------------------------------------------*/
void one_bench(long rdim, long cdim, long degA, long degB)
{
    cout << rdim << "\t" << cdim << "\t" << degA << "\t" << degB;

    double t1,t2;

    double t_div=0.0;
    long nb_iter=0;
    while (t_div<0.5)
    {
        // random A
        Mat<zz_pX> A;
        random(A,rdim,cdim,degA+1);
        // random B of given row degree
        Mat<zz_pX> B;
        do
            random(B,rdim,rdim,degB);
        while (not is_row_reduced(B));

        t1 = GetWallTime();
        Mat<zz_pX> Q,R;
        quo_rem(Q,R,A,B);
        t2 = GetWallTime();

        t_div += t2-t1;
        ++nb_iter;
    }
    t_div /= nb_iter;
    cout << "\t" << t_div;

    if (rdim>1)
    {
        // build some non-uniform rdeg with entries degB/2, degB, 3*degB
        VecLong rdegB(rdim);
        for (long i = 0; i < rdim/3; ++i)
            rdegB[i] = degB/2+1;
        for (long i = rdim/3; i < 2*rdim/3; ++i)
            rdegB[i] = degB+1;
        for (long i = 2*rdim/3; i < rdim; ++i)
            rdegB[i] = 3*degB+1;
        std::shuffle(rdegB.begin(), rdegB.end(), std::mt19937{std::random_device{}()});

        t_div=0.0; nb_iter=0;
        while (t_div<0.2)
        {
            // random A
            Mat<zz_pX> A;
            random(A,rdim,cdim,degA+1);
            // random B of given row degree
            Mat<zz_pX> B;
            do
                random_mat_zz_pX_rdeg(B,rdim,rdim,rdegB);
            while (not is_row_reduced(B));

            t1 = GetWallTime();
            Mat<zz_pX> Q,R;
            quo_rem(Q,R,A,B);
            t2 = GetWallTime();

            t_div += t2-t1;
            ++nb_iter;
        }
        t_div /= nb_iter;
        cout << "\t" << t_div;
    }

    cout << endl;
}

/*------------------------------------------------------------*/
/* run bench on variety of parameters                         */
/*------------------------------------------------------------*/
void run_bench(long nthreads, long nbits, bool fftprime, long rdim=-1, long cdim=-1, long degA=-1, long degB=-1)
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

    cout << "rdim\tcdim\tdegA\tdegB\tuniform\t\tnon-uniform" << endl;

    if (rdim==-1) // then cdim==-1 && order==-1, default case
    {
        VecLong szs =
        {
            1, 2, 3, 5, 10, 20, 30, 50, 100
        };

        VecLong degs =
        {
            20, 50, 75, 99, 150, 200, 500, 1000
        };

        for (long sz : szs)
        {
            for (long deg : degs)
            {
                one_bench(sz, 1, 2*deg, deg);
                if (sz>2)
                {
                    one_bench(sz, 1, sz*deg, deg);
                    one_bench(sz, sz/2, 2*deg, deg);
                    one_bench(sz, sz/2, 8*deg, deg);
                    one_bench(sz, sz, 2*deg, deg);
                    one_bench(sz, sz, 8*deg, deg);
                }
            }
        }
        cout << endl;
    }
    else
        one_bench(rdim,cdim,degA,degB);
}

/*------------------------------------------------------------*/
/* main calls check                                           */
/*------------------------------------------------------------*/
int main(int argc, char ** argv)
{
    std::cout << std::fixed;
    std::cout << std::setprecision(8);

    if (argc!=3 and argc!=7)
        throw std::invalid_argument("Usage: ./time_division nbits fftprime (rdim cdim degA degB)");
    // assume rdim>0 , cdim>0, order>0

    const long nbits = atoi(argv[1]);
    const bool fftprime = (atoi(argv[2])==1);

    if (argc==7)
    {
        const long rdim = atoi(argv[3]);
        const long cdim = atoi(argv[4]);
        const long degA = atoi(argv[5]);
        const long degB = atoi(argv[6]);
        warmup();
        run_bench(1,nbits,fftprime,rdim,cdim,degA,degB);
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
