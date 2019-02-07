#include <iomanip>
#include "util.h"
#include "mat_lzz_pX_utils.h"

NTL_CLIENT

/******************************************************/
/* Times the degree functions for polynomial matrices */
/******************************************************/

using namespace std;

void run_bench(long m, long n, long d)
{
    double tt, t;
    long nb_iter;

    double tconv1, tconv2;

    std::cout << m << "\t" << n << "\t" << d << "\t";


    { // WARMUP pmat --> matp
        Mat<zz_pX> pmat; random(pmat, m, n, d);
        t=0.0, nb_iter=0;
        while (t<0.2)
        {
            tt = GetWallTime();
            Vec<Mat<zz_p>> matp;
            conv(matp, pmat);
            t += GetWallTime()-tt;
            ++nb_iter;
        }
    }

    { // WARMUP NEW pmat --> matp
        Mat<zz_pX> pmat; random(pmat, m, n, d);
        t=0.0, nb_iter=0;
        while (t<0.2)
        {

            tt = GetWallTime();
            Vec<Mat<zz_p>> matp;
            conv_new(matp, pmat);
            t += GetWallTime()-tt;
            ++nb_iter;
        }
    }


    { // pmat --> matp
        Mat<zz_pX> pmat; random(pmat, m, n, d);
        t=0.0, nb_iter=0;
        while (t<0.2)
        {

            tt = GetWallTime();
            Vec<Mat<zz_p>> matp;
            conv(matp, pmat);
            t += GetWallTime()-tt;
            ++nb_iter;
        }
        //std::cout << t / nb_iter << "\t";
        tconv1 = t/nb_iter;
    }

    { // NEW pmat --> matp
        Mat<zz_pX> pmat; random(pmat, m, n, d);
        t=0.0, nb_iter=0;
        while (t<0.2)
        {

            tt = GetWallTime();
            Vec<Mat<zz_p>> matp;
            conv_new(matp, pmat);
            t += GetWallTime()-tt;
            ++nb_iter;
        }
        t /= nb_iter;
        std::cout << t / tconv1 << "\t";
    }


    { // NEWNEW pmat --> matp
        Mat<zz_pX> pmat; random(pmat, m, n, d);
        t=0.0, nb_iter=0;
        while (t<0.2)
        {

            tt = GetWallTime();
            Vec<Mat<zz_p>> matp;
            conv_new2(matp, pmat);
            t += GetWallTime()-tt;
            ++nb_iter;
        }
        t /= nb_iter;
        std::cout << t / tconv1 << "\t";
    }


    { // matp --> pmat
        Vec<Mat<zz_p>> matp; random(matp, m, n, d);
        t=0.0, nb_iter=0;
        while (t<0.2)
        {
            tt = GetWallTime();
            Mat<zz_pX> pmat;
            conv(pmat, matp);
            t += GetWallTime()-tt;
            ++nb_iter;
        }
        tconv2 = t/nb_iter;
        //std::cout << t / nb_iter << "\t";
    }

    { // NEW matp --> pmat
        Vec<Mat<zz_p>> matp; random(matp, m, n, d);
        t=0.0, nb_iter=0;
        while (t<0.2)
        {
            tt = GetWallTime();
            Mat<zz_pX> pmat;
            conv_new(pmat, matp);
            t += GetWallTime()-tt;
            ++nb_iter;
        }
        t /= nb_iter;
        std::cout << t / tconv2 << "\t";
    }


    { // NEW2 matp --> pmat
        Vec<Mat<zz_p>> matp; random(matp, m, n, d);
        t=0.0, nb_iter=0;
        while (t<0.2)
        {
            tt = GetWallTime();
            Mat<zz_pX> pmat;
            conv_new2(pmat, matp);
            t += GetWallTime()-tt;
            ++nb_iter;
        }
        t /= nb_iter;
        std::cout << t / tconv2 << "\t";
    }



    std::cout << std::endl;
}


int main(int argc, char * argv[])
{
    if (argc!=2 && argc!=5)
        throw std::invalid_argument("Usage: ./time_conv nbits (rdim cdim deg)");

    std::cout << std::fixed;
    std::cout << std::setprecision(5);

    const long nbits = atoi(argv[1]);
    zz_p::init(NTL::GenPrime_long(nbits));

    std::cout << "Bench conv, nbits = " << nbits << ", prime p = " << zz_p::modulus() << std::endl;
    //std::cout << "rdim\tcdim\tdeg\tpm->mp(old)\tpm->mp(new)\tpm->mp(new2)\tmp->pm(old)\tmp->pm(new)\tmp->pm(new2)" << std::endl;
    std::cout << "rdim\tcdim\tdeg\tpm->mp1\tpm->mp2\tmp->pm1\tmp->pm2" << std::endl;

    if (argc==5)
    {
        const long m = atoi(argv[2]);
        const long n = atoi(argv[3]);
        const long d = atoi(argv[4]);
        warmup();
        run_bench(m, n, d);
    }
    else
    {
        VecLong rdims = {5,10,20,50,100,200,400,1000,};
        VecLong cdims = {5,10,20,50,100,200,400,1000,};
        VecLong degs = {5,10,20,50,100,200,400,1000,2000,4000,8000,16000,};

        for (size_t m : rdims)
            for (size_t n : cdims)
                for (size_t d : degs)
                    if (m*n*d < 500000000)
                        run_bench(m, n, d);
    }

    return 0;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
