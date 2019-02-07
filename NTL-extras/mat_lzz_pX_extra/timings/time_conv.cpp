#include <iomanip>
#include "util.h"
#include "mat_lzz_pX_utils.h"

#define TIME(a)         \
tt = get_time();    \
a                   \
tt = get_time()-tt; \
cout << tt << "\t";

#define SMALL_SUITE


NTL_CLIENT

/******************************************************/
/* Times the degree functions for polynomial matrices */
/******************************************************/

using namespace std;

void run_bench(long m, long n, long d)
{
    double tt, t;
    long nb_iter;

    std::cout << m << "\t" << n << "\t" << d << "\t";

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
        std::cout << t / nb_iter << "\t";
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
        std::cout << t / nb_iter << "\t";
    }

    std::cout << std::endl;
}


int main(int argc, char * argv[])
{
    if (argc!=2 && argc!=5)
        throw std::invalid_argument("Usage: ./time_conv nbits (rdim cdim deg)");

    std::cout << std::fixed;
    std::cout << std::setprecision(8);

    const long nbits = atoi(argv[1]);
    zz_p::init(NTL::GenPrime_long(nbits));

    std::cout << "Bench conv, nbits = " << nbits << ", prime p = " << zz_p::modulus() << std::endl;
    std::cout << "rdim\tcdim\tdeg\tpm->mp\t\tpm->mp(new)\t\tmp->pm\t\tmp->pm(new)" << std::endl;

    if (argc==5)
    {
        const long m = atoi(argv[2]);
        const long n = atoi(argv[3]);
        const long d = atoi(argv[4]);
        run_bench(m, n, d);
    }
    else
    {
        VecLong rdims = {5,10,20,40,70,100,150,200,400,1000,};
        VecLong cdims = {5,10,20,40,70,100,150,200,400,1000,};
        VecLong degs = {5,10,20,40,70,100,150,200,400,1000,2000,4000,8000,16000,};

        warmup();

        for (size_t i = 0; i < rdims.size(); ++i)
            for (size_t j = 0; j < cdims.size(); ++j)
                for (size_t k = 0; k < degs.size(); ++k)
                    run_bench(rdims[i], cdims[j], degs[k]);
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
