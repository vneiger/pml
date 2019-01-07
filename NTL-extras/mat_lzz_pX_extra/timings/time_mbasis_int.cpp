#include <NTL/lzz_pX.h>
#include <NTL/matrix.h>
#include <NTL/vector.h>
#include <iomanip>
#include <NTL/BasicThreadPool.h>

#include "util.h"
#include "mat_lzz_pX_approximant.h"
#include "mat_lzz_pX_interpolant.h"

std::ostream &operator<<(std::ostream &out, const VecLong &s)
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
void one_bench_mbasis(long rdim, long cdim, long npoints)
{
    VecLong shift(rdim,0);

    double t1,t2;

    long nb_iter=0;

    double t_mbasis_rescomp=0.0;
    while (t_mbasis_rescomp<0.1)
    {
        Vec<Mat<zz_p>> evals;
        Vec<zz_p> pts;
        evals.SetLength(npoints);
        for (long pt = 0; pt < npoints; ++pt)
            random(evals[pt], rdim, cdim);
        random(pts, npoints);

        t1 = GetWallTime();
        Mat<zz_pX> intbas;
        VecLong pivdeg = mbasis_rescomp(intbas,evals,pts,shift,0,npoints);
        t2 = GetWallTime();

        t_mbasis_rescomp += t2-t1;
        ++nb_iter;
    }

    t_mbasis_rescomp /= nb_iter;

    double t_mbasis_resupdate=0.0;
    nb_iter=0;
    while (t_mbasis_resupdate<0.1)
    {
        Vec<Mat<zz_p>> evals;
        Vec<zz_p> pts;
        evals.SetLength(npoints);
        for (long pt = 0; pt < npoints; ++pt)
            random(evals[pt], rdim, cdim);
        random(pts, npoints);

        t1 = GetWallTime();
        Mat<zz_pX> intbas;
        VecLong pivdeg = mbasis_resupdate(intbas,evals,pts,shift,0,npoints);
        t2 = GetWallTime();

        t_mbasis_resupdate += t2-t1;
        ++nb_iter;
    }

    t_mbasis_resupdate /= nb_iter;


    double t_mbasis_generic_rescomp=0.0;
    //nb_iter=0;
    //while (t_mbasis_generic_rescomp<0.1)
    //{
    //    Mat<zz_pX> pmat;
    //    random(pmat, rdim, cdim, npoints);
    //    VecLong pivdeg;

    //    t1 = GetWallTime();
    //    Mat<zz_pX> intbas;
    //    mbasis_generic_2n_n_rescomp(intbas,pmat,npoints);
    //    t2 = GetWallTime();

    //    t_mbasis_generic_rescomp += t2-t1;
    //    ++nb_iter;
    //}

    //t_mbasis_generic_rescomp /= nb_iter;

    double t_mbasis_generic_resupdate=0.0;
    //nb_iter=0;
    //while (t_mbasis_generic_resupdate<0.1)
    //{
    //    Mat<zz_pX> pmat;
    //    random(pmat, rdim, cdim, npoints);
    //    VecLong pivdeg;

    //    t1 = GetWallTime();
    //    Mat<zz_pX> intbas;
    //    mbasis_generic_2n_n_resupdate(intbas,pmat,npoints);
    //    t2 = GetWallTime();

    //    t_mbasis_generic_resupdate += t2-t1;
    //    ++nb_iter;
    //}

    //t_mbasis_generic_resupdate /= nb_iter;



    cout << rdim << "\t" << cdim << "\t" << npoints << "\t" << AvailableThreads();
    cout << "\t" << t_mbasis_rescomp << "\t" << t_mbasis_resupdate << "\t" << t_mbasis_generic_rescomp << "\t" << t_mbasis_generic_resupdate;

    double best = t_mbasis_rescomp;
    //if (t_mbasis_rescomp <= t_mbasis_rescomp && t_mbasis_rescomp <= t_mbasis_generic_rescomp)
    //{
    //    cout << "," << "rescomp";
    //    best = t_mbasis_rescomp;
    //}
    //else if (t_mbasis_resupdate <= t_mbasis_rescomp and t_mbasis_resupdate <= t_mbasis_generic_rescomp)
    //{
    //    cout << "," << "resupdate";
    //    best = t_mbasis_resupdate;
    //}
    //else
    //{
    //    cout << "," << "generic";
    //    best = t_mbasis_generic_rescomp;
    //}
    cout << "\t" << t_mbasis_rescomp / best;
    cout << "\t" << t_mbasis_resupdate / best;
    cout << "\t" << t_mbasis_generic_rescomp / best;
    cout << "\t" << t_mbasis_generic_resupdate / best;
    cout << endl;
}

/*------------------------------------------------------------*/
/* run bench on variety of parameters                         */
/*------------------------------------------------------------*/
void run_bench(long nthreads, long nbits, bool fftprime)
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
        else if (nbits < 65)
        {
            zz_p::UserFFTInit(1139410705724735489); // 60 bits
            cout << zz_p::modulus() << ", bit length = " << 60 << endl;
        }
    }
    else
    {
        cout << "Bench mbasis, random prime p = ";
        zz_p::init(NTL::GenPrime_long(nbits));
        cout << zz_p::modulus() << ", bit length = " << nbits << endl;
    }

    VecLong szs = {2, 4, 8, 16, 32, 64, 128, 256, 512, 1024};

    cout << "rdim\tcdim\torder\tnthread\tr-cmp\tr-upd\tg-rcmp\tg-rupd\tratios versus generic resupdate" << endl;
    for (size_t i=0; i<szs.size(); ++i)
    {
        // for generic, currently we can only use dimensions n=2m

        //long interval = ceil( (double)szs[i] / 20);
        //for (long j=1; j<szs[i]; j+=interval)
        {
            long j=szs[i]/2;
            long max_order=128;
            if (szs[i]==512)
                max_order=64;
            else if (szs[i]==1024)
                max_order=32;
            for (long k=2; k<=max_order; k=2*k)
                one_bench_mbasis(szs[i],j,k);
        }
    }
    cout << endl;
}

/*------------------------------------------------------------*/
/* main calls check                                           */
/*------------------------------------------------------------*/
int main(int argc, char ** argv)
{
    std::cout << std::fixed;
    std::cout << std::setprecision(5);

    // TODO one thread for the moment
    VecLong nthreads = {1}; // {1,2,3,4};
    VecLong nbits = {55};
    std::vector<bool> fftprime = {true};

    warmup();

    for (size_t i = 0; i < nthreads.size(); ++i)
        for (size_t j = 0; j < nbits.size(); ++j)
            for (size_t k = 0; k < fftprime.size(); ++k)
                run_bench(nthreads[i],nbits[j],fftprime[k]);

    return 0;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
