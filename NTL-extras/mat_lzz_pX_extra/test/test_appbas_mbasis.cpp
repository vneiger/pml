#include <NTL/lzz_pX.h>
#include <NTL/matrix.h>
#include <NTL/vector.h>
#include <iomanip>
#include <vector>
#include <numeric>
#include <algorithm>
#include <random>
#include <NTL/BasicThreadPool.h>

#include "util.h"
#include "mat_lzz_pX_extra.h"

NTL_CLIENT

/********************************************
 *  tests the approximant basis algorithms  *
 ********************************************/

std::ostream &operator<<(std::ostream &out, const std::vector<long> &s)
{
    out << "[ ";
    for (auto &i: s)
        out << i << " ";
    return out << "]";
}

int main(int argc, char *argv[])
{
    if (argc!=8)
        throw std::invalid_argument("Usage: ./test_appbas_mbasis rdim cdim order nbits verify nthreads thres");
    long threshold = atoi(argv[7]);

    long rdim = atoi(argv[1]);
    long cdim = atoi(argv[2]);
    long order = atoi(argv[3]);
    long nbits = atoi(argv[4]);
    bool verify = (atoi(argv[5])==1);
    SetNumThreads(atoi(argv[6]));

    if (nbits==0)
        zz_p::FFTInit(0);
    else
        zz_p::init(NTL::GenPrime_long(nbits));

    std::cout << "Testing approximant basis (mbasis) with random input matrix" << std::endl;
    std::cout << "--prime =\t" << zz_p::modulus();
    if (nbits==0) std::cout << "  (FFT prime)";
    std::cout << std::endl;
    std::cout << "--rdim =\t" << rdim << std::endl;
    std::cout << "--cdim =\t" << cdim << std::endl;
    std::cout << "--order =\t" << order << std::endl;

    // declare shifts
    Shift shift1(rdim,0); // uniform [0,...,0]
    Shift shift2(rdim); // increasing [0,1,2,..,rdim-1]
    std::iota(shift2.begin(), shift2.end(),0);
    Shift shift3(rdim); // decreasing [rdim,..,3,2,1]
    for (long i = 0; i < rdim; ++i)
        shift3[i] = rdim - i;
    Shift shift4(rdim); // random shuffle of [0,1,...,rdim-1]
    std::iota(shift4.begin(), shift4.end(),0);
    std::shuffle(shift4.begin(), shift4.end(), std::mt19937{std::random_device{}()});
    Shift shift5(rdim); // Hermite shift
    for (long i = 0; i < rdim; ++i)
        shift5[i] = rdim*cdim*order*i;
    Shift shift6(rdim); // reverse Hermite shift
    for (long i = 0; i < rdim; ++i)
        shift6[i] = rdim*cdim*order*(rdim-1-i);
    Shift shift7(rdim);
    for (long i = 0; i < rdim; ++i)
        if (i>=rdim/2)
            shift7[i] = rdim*cdim*order;

    std::vector<Shift> shifts = {shift1, shift2, shift3, shift4, shift5, shift6, shift7};

    std::cout << "Testing approximant basis computation with random input matrix" << std::endl;
    std::cout << "--prime =\t" << zz_p::modulus() << std::endl;
    std::cout << "--rdim =\t" << rdim << std::endl;
    std::cout << "--cdim =\t" << cdim << std::endl;
    std::cout << "--order <\t" << order << std::endl;
    std::cout << "--nthreads =\t" << AvailableThreads() << std::endl;

    // build random matrix
    double t1,t2,t1w,t2w;
    Mat<zz_pX> pmat;
    t1w = GetWallTime(); t1 = GetTime();
    random_mat_zz_pX(pmat, rdim, cdim, order);
    t2w =  GetWallTime(); t2 = GetTime();
    //std::cout << "Time(random mat creation): " << (t2w-t1w) <<  "s,  " << (t2-t1) << "s\n";

    std::cout << "warming up..." << std::endl;
    warmup();

    for (Shift shift : shifts)
    {
        if (rdim<40)
            std::cout << "--shift =\t" << shift << std::endl;
        else
            std::cout << "--shift =\t" << "<length " << rdim << ">" << std::endl;

        // build random matrix
        Mat<zz_pX> pmat;
        t1w = GetWallTime(); t1 = GetTime();
        random_mat_zz_pX(pmat, rdim, cdim, order);
        t2w = GetWallTime(); t2 = GetTime();

        std::cout << "Time(random mat creation): " << (t2w-t1w) << "s,  " << (t2-t1) << "s\n";

        Mat<zz_p> kerbas;
        std::vector<long> pivdeg;

        std::cout << "warming up..." << std::endl;
        warmup();

        Mat<zz_p> mat;
        mat = random_mat_zz_p(rdim, cdim);
        t1w = GetWallTime(); t1 = GetTime();
        Mat<zz_p> kerbas2;
        kernel(kerbas2,mat);
        t2w = GetWallTime(); t2 = GetTime();
        double ref_kernel_wall = t2w-t1w;
        double ref_kernel = t2-t1;
        std::cout << "Time(kernel same size): " << ref_kernel_wall << "s,  " << ref_kernel << "s\n";

        std::cout << "~~~Testing popov_mbasis1 on constant matrix~~~" << std::endl;
        t1w = GetWallTime(); t1 = GetTime();
        pivdeg = popov_mbasis1(kerbas,coeff(pmat,0),shift);
        t2w = GetWallTime(); t2 = GetTime();
        std::cout << "Time(popov_mbasis1 computation): " << (t2w-t1w) << "s,  " << (t2-t1) << "s\n";
        std::cout << "Ratio versus kernel: " << ((t2w-t1w)/ref_kernel_wall) << ", " << ((t2-t1)/ref_kernel) << std::endl;

        if (verify)
        {
            Mat<zz_pX> appbas1;
            appbas1.SetDims(rdim,rdim);
            long row=0;
            for (long i = 0; i < rdim; ++i) {
                if (pivdeg[i]==0) {
                    for (long j = 0; j < rdim; ++j)
                        appbas1[i][j] = kerbas[row][j];
                    ++row;
                } else {
                    SetX(appbas1[i][i]);
                }
            }

            std::cout << "Verifying Popov approximant basis..." << std::endl;
            bool verif1 = is_approximant_basis(appbas1,pmat,1,shift,POPOV,true,false);
            std::cout << (verif1?"correct":"wrong") << std::endl;

            if (std::max(rdim,cdim)<33) {
                Mat<long> degmat;
                degree_matrix(degmat,appbas1,shift,true);
                std::cout << "Print degree matrix of approx basis..." << std::endl;
                std::cout << degmat << std::endl;
            }
        }

        ref_kernel_wall *= order;
        ref_kernel *= order;

        { // plain mbasis
            std::cout << "~~~Testing mbasis_plain~~~" << std::endl;
            t1w = GetWallTime(); t1 = GetTime();
            Mat<zz_pX> appbas;
            pivdeg = mbasis_plain(appbas,pmat,order,shift);
            t2w = GetWallTime(); t2 = GetTime();
            std::cout << "Time(mbasis computation): " << (t2w-t1w) << "s,  " << (t2-t1) << "s\n";
            std::cout << "Ratio versus kernel: " << ((t2w-t1w)/ref_kernel_wall) << ", " << ((t2-t1)/ref_kernel) << std::endl;

            if (verify)
            {
                std::cout << "Verifying ordered weak Popov approximant basis..." << std::endl;
                t1w = GetWallTime(); t1 = GetTime();
                bool verif = is_approximant_basis(appbas,pmat,order,shift,ORD_WEAK_POPOV,true,false);
                t2w = GetWallTime(); t2 = GetTime();
                std::cout << (verif?"correct":"wrong") << std::endl;
                std::cout << "Time(verification): " << (t2w-t1w) << "s,  " << (t2-t1) << "s\n";

                if (std::max(rdim,cdim)<33) {
                    Mat<long> degmat;
                    degree_matrix(degmat,appbas,shift,true);
                    std::cout << "Print degree matrix of approx basis..." << std::endl;
                    std::cout << degmat << std::endl;
                }
            }
        }

        // mbasis Vec<Mat<zz_p>> version
        {
            std::cout << "~~~Testing mbasis ~~~" << std::endl;
            t1w = GetWallTime(); t1 = GetTime();
            Mat<zz_pX> appbas;
            pivdeg = mbasis(appbas,pmat,order,shift);
            t2w = GetWallTime(); t2 = GetTime();

            std::cout << "Time(mbasis computation): " << (t2w-t1w) << "s,  " << (t2-t1) << "s\n";
            std::cout << "Ratio versus kernel: " << ((t2w-t1w)/ref_kernel_wall) << ", " << ((t2-t1)/ref_kernel) << std::endl;

            if (verify)
            {
                std::cout << "Verifying ordered weak Popov approximant basis..." << std::endl;
                t1w = GetWallTime(); t1 = GetTime();
                bool verif = is_approximant_basis(appbas,pmat,order,shift,ORD_WEAK_POPOV,true,false);
                t2w = GetWallTime(); t2 = GetTime();
                std::cout << (verif?"correct":"wrong") << std::endl;
                std::cout << "Time(verification): " << (t2w-t1w) << "s,  " << (t2-t1) << "s\n";

                if (std::max(rdim,cdim)<33) {
                    Mat<long> degmat;
                    degree_matrix(degmat,appbas,shift,true);
                    std::cout << "Print degree matrix of approx basis..." << std::endl;
                    std::cout << degmat << std::endl;
                }
            }
        }

        { // mbasis_mix
            std::cout << "~~~Testing mbasis_mix~~~" << std::endl;
            t1w = GetWallTime(); t1 = GetTime();
            Mat<zz_pX> appbas;
            pivdeg = mbasis_mix(appbas,pmat,order,shift,threshold);
            t2w = GetWallTime(); t2 = GetTime();

            std::cout << "Time(mbasis_mix computation): " << (t2w-t1w) << "s,  " << (t2-t1) << "s\n";
            std::cout << "Ratio versus kernel: " << ((t2w-t1w)/ref_kernel_wall) << ", " << ((t2-t1)/ref_kernel) << std::endl;

            if (verify)
            {
                std::cout << "Verifying ordered weak Popov approximant basis..." << std::endl;
                bool verif = is_approximant_basis(appbas,pmat,order,shift,ORD_WEAK_POPOV,true,false);
                std::cout << (verif?"correct":"wrong") << std::endl;
                std::cout << "Time(verification): " << (t2w-t1w) << "s,  " << (t2-t1) << "s\n";

                if (std::max(rdim,cdim)<33) {
                    Mat<long> degmat;
                    degree_matrix(degmat,appbas,shift,true);
                    std::cout << "Print degree matrix of approx basis..." << std::endl;
                    std::cout << degmat << std::endl;
                }
            }
        }


        // mbasis_generic
        //{
        //    std::cout << "~~~Testing mbasis_generic~~~" << std::endl;
        //    std::cout << "WARNING: here coldim must be 1; order must be multiple of rowdim; shift must be uniform" << std::endl;
        //    t1w = GetWallTime(); t1 = GetTime();
        //    Mat<zz_pX> appbas;
        //    pivdeg = mbasis_generic(appbas,pmat,order,shift);
        //    t2w = GetWallTime(); t2 = GetTime();

        //    std::cout << "Time(mbasis_generic computation): " << (t2w-t1w) << "s,  " << (t2-t1) << "s\n";
        //    std::cout << "Ratio versus kernel: " << ((t2w-t1w)/ref_kernel_wall) << ", " << ((t2-t1)/ref_kernel) << std::endl;

        //    if (verify)
        //    {
        //        std::cout << "Verifying Popov approximant basis..." << std::endl;
        //        t1w = GetWallTime(); t1 = GetTime();
        //        bool verif = is_approximant_basis(appbas,pmat,order,shift,POPOV,true,false);
        //        t2w = GetWallTime(); t2 = GetTime();
        //        std::cout << (verif?"correct":"wrong") << std::endl;
        //        std::cout << "Time(verification): " << (t2w-t1w) << "s,  " << (t2-t1) << "s\n";

        //        if (std::max(rdim,cdim)<33) {
        //            Mat<long> degmat;
        //            degree_matrix(degmat,appbas,shift,true);
        //            std::cout << "Print degree matrix of approx basis..." << std::endl;
        //            std::cout << degmat << std::endl;
        //        }
        //    }
        //}

        // mbasis Popov output
        {
            std::cout << "~~~Testing popov_mbasis~~~" << std::endl;
            t1w = GetWallTime(); t1 = GetTime();
            Mat<zz_pX> appbas;
            pivdeg = popov_mbasis(appbas,pmat,order,shift);
            t2w = GetWallTime(); t2 = GetTime();

            std::cout << "Time(popov_mbasis computation): " << (t2w-t1w) << "s,  " << (t2-t1) << "s\n";
            std::cout << "Ratio versus kernel: " << ((t2w-t1w)/ref_kernel_wall) << ", " << ((t2-t1)/ref_kernel) << std::endl;

            if (verify)
            {
                std::cout << "Verifying Popov approximant basis..." << std::endl;
                t1w = GetWallTime(); t1 = GetTime();
                bool verif = is_approximant_basis(appbas,pmat,order,shift,POPOV,true,false);
                t2w = GetWallTime(); t2 = GetTime();
                std::cout << (verif?"correct":"wrong") << std::endl;
                std::cout << "Time(verification): " << (t2w-t1w) << "s,  " << (t2-t1) << "s\n";

                if (std::max(rdim,cdim)<33) {
                    Mat<long> degmat;
                    degree_matrix(degmat,appbas,shift,true);
                    std::cout << "Print degree matrix of approx basis..." << std::endl;
                    std::cout << degmat << std::endl;
                }
            }
        }
        std::cout << std::endl;
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
