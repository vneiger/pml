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
#include "mat_lzz_pX_approximant.h"

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
    if (argc!=7)
        throw std::invalid_argument("Usage: ./test_appbas_mbasis rdim cdim order nbits verify nthreads");

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
        shift5[i] = cdim*order*i;
    Shift shift6(rdim); // reverse Hermite shift
    for (long i = 0; i < rdim; ++i)
        shift6[i] = (cdim*order+1)*(rdim-1-i);
    Shift shift7(rdim);
    for (long i = 0; i < rdim; ++i)
        if (i>=rdim/2)
            shift7[i] = cdim*order;

    std::vector<Shift> shifts = {shift1, shift2, shift3, shift4, shift5, shift6, shift7};

    std::cout << "Testing approximant basis computation (mbasis) with random input matrix" << std::endl;
    std::cout << "--prime =\t" << zz_p::modulus() << std::endl;
    std::cout << "--rdim =\t" << rdim << std::endl;
    std::cout << "--cdim =\t" << cdim << std::endl;
    std::cout << "--order <\t" << order << std::endl;
    std::cout << "--nthreads =\t" << AvailableThreads() << std::endl;

    double t1,t2,t1w,t2w;

    warmup();

    for (Shift shift : shifts)
    {
        if (rdim<40)
            std::cout << "--shift =\t" << shift << std::endl;
        else
            std::cout << "--shift =\t" << "<length " << rdim << ">" << std::endl;

        std::vector<long> pivdeg;
        double ref_kernel_wall=0.0,ref_kernel=0.0;

        long nb_iter=0;
        while (ref_kernel_wall<0.1)
        {
            Mat<zz_p> mat;
            mat = random_mat_zz_p(rdim, cdim);
            t1w = GetWallTime(); t1 = GetTime();
            Mat<zz_p> kerbas2;
            kernel(kerbas2,mat);
            t2w = GetWallTime(); t2 = GetTime();
            ref_kernel_wall += t2w-t1w;
            ref_kernel += t2-t1;
            ++nb_iter;
        }
        ref_kernel /= nb_iter;
        ref_kernel_wall /= nb_iter;
        std::cout << "Time(kernel same size): " << ref_kernel_wall << "s,  " << ref_kernel << "s\n";

        std::cout << "~~~Testing popov_mbasis1 on constant matrix~~~" << std::endl;

        Mat<zz_pX> pmat_copy;
        Mat<zz_p> kerbas_copy;
        double t_mbasis1w=0.0,t_mbasis1=0.0;
        nb_iter=0;
        while (t_mbasis1w<0.1)
        {
            // build random matrix
            Mat<zz_pX> pmat;
            random(pmat, rdim, cdim, 1);
            t1w = GetWallTime(); t1 = GetTime();
            Mat<zz_p> kerbas;
            pivdeg = popov_mbasis1(kerbas,coeff(pmat,0),shift);
            t2w = GetWallTime(); t2 = GetTime();
            t_mbasis1w += t2w-t1w;
            t_mbasis1 += t2-t1;
            ++nb_iter;
            if (nb_iter==1)
            {
                pmat_copy = pmat;
                kerbas_copy = kerbas;
            }
        }
        t_mbasis1 /= nb_iter;
        t_mbasis1w /= nb_iter;
        std::cout << "Time(popov_mbasis1 computation): " << t_mbasis1w << "s,  " << t_mbasis1 << "s\n";
        std::cout << "Ratio versus kernel: " << (t_mbasis1w/ref_kernel_wall) << ", " << (t_mbasis1/ref_kernel) << std::endl;

        if (verify)
        {
            std::cout << "Verifying Popov approximant basis..." << std::endl;
            Mat<zz_pX> appbas1;
            appbas1.SetDims(rdim,rdim);
            long row=0;
            for (long i = 0; i < rdim; ++i)
            {
                if (pivdeg[i]==0)
                {
                    for (long j = 0; j < rdim; ++j)
                        appbas1[i][j] = kerbas_copy[row][j];
                    ++row;
                }
                else
                    SetX(appbas1[i][i]);
            }

            bool verif1 = is_approximant_basis(appbas1,pmat_copy,1,shift,POPOV,true,false);
            std::cout << (verif1?"correct":"wrong") << std::endl;

            if (std::max(rdim,cdim)<33) {
                Mat<long> degmat;
                degree_matrix(degmat,appbas1,shift,true);
                std::cout << "Print degree matrix of approx basis..." << std::endl;
                std::cout << degmat << std::endl;
            }
        }

        { // plain mbasis
            if (verify) std::cout << "~~~Testing mbasis_plain~~~" << std::endl;
            Mat<zz_pX> appbas_copy, pmat_copy;
            nb_iter=0;
            double t_mbasisw=0.0, t_mbasis=0.0;
            while (t_mbasisw<0.1)
            {
                Mat<zz_pX> pmat;
                random(pmat, rdim, cdim, order);
                t1w = GetWallTime(); t1 = GetTime();
                Mat<zz_pX> appbas;
                pivdeg = mbasis_plain(appbas,pmat,order,shift);
                t2w = GetWallTime(); t2 = GetTime();
                t_mbasisw += t2w-t1w;
                t_mbasis += t2-t1;
                ++nb_iter;
                if (nb_iter==1)
                {
                    pmat_copy = pmat;
                    appbas_copy = appbas;
                }
            }
            t_mbasis /= nb_iter;
            t_mbasisw /= nb_iter;
            std::cout << "Time(mbasis_plain): " << t_mbasisw << "s,  " << t_mbasis << "s\n";
            //std::cout << "Ratio versus kernel: " << (t_mbasisw/ref_kernel_wall) << ", " << (t_mbasis/ref_kernel) << std::endl;

            if (verify)
            {
                std::cout << "Verifying ordered weak Popov approximant basis..." << std::endl;
                t1w = GetWallTime(); t1 = GetTime();
                bool verif = is_approximant_basis(appbas_copy,pmat_copy,order,shift,ORD_WEAK_POPOV,true,false);
                t2w = GetWallTime(); t2 = GetTime();
                std::cout << (verif?"correct":"wrong") << std::endl;
                std::cout << "Time(verification): " << (t2w-t1w) << "s,  " << (t2-t1) << "s\n";

                if (std::max(rdim,cdim)<33) {
                    Mat<long> degmat;
                    degree_matrix(degmat,appbas_copy,shift,true);
                    std::cout << "Print degree matrix of approx basis..." << std::endl;
                    std::cout << degmat << std::endl;
                }
            }
        }

        { // mbasis Vec<Mat<zz_p>> version
            if (verify) std::cout << "~~~Testing mbasis_rescomp~~~" << std::endl;
            Mat<zz_pX> appbas_copy, pmat_copy;
            nb_iter=0;
            double t_mbasisw=0.0, t_mbasis=0.0;
            while (t_mbasisw<0.1)
            {
                Mat<zz_pX> pmat;
                random(pmat, rdim, cdim, order);
                t1w = GetWallTime(); t1 = GetTime();
                Mat<zz_pX> appbas;
                pivdeg = mbasis_rescomp(appbas,pmat,order,shift);
                t2w = GetWallTime(); t2 = GetTime();
                t_mbasisw += t2w-t1w;
                t_mbasis += t2-t1;
                ++nb_iter;
                if (nb_iter==1)
                {
                    pmat_copy = pmat;
                    appbas_copy = appbas;
                }
            }
            t_mbasis /= nb_iter;
            t_mbasisw /= nb_iter;
            std::cout << "Time(mbasis_rescomp): " << t_mbasisw << "s,  " << t_mbasis << "s\n";
            //std::cout << "Ratio versus kernel: " << (t_mbasisw/ref_kernel_wall) << ", " << (t_mbasis/ref_kernel) << std::endl;

            if (verify)
            {
                std::cout << "Verifying ordered weak Popov approximant basis..." << std::endl;
                t1w = GetWallTime(); t1 = GetTime();
                bool verif = is_approximant_basis(appbas_copy,pmat_copy,order,shift,ORD_WEAK_POPOV,true,false);
                t2w = GetWallTime(); t2 = GetTime();
                std::cout << (verif?"correct":"wrong") << std::endl;
                std::cout << "Time(verification): " << (t2w-t1w) << "s,  " << (t2-t1) << "s\n";

                if (std::max(rdim,cdim)<33) {
                    Mat<long> degmat;
                    degree_matrix(degmat,appbas_copy,shift,true);
                    std::cout << "Print degree matrix of approx basis..." << std::endl;
                    std::cout << degmat << std::endl;
                }
            }
        }

        { // mbasis Vec<Mat<zz_p>> version, v2
            if (verify) std::cout << "~~~Testing mbasis_rescomp_v2~~~" << std::endl;
            Mat<zz_pX> appbas_copy, pmat_copy;
            nb_iter=0;
            double t_mbasisw=0.0, t_mbasis=0.0;
            while (t_mbasisw<0.1)
            {
                Mat<zz_pX> pmat;
                random(pmat, rdim, cdim, order);
                t1w = GetWallTime(); t1 = GetTime();
                Mat<zz_pX> appbas;
                pivdeg = mbasis_rescomp_v2(appbas,pmat,order,shift);
                t2w = GetWallTime(); t2 = GetTime();
                t_mbasisw += t2w-t1w;
                t_mbasis += t2-t1;
                ++nb_iter;
                if (nb_iter==1)
                {
                    pmat_copy = pmat;
                    appbas_copy = appbas;
                }
            }
            t_mbasis /= nb_iter;
            t_mbasisw /= nb_iter;
            std::cout << "Time(mbasis_rescomp_v2): " << t_mbasisw << "s,  " << t_mbasis << "s\n";
            //std::cout << "Ratio versus kernel: " << (t_mbasisw/ref_kernel_wall) << ", " << (t_mbasis/ref_kernel) << std::endl;

            if (verify)
            {
                std::cout << "Verifying ordered weak Popov approximant basis..." << std::endl;
                t1w = GetWallTime(); t1 = GetTime();
                bool verif = is_approximant_basis(appbas_copy,pmat_copy,order,shift,ORD_WEAK_POPOV,true,false);
                t2w = GetWallTime(); t2 = GetTime();
                std::cout << (verif?"correct":"wrong") << std::endl;
                std::cout << "Time(verification): " << (t2w-t1w) << "s,  " << (t2-t1) << "s\n";

                if (std::max(rdim,cdim)<33) {
                    Mat<long> degmat;
                    degree_matrix(degmat,appbas_copy,shift,true);
                    std::cout << "Print degree matrix of approx basis..." << std::endl;
                    std::cout << degmat << std::endl;
                }
            }
        }

        { // mbasis_resupdate
            if (verify) std::cout << "~~~Testing mbasis_resupdate~~~" << std::endl;
            Mat<zz_pX> appbas_copy, pmat_copy;
            nb_iter=0;
            double t_mbasisw=0.0, t_mbasis=0.0;
            while (t_mbasisw<0.1)
            {
                Mat<zz_pX> pmat;
                random(pmat, rdim, cdim, order);
                t1w = GetWallTime(); t1 = GetTime();
                Mat<zz_pX> appbas;
                pivdeg = mbasis_resupdate(appbas,pmat,order,shift);
                t2w = GetWallTime(); t2 = GetTime();
                t_mbasisw += t2w-t1w;
                t_mbasis += t2-t1;
                ++nb_iter;
                if (nb_iter==1)
                {
                    pmat_copy = pmat;
                    appbas_copy = appbas;
                }
            }
            t_mbasis /= nb_iter;
            t_mbasisw /= nb_iter;
            std::cout << "Time(mbasis_resupdate): " << t_mbasisw << "s,  " << t_mbasis << "s\n";
            //std::cout << "Ratio versus kernel: " << (t_mbasisw/ref_kernel_wall) << ", " << (t_mbasis/ref_kernel) << std::endl;

            if (verify)
            {
                std::cout << "Verifying ordered weak Popov approximant basis..." << std::endl;
                bool verif = is_approximant_basis(appbas_copy,pmat_copy,order,shift,ORD_WEAK_POPOV,true,false);
                std::cout << (verif?"correct":"wrong") << std::endl;
                std::cout << "Time(verification): " << (t2w-t1w) << "s,  " << (t2-t1) << "s\n";

                if (std::max(rdim,cdim)<33) {
                    Mat<long> degmat;
                    degree_matrix(degmat,appbas_copy,shift,true);
                    std::cout << "Print degree matrix of approx basis..." << std::endl;
                    std::cout << degmat << std::endl;
                }
            }
        }

        { // mbasis
            if (verify) std::cout << "~~~Testing mbasis~~~" << std::endl;
            Mat<zz_pX> appbas_copy, pmat_copy;
            nb_iter=0;
            double t_mbasisw=0.0, t_mbasis=0.0;
            while (t_mbasisw<0.1)
            {
                Mat<zz_pX> pmat;
                random(pmat, rdim, cdim, order);
                t1w = GetWallTime(); t1 = GetTime();
                Mat<zz_pX> appbas;
                pivdeg = mbasis(appbas,pmat,order,shift);
                t2w = GetWallTime(); t2 = GetTime();
                t_mbasisw += t2w-t1w;
                t_mbasis += t2-t1;
                ++nb_iter;
                if (nb_iter==1)
                {
                    pmat_copy = pmat;
                    appbas_copy = appbas;
                }
            }
            t_mbasis /= nb_iter;
            t_mbasisw /= nb_iter;
            std::cout << "Time(mbasis): " << t_mbasisw << "s,  " << t_mbasis << "s\n";
            //std::cout << "Ratio versus kernel: " << (t_mbasisw/ref_kernel_wall) << ", " << (t_mbasis/ref_kernel) << std::endl;

            if (verify)
            {
                std::cout << "Verifying ordered weak Popov approximant basis..." << std::endl;
                bool verif = is_approximant_basis(appbas_copy,pmat_copy,order,shift,ORD_WEAK_POPOV,true,false);
                std::cout << (verif?"correct":"wrong") << std::endl;
                std::cout << "Time(verification): " << (t2w-t1w) << "s,  " << (t2-t1) << "s\n";

                if (std::max(rdim,cdim)<33) {
                    Mat<long> degmat;
                    degree_matrix(degmat,appbas_copy,shift,true);
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

        { // popov_mbasis
            if (verify) std::cout << "~~~Testing popov_mbasis~~~" << std::endl;
            Mat<zz_pX> appbas_copy, pmat_copy;
            nb_iter=0;
            double t_mbasisw=0.0, t_mbasis=0.0;
            while (t_mbasisw<0.1)
            {
                Mat<zz_pX> pmat;
                random(pmat, rdim, cdim, order);
                t1w = GetWallTime(); t1 = GetTime();
                Mat<zz_pX> appbas;
                pivdeg = popov_mbasis(appbas,pmat,order,shift);
                t2w = GetWallTime(); t2 = GetTime();
                t_mbasisw += t2w-t1w;
                t_mbasis += t2-t1;
                ++nb_iter;
                if (nb_iter==1)
                {
                    pmat_copy = pmat;
                    appbas_copy = appbas;
                }
            }
            t_mbasis /= nb_iter;
            t_mbasisw /= nb_iter;
            std::cout << "Time(popov_mbasis): " << t_mbasisw << "s,  " << t_mbasis << "s\n";
            //std::cout << "Ratio versus kernel: " << (t_mbasisw/ref_kernel_wall) << ", " << (t_mbasis/ref_kernel) << std::endl;

            if (verify)
            {
                std::cout << "Verifying ordered weak Popov approximant basis..." << std::endl;
                bool verif = is_approximant_basis(appbas_copy,pmat_copy,order,shift,ORD_WEAK_POPOV,true,false);
                std::cout << (verif?"correct":"wrong") << std::endl;
                std::cout << "Time(verification): " << (t2w-t1w) << "s,  " << (t2-t1) << "s\n";

                if (std::max(rdim,cdim)<33) {
                    Mat<long> degmat;
                    degree_matrix(degmat,appbas_copy,shift,true);
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
