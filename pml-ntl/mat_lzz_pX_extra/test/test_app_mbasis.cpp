#include "mat_lzz_pX_approximant.h"
#include "test_examples.h"

NTL_CLIENT

std::ostream &operator<<(std::ostream &out, const VecLong &s)
{
    out << "[ ";
    for (auto &i: s)
        out << i << " ";
    return out << "]";
}

int main(int argc, char *argv[])
{
    if (argc!=3)
        throw std::invalid_argument("Usage: ./test_app_mbasis nbits verbose\n\t--nbits: 0 for FFT prime, between 3 and 60 for normal prime\n\t--verbose: 0/1");

    const long nbits = atoi(argv[1]);
    const bool verbose = (atoi(argv[2])==1);

    if (nbits==0)
        zz_p::FFTInit(0);
    else if (nbits==2)
        //zz_p::init(3); // make sure modulus!=2
        zz_p::init(2);
    else
        zz_p::init(NTL::GenPrime_long(nbits));

    if (!verbose)
    {
        std::cout << "Launching tests on many orders / shifts / dimensions / types of matrices." << std::endl;
        std::cout << "This may take several minutes." << std::endl;
    }

    // build couple (test_matrices, test_shifts)
    std::pair<std::vector<Mat<zz_pX>>, std::vector<std::vector<VecLong>>>
    test_examples = build_test_examples();

    std::cout << "Testing approximant basis computation (mbasis)." << std::endl;
    std::cout << "--prime =\t" << zz_p::modulus() << std::endl;

    VecLong pivdeg; 
    Mat<zz_pX> appbas;

    size_t i=0;
    size_t inst=0;
    for (auto pmat = test_examples.first.begin(); pmat!= test_examples.first.end(); ++pmat, ++i)
    {
        const long rdim = pmat->NumRows();
        const long cdim = pmat->NumCols();
        const long d = deg(*pmat);
        VecLong orders;
        if (d <= 1)
            orders = {1, 2, 5}; // recall the order must be (strictly) positive
        else
            orders = {d/2,d+1,2*d};

        for (long order : orders)
        {
            for (VecLong shift : test_examples.second[i])
            {
                ++inst;
                if (verbose)
                    std::cout << std::endl << "instance number " << inst << ":" << std::endl;

                if (verbose)
                {
                    std::cout << "--rdim =\t" << rdim << std::endl;
                    std::cout << "--cdim =\t" << cdim << std::endl;
                    std::cout << "--deg =\t" << d << std::endl;
                    std::cout << "--order =\t" << order << std::endl;
                    std::cout << "--shift =\t" << shift << std::endl;
                }

                // POPOV_MBASIS1
                if (verbose)
                    std::cout << "Computation popov_mbasis1... ";

                Mat<zz_p> kerbas;

                pivdeg = popov_mbasis1(kerbas,coeff(*pmat,0),shift);

                if (verbose)
                    std::cout << "OK. Testing... ";

                // build approx basis from kerbas
                clear(appbas);
                appbas.SetDims(rdim,rdim);
                long row=0;
                for (long r = 0; r < rdim; ++r)
                {
                    if (pivdeg[r]==0)
                    {
                        for (long j = 0; j < rdim; ++j)
                            appbas[r][j] = kerbas[row][j];
                        ++row;
                    }
                    else
                        SetX(appbas[r][r]);
                }

                if (not is_approximant_basis(appbas,*pmat,1,shift,POPOV,true))
                {
                    std::cout << "Error in popov_mbasis1." << std::endl;
                    std::cout << "--rdim =\t" << rdim << std::endl;
                    std::cout << "--cdim =\t" << cdim << std::endl;
                    std::cout << "--deg =\t" << d << std::endl;
                    std::cout << "--order =\t" << order << std::endl;
                    std::cout << "--shift =\t" << shift << std::endl;
                    std::cout << zz_p::modulus() << std::endl;
                    std::cout << "Input: " << std::endl << coeff(*pmat,0) << std::endl;
                    std::cout << "Approx basis: " << std::endl << appbas << std::endl;
                    std::cout << "Kernel : " << std::endl << kerbas << std::endl;
                    std::cout << std::endl << pivdeg << std::endl;
                    return 0;
                }
                if (verbose)
                    std::cout << "OK." << std::endl;

                // MBASIS1
                if (verbose)
                    std::cout << "Computation mbasis1... ";

                pivdeg = mbasis1(kerbas,coeff(*pmat,0),shift);

                if (verbose)
                    std::cout << "OK. Testing... ";

                // build approx basis from kerbas
                clear(appbas);
                appbas.SetDims(rdim,rdim);
                row=0;
                for (long r = 0; r < rdim; ++r)
                {
                    if (pivdeg[r]==0)
                    {
                        for (long j = 0; j < rdim; ++j)
                            appbas[r][j] = kerbas[row][j];
                        ++row;
                    }
                    else
                        SetX(appbas[r][r]);
                }

                if (not is_approximant_basis(appbas,*pmat,1,shift,ORD_WEAK_POPOV,true))
                {
                    std::cout << "Error in mbasis1." << std::endl;
                    std::cout << "--rdim =\t" << rdim << std::endl;
                    std::cout << "--cdim =\t" << cdim << std::endl;
                    std::cout << "--deg =\t" << d << std::endl;
                    std::cout << "--order =\t" << order << std::endl;
                    std::cout << "--shift =\t" << shift << std::endl;
                    std::cout << zz_p::modulus() << std::endl;
                    std::cout << "Input: " << std::endl << coeff(*pmat,0) << std::endl;
                    std::cout << "Approx basis: " << std::endl << appbas << std::endl;
                    std::cout << "Kernel : " << std::endl << kerbas << std::endl;
                    std::cout << std::endl << pivdeg << std::endl;
                    return 0;
                }
                if (verbose)
                    std::cout << "OK." << std::endl;

                // 2 x 2 specific
                if (rdim==2 && cdim==1)
                {
                    if (verbose)
                        std::cout << "Computation with 2x2 specific code... ";

                    long s0 = shift[0];
                    long s1 = shift[1];
                    zz_pX f0 = (*pmat)[0][0];
                    zz_pX f1 = (*pmat)[1][0];
                    zz_pX p00,p01,p10,p11;
                    appbas_iterative_2x1(p00,p01,p10,p11,f0,f1,order,s0,s1);

                    if (verbose)
                        std::cout << "OK. Testing... ";
                    appbas[0][0] = p00;
                    appbas[0][1] = p01;
                    appbas[1][0] = p10;
                    appbas[1][1] = p11;
                    if (not is_approximant_basis(appbas,*pmat,order,shift,ORD_WEAK_POPOV,true))
                    {
                        std::cout << "Error in appbas_iterative_dim2." << std::endl;
                        std::cout << "--rdim =\t" << rdim << std::endl;
                        std::cout << "--cdim =\t" << cdim << std::endl;
                        std::cout << "--deg =\t" << d << std::endl;
                        std::cout << "--order =\t" << order << std::endl;
                        std::cout << "--shift =\t" << shift << std::endl;
                        std::cout << zz_p::modulus() << std::endl;
                        std::cout << *pmat << std::endl;
                        std::cout << appbas << std::endl;
                        return 0;
                    }
                    if (verbose)
                        std::cout << "OK." << std::endl;
                }

                // PLAIN MBASIS
                if (verbose)
                    std::cout << "Computation mbasis_plain... ";

                VecLong rdeg_plain(shift);
                mbasis_plain(appbas,*pmat,order,rdeg_plain);

                if (verbose)
                    std::cout << "OK. Testing... ";
                if (not is_approximant_basis(appbas,*pmat,order,shift,ORD_WEAK_POPOV,true))
                {
                    std::cout << "Error in mbasis_plain." << std::endl;
                    std::cout << "--rdim =\t" << rdim << std::endl;
                    std::cout << "--cdim =\t" << cdim << std::endl;
                    std::cout << "--deg =\t" << d << std::endl;
                    std::cout << "--order =\t" << order << std::endl;
                    std::cout << "--shift =\t" << shift << std::endl;
                    std::cout << zz_p::modulus() << std::endl;
                    std::cout << *pmat << std::endl;
                    std::cout << appbas << std::endl;
                    return 0;
                }
                if (verbose)
                    std::cout << "OK." << std::endl;

                // MBASIS_RESCOMP
                if (verbose)
                    std::cout << "Computation mbasis_rescomp... ";

                VecLong rdeg_rescomp(shift);
                mbasis_rescomp(appbas,*pmat,order,rdeg_rescomp);

                if (verbose)
                    std::cout << "OK. Testing... ";
                if (not is_approximant_basis(appbas,*pmat,order,shift,ORD_WEAK_POPOV,true))
                {
                    std::cout << "Error in mbasis_rescomp." << std::endl;
                    std::cout << "--rdim =\t" << rdim << std::endl;
                    std::cout << "--cdim =\t" << cdim << std::endl;
                    std::cout << "--deg =\t" << d << std::endl;
                    std::cout << "--order =\t" << order << std::endl;
                    std::cout << "--shift =\t" << shift << std::endl;
                    std::cout << zz_p::modulus() << std::endl;
                    std::cout << *pmat << std::endl;
                    std::cout << appbas << std::endl;
                    return 0;
                }
                if (verbose)
                    std::cout << "OK." << std::endl;

                // MBASIS_RESUPDATE
                if (verbose)
                    std::cout << "Computation mbasis_resupdate... ";

                VecLong rdeg_resupdate(shift);
                mbasis_resupdate(appbas,*pmat,order,rdeg_resupdate);

                if (verbose)
                    std::cout << "OK. Testing...";
                if (not is_approximant_basis(appbas,*pmat,order,shift,ORD_WEAK_POPOV,true))
                {
                    std::cout << "Error in mbasis_resupdate." << std::endl;
                    std::cout << "--rdim =\t" << rdim << std::endl;
                    std::cout << "--cdim =\t" << cdim << std::endl;
                    std::cout << "--deg =\t" << d << std::endl;
                    std::cout << "--order =\t" << order << std::endl;
                    std::cout << "--shift =\t" << shift << std::endl;
                    std::cout << zz_p::modulus() << std::endl;
                    std::cout << *pmat << std::endl;
                    std::cout << appbas << std::endl;
                    return 0;
                }
                if (verbose)
                    std::cout << "OK." << std::endl;

                // POPOV_MBASIS
                if (verbose)
                    std::cout << "Computation popov_mbasis... ";

                VecLong rdeg_popov(shift);
                popov_mbasis(appbas,*pmat,order,rdeg_popov);

                if (verbose)
                    std::cout << "OK. Testing...";
                if (not is_approximant_basis(appbas,*pmat,order,shift,POPOV,true))
                {
                    std::cout << "Error in popov_mbasis." << std::endl;
                    std::cout << "--rdim =\t" << rdim << std::endl;
                    std::cout << "--cdim =\t" << cdim << std::endl;
                    std::cout << "--deg =\t" << d << std::endl;
                    std::cout << "--order =\t" << order << std::endl;
                    std::cout << "--shift =\t" << shift << std::endl;
                    std::cout << zz_p::modulus() << std::endl;
                    std::cout << *pmat << std::endl;
                    std::cout << appbas << std::endl;
                    return 0;
                }
                if (verbose)
                    std::cout << "OK." << std::endl;
            }
        }
    }
    std::cout << inst << " instances processed with success." << std::endl;
    return 0;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
