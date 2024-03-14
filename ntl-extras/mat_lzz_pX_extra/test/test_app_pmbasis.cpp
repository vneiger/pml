#include "mat_lzz_pX_approximant.h"
#include "test_examples.h"

#define VERBOSE

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
        throw std::invalid_argument("Usage: ./test_app_pmbasis nbits verbose\n\t--nbits: 0 for FFT prime, between 3 and 60 for normal prime\n\t--verbose: 0/1");

    const long nbits = atoi(argv[1]);
    const bool verbose = (atoi(argv[2])==1);

    if (nbits==0)
        zz_p::FFTInit(0);
    else if (nbits==2)
        zz_p::init(3); // make sure modulus!=2
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

    std::cout << "Testing approximant basis computation (pmbasis)." << std::endl;
    std::cout << "--prime =\t" << zz_p::modulus() << std::endl;

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

                // PMBASIS
                if (verbose)
                    std::cout << "Computation pmbasis... ";

                VecLong rdeg_pmbasis(shift);

                pmbasis(appbas,*pmat,order,rdeg_pmbasis);

                if (verbose)
                    std::cout << "OK. Testing... ";
                if (not is_approximant_basis(appbas,*pmat,order,shift,ORD_WEAK_POPOV,true))
                {
                    std::cout << "Error in pmbasis." << std::endl;
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

                // POPOV_PMBASIS
                if (verbose)
                    std::cout << "Computation popov_pmbasis... ";

                VecLong rdeg_popov_pmbasis(shift);
                popov_pmbasis(appbas,*pmat,order,rdeg_popov_pmbasis);

                if (verbose)
                    std::cout << "OK. Testing... ";
                if (not is_approximant_basis(appbas,*pmat,order,shift,POPOV,true))
                {
                    std::cout << "Error in popov_pmbasis." << std::endl;
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

                // PMBASIS_2x1
                if (rdim==2 && cdim==1)
                {
                    if (verbose)
                        std::cout << "Computation pmbasis_2x1... ";

                    long s0=shift[0];
                    long s1=shift[1];
                    zz_pX p00,p01,p10,p11;
                    zz_pX f0 = (*pmat)[0][0];
                    zz_pX f1 = (*pmat)[1][0];

                    pmbasis_2x1(p00,p01,p10,p11,f0,f1,order,s0,s1);

                    if (verbose)
                        std::cout << "OK. Testing... ";
                    appbas[0][0] = p00;
                    appbas[0][1] = p01;
                    appbas[1][0] = p10;
                    appbas[1][1] = p11;
                    if (not is_approximant_basis(appbas,*pmat,order,shift,ORD_WEAK_POPOV,true))
                    {
                        std::cout << "Error in pmbasis." << std::endl;
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
