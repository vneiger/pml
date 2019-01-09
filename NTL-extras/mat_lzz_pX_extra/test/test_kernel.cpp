#include "mat_lzz_pX_kernel.h"
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
        throw std::invalid_argument("Usage: ./test_kernel nbits verbose\n\t--nbits: 0 for FFT prime, between 3 and 60 for normal prime\n\t--verbose: 0/1");

    const long nbits = atoi(argv[1]);
    const bool verbose = (atoi(argv[2])==1);

    if (nbits==0)
        zz_p::FFTInit(0);
    else
        zz_p::init(NTL::GenPrime_long(nbits));

    if (!verbose)
    {
        std::cout << "Launching tests on many shifts / dimensions / types of matrices." << std::endl;
        std::cout << "If nothing gets printed below, this means all tests were passed." << std::endl;
        std::cout << "This may take several minutes." << std::endl;
    }

    // build couple (test_matrices, test_shifts)
    std::pair<std::vector<Mat<zz_pX>>, std::vector<std::vector<VecLong>>>
    test_examples = build_test_examples();

    std::cout << "Testing kernel basis computation." << std::endl;
    std::cout << "--prime =\t" << zz_p::modulus() << std::endl;

    VecLong pivdeg; 
    Mat<zz_pX> kerbas;

    size_t i=0;
    for (auto pmat = test_examples.first.begin(); pmat!= test_examples.first.end(); ++pmat, ++i)
    {
        if (verbose)
            std::cout << i << std::endl;
        
        const long rdim = pmat->NumRows();
        const long cdim = pmat->NumCols();
        const long d = deg(*pmat);

        for (VecLong shift : test_examples.second[i])
        {
            if (verbose)
            {
                std::cout << "--rdim =\t" << rdim << std::endl;
                std::cout << "--cdim =\t" << cdim << std::endl;
                std::cout << "--deg =\t" << d << std::endl;
                std::cout << "--shift =\t" << shift << std::endl;
            }

            { // direct via approximation
                if (verbose)
                    std::cout << "Computation of the kernel via approximation... ";

                Mat<zz_pX> copy_pmat(*pmat);
                Mat<zz_pX> kerbas;
                VecLong pivind;
                VecLong rdeg(shift);
                kernel_basis_via_approximation(kerbas,pivind,copy_pmat,rdeg);

                if (verbose)
                    std::cout << "OK. Testing... ";
                if (not is_kernel_basis(kerbas,*pmat,shift,ORD_WEAK_POPOV,true))
                {
                    std::cout << "Error in kernel_basis_via_approximation." << std::endl;
                    std::cout << "--rdim =\t" << rdim << std::endl;
                    std::cout << "--cdim =\t" << cdim << std::endl;
                    std::cout << "--deg =\t" << d << std::endl;
                    std::cout << "--shift =\t" << shift << std::endl;
                    std::cout << zz_p::modulus() << std::endl;
                    std::cout << "Input: " << std::endl << *pmat << std::endl;
                    std::cout << "Kernel : " << std::endl << kerbas << std::endl;
                    std::cout << "pivot degree: " << std::endl << pivdeg << std::endl;
                    return 0;
                }
                if (verbose)
                    std::cout << "OK." << std::endl;
            }

            { // zls - approximation
                if (verbose)
                    std::cout << "Computation of the kernel via ZLS algorithm (using approximation)... ";

                Mat<zz_pX> copy_pmat(*pmat);
                Mat<zz_pX> kerbas;
                VecLong rdeg(shift);
                VecLong pivind;
                kernel_basis_zls_via_approximation(kerbas,pivind,copy_pmat,rdeg);

                if (verbose)
                    std::cout << "OK. Testing... ";
                if (not is_kernel_basis(kerbas,*pmat,shift,ORD_WEAK_POPOV,true))
                {
                    std::cout << "Error in kernel_basis_zls_via_approximation." << std::endl;
                    std::cout << "--rdim =\t" << rdim << std::endl;
                    std::cout << "--cdim =\t" << cdim << std::endl;
                    std::cout << "--deg =\t" << d << std::endl;
                    std::cout << "--shift =\t" << shift << std::endl;
                    std::cout << zz_p::modulus() << std::endl;
                    std::cout << "Input: " << std::endl << *pmat << std::endl;
                    std::cout << "Kernel : " << std::endl << kerbas << std::endl;
                    std::cout << "pivot degree: " << std::endl << pivdeg << std::endl;
                    return 0;
                }
                if (verbose)
                    std::cout << "OK." << std::endl;
            }

            { // zls - interpolation
                if (verbose)
                    std::cout << "Computation of the kernel via ZLS algorithm (using interpolation)... ";

                Mat<zz_pX> copy_pmat(*pmat);
                VecLong rdeg(shift);
                Mat<zz_pX> kerbas;
                VecLong pivind;
                kernel_basis_zls_via_interpolation(kerbas,pivind,copy_pmat,rdeg);

                if (verbose)
                    std::cout << "OK. Testing... ";
                if (not is_kernel_basis(kerbas,*pmat,shift,ORD_WEAK_POPOV,true))
                {
                    std::cout << "Error in kernel_basis_zls_via_interpolation." << std::endl;
                    std::cout << "--rdim =\t" << rdim << std::endl;
                    std::cout << "--cdim =\t" << cdim << std::endl;
                    std::cout << "--deg =\t" << d << std::endl;
                    std::cout << "--shift =\t" << shift << std::endl;
                    std::cout << zz_p::modulus() << std::endl;
                    std::cout << "Input: " << std::endl << *pmat << std::endl;
                    std::cout << "Kernel : " << std::endl << kerbas << std::endl;
                    std::cout << "pivot degree: " << std::endl << pivdeg << std::endl;
                    return 0;
                }
                if (verbose)
                    std::cout << "OK." << std::endl;
            }
        }
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
