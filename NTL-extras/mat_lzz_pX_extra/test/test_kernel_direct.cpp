#include "mat_lzz_pX_kernel.h"
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
    if (argc!=2)
        throw std::invalid_argument("Usage: ./test_kernel_direct nbits");

    long nbits = atoi(argv[1]);

    if (nbits==0)
        zz_p::FFTInit(0);
    else
        zz_p::init(NTL::GenPrime_long(nbits));

    // build couple (test_matrices, test_shifts)
    std::pair<std::vector<Mat<zz_pX>>, std::vector<std::vector<VecLong>>>
    test_examples = build_test_examples();

    std::cout << "Testing kernel basis computation (direct via approximation/interpolation)." << std::endl;
    std::cout << "--prime =\t" << zz_p::modulus() << std::endl;

    VecLong pivdeg; 
    Mat<zz_pX> kerbas;

    size_t i=0;
    for (auto pmat = test_examples.first.begin(); pmat!= test_examples.first.end(); ++pmat, ++i)
    {
#ifdef VERBOSE
        std::cout << i << std::endl;
#endif // VERBOSE
        
        long rdim = pmat->NumRows();
        long cdim = pmat->NumCols();
        long d = deg(*pmat);
        for (VecLong shift : test_examples.second[i])
        {
#ifdef VERBOSE
            std::cout << "--rdim =\t" << rdim << std::endl;
            std::cout << "--cdim =\t" << cdim << std::endl;
            std::cout << "--deg =\t" << d << std::endl;
            std::cout << "--shift =\t" << shift << std::endl;
#endif // VERBOSE

            { // direct via approximation
#ifdef VERBOSE
                std::cout << "Computation of the kernel via approximation... ";
#endif // VERBOSE
                Mat<zz_pX> kerbas;
                VecLong pivind;
                VecLong rdeg(shift);
                kernel_basis_via_approximation(kerbas,pivind,*pmat,rdeg);
#ifdef VERBOSE
                std::cout << "OK. Testing... ";
#endif // VERBOSE
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
#ifdef VERBOSE
                std::cout << "OK." << std::endl;
#endif // VERBOSE
            }

            if (0)
            { // zls - approximation
#ifdef VERBOSE
                std::cout << "Computation of the kernel via ZLS algorithm (using approximation)... ";
#endif // VERBOSE
                Mat<zz_pX> kerbas;
                VecLong rdeg(shift);
                VecLong pivind;
                kernel_basis_zls_via_approximation(kerbas,pivind,*pmat,rdeg);
#ifdef VERBOSE
                std::cout << "OK. Testing... ";
#endif // VERBOSE
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
#ifdef VERBOSE
                std::cout << "OK." << std::endl;
#endif // VERBOSE
            }

            if (0)
            { // zls - interpolation
#ifdef VERBOSE
                std::cout << "Computation of the kernel via ZLS algorithm (using interpolation)... ";
#endif // VERBOSE
                Mat<zz_pX> kerbas;
                pivdeg = kernel_basis_zls_via_interpolation(kerbas,*pmat,shift);
#ifdef VERBOSE
                std::cout << "OK. Testing... ";
#endif // VERBOSE
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
#ifdef VERBOSE
                std::cout << "OK." << std::endl;
#endif // VERBOSE
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
