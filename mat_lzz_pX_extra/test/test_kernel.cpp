#include "mat_lzz_pX_kernel.h"
#include "mat_lzz_pX_forms.h" // for degree_matrix, row_pivots
#include "test_examples.h"

NTL_CLIENT

static std::ostream &operator<<(std::ostream &out, const VecLong &s)
{
    out << "[ ";
    for (auto &i: s)
        out << i << " ";
    return out << "]";
}

//void get_mat(Mat<zz_pX> & pmat, std::vector<std::vector<std::vector<long>>> std_pmat)
//{
//    for (size_t i = 0; i < std_pmat.size(); ++i)
//        for (size_t j = 0; j < std_pmat[i].size(); ++j)
//            for (size_t k = 0; k < std_pmat[i][j].size(); ++k)
//                SetCoeff(pmat[i][j], k, std_pmat[i][j][k]);
//}

int main(int argc, char *argv[])
{
    if (argc!=3)
        throw std::invalid_argument("Usage: ./test_kernel nbits verbose\n\t--nbits: 0 for FFT prime, between 3 and 60 for normal prime\n\t--verbose: 0/1");

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
        std::cout << "Launching tests on many shifts / dimensions / types of matrices." << std::endl;
        std::cout << "This may take several minutes." << std::endl;
    }

    // build couple (test_matrices, test_shifts)
    // without large amplitude shifts (otherwise, too long)
    std::pair<std::vector<Mat<zz_pX>>, std::vector<std::vector<VecLong>>>
    test_examples = build_test_examples(false);

    std::cout << "Testing kernel basis computation." << std::endl;
    std::cout << "--prime =\t" << zz_p::modulus() << std::endl;

    VecLong pivdeg; 
    Mat<zz_pX> kerbas;

    size_t i=0;
    size_t inst=0;
    for (auto pmat = test_examples.first.begin(); pmat!= test_examples.first.end(); ++pmat, ++i)
    {
        const long rdim = pmat->NumRows();
        const long cdim = pmat->NumCols();
        const long d = deg(*pmat);

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
                std::cout << "--shift =\t" << shift << std::endl;
            }

            { // direct via approximation
                if (verbose)
                    std::cout << "Computation of the kernel via approximation... ";

                Mat<zz_pX> kerbas;
                VecLong pivind;
                VecLong rdeg(shift);
                kernel_basis_via_approximation(kerbas,pivind,*pmat,rdeg);

                if (verbose)
                    std::cout << "OK. Testing... ";

                VecLong pivind_test, pivdeg_test;
                row_pivots(pivind_test, pivdeg_test, kerbas, shift);

                if (pivind != pivind_test || not is_kernel_basis(kerbas,*pmat,shift,ORD_WEAK_POPOV,true))
                {
                    std::cout << "Error in kernel_basis_via_approximation." << std::endl;
                    std::cout << "--rdim =\t" << rdim << std::endl;
                    std::cout << "--cdim =\t" << cdim << std::endl;
                    std::cout << "--deg =\t" << d << std::endl;
                    std::cout << "--shift =\t" << shift << std::endl;
                    std::cout << "--modulus = \t" << zz_p::modulus() << std::endl;
                    std::cout << "Input: " << std::endl << *pmat << std::endl;
                    std::cout << "Kernel : " << std::endl << kerbas << std::endl;
                    std::cout << "Degree matrix: " << std::endl << degree_matrix(kerbas) << std::endl;
                    std::cout << "Pivot index: " << std::endl << pivind << std::endl;
                    return 0;
                }
                if (verbose)
                    std::cout << "OK." << std::endl;
            }

            { // direct via interpolation
                if (verbose)
                    std::cout << "Computation of the kernel via interpolation... ";

                Mat<zz_pX> kerbas;
                VecLong pivind;
                VecLong rdeg(shift);
                kernel_basis_via_interpolation(kerbas,pivind,*pmat,rdeg);

                if (verbose)
                    std::cout << "OK. Testing... ";

                VecLong pivind_test, pivdeg_test;
                row_pivots(pivind_test, pivdeg_test, kerbas, shift);

                if (pivind != pivind_test || not is_kernel_basis(kerbas,*pmat,shift,ORD_WEAK_POPOV,true))
                {
                    std::cout << "Error in kernel_basis_via_interpolation." << std::endl;
                    std::cout << "--rdim =\t" << rdim << std::endl;
                    std::cout << "--cdim =\t" << cdim << std::endl;
                    std::cout << "--deg =\t" << d << std::endl;
                    std::cout << "--shift =\t" << shift << std::endl;
                    std::cout << "--modulus = \t" << zz_p::modulus() << std::endl;
                    std::cout << "Input: " << std::endl << *pmat << std::endl;
                    std::cout << "Kernel : " << std::endl << kerbas << std::endl;
                    std::cout << "Degree matrix: " << std::endl << degree_matrix(kerbas) << std::endl;
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
                kernel_basis_zls_via_approximation(kerbas,copy_pmat,rdeg);

                if (verbose)
                    std::cout << "OK. Testing... ";

                if (not is_kernel_basis(kerbas,*pmat,shift,ORD_WEAK_POPOV,true))
                {
                    std::cout << "Error in kernel_basis_zls_via_approximation." << std::endl;
                    std::cout << "--rdim =\t" << rdim << std::endl;
                    std::cout << "--cdim =\t" << cdim << std::endl;
                    std::cout << "--deg =\t" << d << std::endl;
                    std::cout << "--shift =\t" << shift << std::endl;
                    std::cout << "--modulus = \t" << zz_p::modulus() << std::endl;
                    std::cout << "Input: " << std::endl << *pmat << std::endl;
                    std::cout << "Kernel : " << std::endl << kerbas << std::endl;
                    std::cout << "Degree matrix: " << std::endl << degree_matrix(kerbas) << std::endl;
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
                kernel_basis_zls_via_interpolation(kerbas,copy_pmat,rdeg);

                if (verbose)
                    std::cout << "OK. Testing... ";

                if (not is_kernel_basis(kerbas,*pmat,shift,ORD_WEAK_POPOV,true))
                {
                    std::cout << "Error in kernel_basis_zls_via_interpolation." << std::endl;
                    std::cout << "--rdim =\t" << rdim << std::endl;
                    std::cout << "--cdim =\t" << cdim << std::endl;
                    std::cout << "--deg =\t" << d << std::endl;
                    std::cout << "--shift =\t" << shift << std::endl;
                    std::cout << "--modulus = \t" << zz_p::modulus() << std::endl;
                    std::cout << "Input: " << std::endl << *pmat << std::endl;
                    std::cout << "Kernel : " << std::endl << kerbas << std::endl;
                    std::cout << "Degree matrix: " << std::endl << degree_matrix(kerbas) << std::endl;
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
