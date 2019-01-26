#include <NTL/lzz_pX.h>
#include <NTL/matrix.h>
#include <NTL/vector.h>
#include <iomanip>

#include "util.h"
#include "mat_lzz_pX_utils.h"
#include "mat_lzz_pX_arith.h"
#include "mat_lzz_pX_multiply.h"
#include "mat_lzz_pX_linsolve.h"

#include "test_examples.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* checks an (sz,sz) matrix in degree < deg                   */
/*------------------------------------------------------------*/
void one_check_via_series(long sz, long deg)
{
    Mat<zz_pX> A;
    Vec<zz_pX> b, b2, u, res;
    zz_pX den;

    do
        random(A, sz, sz, deg);
    while (determinant(coeff(A,0))==0);
    random(b, sz, deg);

    linsolve_via_series(u, den, A, b);

    if (A*u != den*b)
        Error("non-zero residue in linsolve via series");
}

/*------------------------------------------------------------*/
/* for a give prime, checks some (size, degree)               */
/*------------------------------------------------------------*/
void all_checks_via_series()
{
    VecLong szs =
    {
        1, 2, 3, 5, 10, 20
    };

    VecLong degs =
    {
        1, 2, 3, 4, 5, 10, 15, 20, 25, 50, 60, 70, 100, 150, 200, 300, 400
    };

    for (size_t si = 0; si < szs.size(); si++)
        for (size_t di = 0; di < degs.size(); di++)
            one_check_via_series(szs[si], degs[di]);

    std::cout << "Solve via series: all instances processed with success" << std::endl;
}

void all_checks_via_kernel(bool verbose)
{
    // build couple (test_matrices, test_shifts)
    // Note: shifts will not be used here, so don't build the large amplitude ones
    std::pair<std::vector<Mat<zz_pX>>, std::vector<std::vector<VecLong>>>
    test_examples = build_test_examples(false);

    Vec<zz_pX> u;
    zz_pX den;

    size_t i=0;
    size_t inst=0;
    for (auto pmat = test_examples.first.begin(); pmat!= test_examples.first.end(); ++pmat, ++i)
    {
        const long rdim = pmat->NumRows();
        const long cdim = pmat->NumCols();
        const long d = deg(*pmat);

        // for each matrix, test for four vectors b:
        std::vector<Vec<zz_pX>> rhs(4);
        // first vector b: zero
        rhs[0].SetLength(rdim);
        // second vector b: random combination of columns of pmat (den==1)
        rhs[1].SetLength(rdim);
        for (long j = 0; j < cdim; ++j)
        {
            zz_pX comb;
            random(comb, max(d/2,1));
            for (long i = 0; i < rdim; ++i)
                rhs[1][i] += comb * (*pmat)[i][j];
        }
        // third vector b: random of degree deg(pmat)
        random(rhs[2], rdim, d+1);
        // fourth vector b: random of degree 3*deg(pmat)
        // (if d==-1, if pmat==0, just use the three first)
        if (d>=0)
            random(rhs[3], rdim, 3*d+1);
        else
            rhs.pop_back();

        for (size_t i_b=0; i_b<rhs.size(); ++i_b)
        {
            ++inst;
            if (verbose)
            {
                std::cout << std::endl << "instance number " << inst;
                switch (i_b)
                {
                case 0:
                    std::cout << " (rhs=zero):" << std::endl;
                    break;
                case 1:
                    std::cout << " (rhs=comb):" << std::endl;
                    break;
                case 2:
                    std::cout << " (rhs=rd_d):" << std::endl;
                    break;
                case 3:
                    std::cout << " (rhs=rd_3d):" << std::endl;
                    break;
                }
            }

            Vec<zz_pX> b = rhs[i_b];

            if (verbose)
            {
                std::cout << "--rdim =\t" << rdim << std::endl;
                std::cout << "--cdim =\t" << cdim << std::endl;
                std::cout << "--deg =\t" << d << std::endl;
            }

            if (verbose)
                std::cout << "Computation of system solution... ";

            long success = linsolve_via_kernel(u, den, *pmat, b);

            if (verbose)
                std::cout << "OK. Testing... ";

            // TODO check there is no solution if it returned 0 !!
            // TODO check irreducible solution !!
            if (success!=0)
            {
                if ((*pmat)*u != den*b)
                {
                    std::cout << "Error in linsolve_via_kernel." << std::endl;
                    std::cout << "--rdim =\t" << rdim << std::endl;
                    std::cout << "--cdim =\t" << cdim << std::endl;
                    std::cout << "--deg =\t" << d << std::endl;
                    std::cout << "--modulus = \t" << zz_p::modulus() << std::endl;
                    std::cout << "Input matrix: " << std::endl << *pmat << std::endl;
                    std::cout << "Input vector: " << std::endl << b << std::endl;
                    std::cout << "Output vector: " << std::endl << u << std::endl;
                    std::cout << "Output denominator: " << std::endl << den << std::endl;
                    return;
                }
            }
            if (verbose)
                std::cout << "OK." << std::endl;
        }
    }
    std::cout << "Solve via kernel: all " << inst << " instances processed with success." << std::endl;
}


/*------------------------------------------------------------*/
/* checks some primes                                         */
/*------------------------------------------------------------*/
void check(bool verbose)
{
    zz_p::FFTInit(0);
    zz_pX rd;
    std::cout << "--prime =\t" << zz_p::modulus() << std::endl;
    std::cout << "Testing system solving via series..." << std::endl;
    all_checks_via_series();
    std::cout << "Testing system solving via kernel basis..." << std::endl;
    all_checks_via_kernel(verbose);

    zz_p::UserFFTInit(786433);
    std::cout << "--prime =\t" << zz_p::modulus() << std::endl;
    std::cout << "Testing system solving via series..." << std::endl;
    all_checks_via_series();
    std::cout << "Testing system solving via kernel basis..." << std::endl;
    all_checks_via_kernel(verbose);

    zz_p::init(288230376151711813);
    std::cout << "--prime =\t" << zz_p::modulus() << std::endl;
    std::cout << "Testing system solving via series..." << std::endl;
    all_checks_via_series();
    std::cout << "Testing system solving via kernel basis..." << std::endl;
    all_checks_via_kernel(verbose);

    zz_p::init(786433);
    std::cout << "--prime =\t" << zz_p::modulus() << std::endl;
    std::cout << "Testing system solving via series..." << std::endl;
    all_checks_via_series();
    std::cout << "Testing system solving via kernel basis..." << std::endl;
    all_checks_via_kernel(verbose);
}  

/*------------------------------------------------------------*/
/* main calls check                                           */
/*------------------------------------------------------------*/
int main(int argc, char * argv[])
{
    if (argc>2)
        LogicError("Usage: ./test_linsolve OR ./test_linsolve verbose");

    bool verbose = (argc==2 && (atoi(argv[1])==1));
    check(verbose);

    return 0;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
