#include "mat_lzz_pX_interpolant.h"
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

    std::cout << "Testing interpolant basis computation (pmbasis)." << std::endl;
    std::cout << "--prime =\t" << zz_p::modulus() << std::endl;

    Mat<zz_pX> intbas;

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

                if (order<zz_p::modulus()) // to make sure we can find distinct points
                { // PMBASIS general points
                    if (verbose)
                        std::cout << "Computation pmbasis... ";

                    Vec<zz_p> pts(INIT_SIZE, order);
                    for (long k = 0; k < order; ++k)
                        pts[k]=k;
                    std::shuffle(pts.begin(), pts.end(), std::mt19937{std::random_device{}()});
                    VecLong rdeg_pmbasis(shift);
                    pmbasis(intbas,*pmat,pts,rdeg_pmbasis);

                    if (verbose)
                        std::cout << "OK. Testing... ";
                    zz_pX_Multipoint_General ev(pts);
                    Vec<Mat<zz_p>> evals;
                    ev.evaluate_matrix(evals, *pmat);
                    if (not is_interpolant_basis(intbas,evals,pts,shift,ORD_WEAK_POPOV,false))
                    {
                        std::cout << "Error in pmbasis ----->>>" << std::endl;
                        std::cout << "--rdim =\t" << rdim << std::endl;
                        std::cout << "--cdim =\t" << cdim << std::endl;
                        std::cout << "--deg =\t" << d << std::endl;
                        std::cout << "--order =\t" << order << std::endl;
                        std::cout << "--shift =\t" << shift << std::endl;
                        std::cout << "--modulus = \t" << zz_p::modulus() << std::endl;
                        std::cout << "\nPOINTS:\n" << pts << std::endl;
                        std::cout << "\nINPUT MATRIX:\n" << *pmat << std::endl;
                        std::cout << "\nOUTPUT MATRIX:\n" << intbas << std::endl;
                        std::cout << "<<<------- (error in pmbasis)" << std::endl;
                        return 0;
                    }
                    if (verbose)
                        std::cout << "OK." << std::endl;
                }

                if (order<zz_p::modulus()) // to make sure we can find distinct points
                { // POPOV_PMBASIS
                    if (verbose)
                        std::cout << "Computation popov_pmbasis... ";

                    VecLong rdeg_popov_pmbasis(shift);
                    Vec<zz_p> pts(INIT_SIZE, order);
                    for (long k = 0; k < order; ++k)
                        pts[k]=k;
                    std::shuffle(pts.begin(), pts.end(), std::mt19937{std::random_device{}()});
                    popov_pmbasis(intbas,*pmat,pts,rdeg_popov_pmbasis);

                    if (verbose)
                        std::cout << "OK. Testing... ";
                    zz_pX_Multipoint_General ev(pts);
                    Vec<Mat<zz_p>> evals;
                    ev.evaluate_matrix(evals, *pmat);
                    if (not is_interpolant_basis(intbas,evals,pts,shift,POPOV,false))
                    {
                        std::cout << "Error in popov_pmbasis ----->>>" << std::endl;
                        std::cout << "--rdim =\t" << rdim << std::endl;
                        std::cout << "--cdim =\t" << cdim << std::endl;
                        std::cout << "--deg =\t" << d << std::endl;
                        std::cout << "--order =\t" << order << std::endl;
                        std::cout << "--shift =\t" << shift << std::endl;
                        std::cout << "--modulus = \t" << zz_p::modulus() << std::endl;
                        std::cout << "\nPOINTS:\n" << pts << std::endl;
                        std::cout << "\nINPUT MATRIX:\n" << *pmat << std::endl;
                        std::cout << "\nOUTPUT MATRIX:\n" << intbas << std::endl;
                        std::cout << "<<<------- (error in popov_pmbasis)" << std::endl;
                        return 0;
                    }
                    if (verbose)
                        std::cout << "OK." << std::endl;
                }

                // d<order is a requirement of pmbasis_geometric
                // 2*order+1 < zz_p::modulus is to make sure we can find geometric points
                if (d<order && 2*order+1<zz_p::modulus())
                { // PMBASIS GEOMETRIC
                    if (verbose)
                        std::cout << "Computation pmbasis (geometric)... ";

                    VecLong rdeg_pmbasis_geom(shift);
                    Vec<zz_p> pts;
                    zz_p r;
                    // geometric in degree 'order' requires an element order at least 2*deg+1
                    element_of_order(r, 2*order+1); 
                    if (IsZero(r))
                        LogicError("pmbasis-geometric test failed: did not find element of sufficiently large order after 100 trials... may be working over a too small field?");

                    pmbasis_geometric(intbas,*pmat,r,order,rdeg_pmbasis_geom,pts);

                    if (verbose)
                        std::cout << "OK. Testing... ";
                    zz_pX_Multipoint_Geometric ev(r, order+1);
                    Vec<Mat<zz_p>> evals;
                    ev.evaluate_matrix(evals, *pmat);
                    if (not is_interpolant_basis(intbas,evals,pts,shift,ORD_WEAK_POPOV,false))
                    {
                        std::cout << "Error in pmbasis_geometric ----->>>" << std::endl;
                        std::cout << "--rdim =\t" << rdim << std::endl;
                        std::cout << "--cdim =\t" << cdim << std::endl;
                        std::cout << "--deg =\t" << d << std::endl;
                        std::cout << "--order =\t" << order << std::endl;
                        std::cout << "--shift =\t" << shift << std::endl;
                        std::cout << "--modulus = \t" << zz_p::modulus() << std::endl;
                        std::cout << "\nPOINTS:\n" << pts << std::endl;
                        std::cout << "\nINPUT MATRIX:\n" << *pmat << std::endl;
                        std::cout << "\nOUTPUT MATRIX:\n" << intbas << std::endl;
                        std::cout << "<<<------- (error in pmbasis_geometric)" << std::endl;
                        return 0;
                    }
                    if (verbose)
                        std::cout << "OK." << std::endl;
                }
            }
        }
    }
    std::cout << inst << " instances processed with success (ignoring instances where the number of points is too large for the field size)." << std::endl;
    return 0;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
