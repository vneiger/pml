#include "lzz_pX_CRT.h"
#include "mat_lzz_pX_interpolant.h"
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
        throw std::invalid_argument("Usage: ./test_int_mbasis nbits verbose\n\t--nbits: 0 for FFT prime, between 3 and 60 for normal prime\n\t--verbose: 0/1");

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

    std::cout << "Testing interpolant basis computation (mbasis)." << std::endl;
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

                { // MBASIS_RESCOMP
                    if (verbose)
                        std::cout << "Computation mbasis_rescomp... ";

                    VecLong rdeg_rescomp(shift);
                    const Vec<zz_p> pts = random_vec_zz_p(order);
                    zz_pX_Multipoint_General ev(pts);
                    Vec<Mat<zz_p>> evals;
                    ev.evaluate_matrix(evals, *pmat);
                    mbasis_rescomp(intbas,evals,pts,rdeg_rescomp,0,order); // does not modify evals

                    if (verbose)
                        std::cout << "OK. Testing... ";
                    if (not is_interpolant_basis(intbas,evals,pts,shift,ORD_WEAK_POPOV,false))
                    {
                        std::cout << "Error in mbasis_rescomp." << std::endl;
                        std::cout << "--rdim =\t" << rdim << std::endl;
                        std::cout << "--cdim =\t" << cdim << std::endl;
                        std::cout << "--deg =\t" << d << std::endl;
                        std::cout << "--order =\t" << order << std::endl;
                        std::cout << "--shift =\t" << shift << std::endl;
                        std::cout << "--modulus = \t" << zz_p::modulus() << std::endl;
                        std::cout << pts << std::endl;
                        std::cout << *pmat << std::endl;
                        std::cout << intbas << std::endl;
                        return 0;
                    }
                    if (verbose)
                        std::cout << "OK." << std::endl;
                }

                {
                    // MBASIS_RESUPDATE
                    if (verbose)
                        std::cout << "Computation mbasis_resupdate... ";

                    VecLong rdeg_resupdate(shift);
                    const Vec<zz_p> pts = random_vec_zz_p(order);
                    zz_pX_Multipoint_General ev(pts);
                    Vec<Mat<zz_p>> evals;
                    ev.evaluate_matrix(evals, *pmat);
                    Vec<Mat<zz_p>> copy_evals(evals);
                    mbasis_resupdate(intbas,copy_evals,pts,rdeg_resupdate,0,order);

                    if (verbose)
                        std::cout << "OK. Testing...";
                    if (not is_interpolant_basis(intbas,evals,pts,shift,ORD_WEAK_POPOV,false))
                    {
                        std::cout << "Error in mbasis_resupdate." << std::endl;
                        std::cout << "--rdim =\t" << rdim << std::endl;
                        std::cout << "--cdim =\t" << cdim << std::endl;
                        std::cout << "--deg =\t" << d << std::endl;
                        std::cout << "--order =\t" << order << std::endl;
                        std::cout << "--shift =\t" << shift << std::endl;
                        std::cout << "--modulus = \t" << zz_p::modulus() << std::endl;
                        std::cout << pts << std::endl;
                        std::cout << *pmat << std::endl;
                        std::cout << intbas << std::endl;
                        return 0;
                    }
                    if (verbose)
                        std::cout << "OK." << std::endl;
                }

                { // POPOV_MBASIS
                    if (verbose)
                        std::cout << "Computation popov_mbasis... ";

                    VecLong rdeg_popov(shift);
                    const Vec<zz_p> pts = random_vec_zz_p(order);
                    zz_pX_Multipoint_General ev(pts);
                    Vec<Mat<zz_p>> evals;
                    ev.evaluate_matrix(evals, *pmat);
                    Vec<Mat<zz_p>> copy_evals(evals);
                    popov_mbasis(intbas,copy_evals,pts,rdeg_popov);

                    if (verbose)
                        std::cout << "OK. Testing...";
                    if (not is_interpolant_basis(intbas,evals,pts,shift,POPOV,false))
                    {
                        std::cout << "Error in popov_mbasis." << std::endl;
                        std::cout << "--rdim =\t" << rdim << std::endl;
                        std::cout << "--cdim =\t" << cdim << std::endl;
                        std::cout << "--deg =\t" << d << std::endl;
                        std::cout << "--order =\t" << order << std::endl;
                        std::cout << "--shift =\t" << shift << std::endl;
                        std::cout << "--modulus = \t" << zz_p::modulus() << std::endl;
                        std::cout << pts << std::endl;
                        std::cout << *pmat << std::endl;
                        std::cout << intbas << std::endl;
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
