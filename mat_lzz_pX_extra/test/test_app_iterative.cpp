#include <NTL/ZZ.h>
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
        throw std::invalid_argument("Usage: ./test_app_iterative nbits verbose\n\t--nbits: 0 for FFT prime, between 3 and 60 for normal prime\n\t--verbose: 0/1");

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

    std::cout << "Testing approximant basis computation (iterative)." << std::endl;
    std::cout << "--prime =\t" << zz_p::modulus() << std::endl;

    Mat<zz_pX> appbas;

    size_t i=0;
    size_t inst=0;
    for (auto pmat = test_examples.first.begin(); pmat!= test_examples.first.end(); ++pmat, ++i)
    {
        const long rdim = pmat->NumRows();
        const long cdim = pmat->NumCols();
        const long d = deg(*pmat);
        std::vector<VecLong> orders;
        VecLong help_order;
        if (d <= 1)
            help_order = {1, 2, 5}; // recall the order must be (strictly) positive
        else
            help_order = {d/2,d+1,2*d};

        //uniform orders
        orders.push_back(VecLong(cdim, help_order[0]));
        orders.push_back(VecLong(cdim, help_order[1]));
        orders.push_back(VecLong(cdim, help_order[2]));

        // non-uniform orders
        VecLong ord(cdim);

        // increasing
        ord[0] = help_order[0];
        for (long k = 1; k < cdim; ++k)
            ord[k] = ord[k-1]+1;
        orders.push_back(ord);

        // decreasing
        ord[0] = help_order[2];
        for (long k = 1; k < cdim; ++k)
            ord[k] = std::max(ord[k-1]-1, (long)1);
        orders.push_back(ord);

        // increase-decrease
        ord[0] = help_order[0];
        for (long k = 1; k < cdim/2; ++k)
            ord[k] = ord[k-1]+1;
        for (long k = cdim/2+1; k < cdim; ++k)
            ord[k] = std::max(ord[k-1]-1, (long)1);
        orders.push_back(ord);

        // decrease-increase
        ord[0] = help_order[2];
        for (long k = 1; k < cdim/2; ++k)
            ord[k] = std::max(ord[k-1]-1, (long)1);
        for (long k = cdim/2+1; k < cdim; ++k)
            ord[k] = ord[k-1]+1;
        orders.push_back(ord);

        // random
        VectorRandomBnd(ord.size(),ord.data(),help_order[2]);
        for (long k = 0; k < cdim; ++k) // make sure all > 0
            if (ord[k]==0) ord[k]=1;
        orders.push_back(ord);

        for (VecLong order : orders)
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

                // APPBAS_ITERATIVE, order-wise
                if (verbose)
                    std::cout << "Computation appbas_iterative (order-wise)... ";

                VecLong rdeg_iter1(shift);
                appbas_iterative(appbas,*pmat,order,rdeg_iter1,true);

                if (verbose)
                    std::cout << "OK. Testing... ";
                if (not is_approximant_basis(appbas,*pmat,order,shift,ORD_WEAK_POPOV,true))
                {
                    std::cout << "Error in appbas_iterative (order-wise)." << std::endl;
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

                // APPBAS_ITERATIVE, column-wise
                if (verbose)
                    std::cout << "Computation appbas_iterative (column-wise)... ";

                VecLong rdeg_iter2(shift);
                appbas_iterative(appbas,*pmat,order,rdeg_iter2,false);

                if (verbose)
                    std::cout << "OK. Testing... ";
                if (not is_approximant_basis(appbas,*pmat,order,shift,ORD_WEAK_POPOV,true))
                {
                    std::cout << "Error in appbas_iterative (column-wise)." << std::endl;
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

                // POPOV_APPBAS_ITERATIVE, order-wise
                if (verbose)
                    std::cout << "Computation popov_appbas_iterative (order-wise)... ";

                VecLong rdeg_popov_iter1(shift);
                popov_appbas_iterative(appbas,*pmat,order,rdeg_popov_iter1,true);

                if (verbose)
                    std::cout << "OK. Testing... ";
                if (not is_approximant_basis(appbas,*pmat,order,shift,POPOV,true))
                {
                    std::cout << "Error in popov_appbas_iterative (order-wise)." << std::endl;
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

                // POPOV_APPBAS_ITERATIVE, column-wise
                if (verbose)
                    std::cout << "Computation popov_appbas_iterative (column-wise)... ";

                VecLong rdeg_popov_iter2(shift);
                popov_appbas_iterative(appbas,*pmat,order,rdeg_popov_iter2,false);

                if (verbose)
                    std::cout << "OK. Testing... ";
                if (not is_approximant_basis(appbas,*pmat,order,shift,POPOV,true))
                {
                    std::cout << "Error in popov_appbas_iterative (column-wise)." << std::endl;
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
