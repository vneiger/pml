#include "mat_lzz_pX_extra.h"

NTL_CLIENT

int main(int argc, char *argv[])
{
    //if (argc!=3)
    //    throw std::invalid_argument("Usage: ./test_reduction nbits verbose\n\t--nbits: 0 for FFT prime, between 3 and 60 for normal prime\n\t--verbose: 0/1");

    const long nbits = atoi(argv[1]);
    const long dim = atoi(argv[2]);
    const long deg = atoi(argv[3]);

    if (nbits==0)
        zz_p::FFTInit(0);
    else if (nbits==2)
        zz_p::init(3); // make sure modulus!=2
    else
        zz_p::init(NTL::GenPrime_long(nbits));

    // random matrix, degree deg
    Mat<zz_pX> pmat;
    random(pmat, dim, dim, deg);

    // random unit lower triangular, degree deg
    Mat<zz_pX> trans;
    random(trans, dim, dim, deg);
    for (long i = 0; i < dim; ++i)
    {
        set(trans[i][i]);
        for (long j = i+1; j < dim; ++j)
            clear(trans[i][j]);
    }

    // artificially non-reduced pmat
    Mat<zz_pX> trans_pmat;
    multiply(trans_pmat, trans, pmat);
    std::cout << "degree input:\n" << degree_matrix(trans_pmat) << std::endl;

    double t=GetWallTime();
    Mat<zz_pX> reduced;
    reduced_form_gjv(reduced, trans_pmat);
    std::cout << GetWallTime()-t << std::endl;

    std::cout << degree_matrix(reduced) << std::endl;
    std::cout << is_row_reduced(reduced) << std::endl;
    std::cout << is_row_ordered_weak_popov(reduced) << std::endl;

    return 0;
}


// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
