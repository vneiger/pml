#include <iomanip>
#include "mat_lzz_pX_extra.h"

NTL_CLIENT

int main()
{
    //zz_p::FFTInit(0);
    zz_p::UserFFTInit(786433); // 20 bits
    //zz_p::init(GenPrime_long(60));

    std::cout << std::fixed;
    std::cout << std::setprecision(8);

    std::cout << "Basis reduction, FFT prime (FFTInit(0))" << std::endl;
    std::cout << "dim\tdeg\tinvtrunc\tHOL\t\treconstruct\ttotal\t\towP?\tright deg?" << std::endl;

    VecLong szs = {32, 64, 128, 256, 512};
    //VecLong szs = {2, 4, 8, 16, 32, 64, 128, 256, 512};

    for (size_t i=0; i<szs.size(); ++i)
    {
        const long dim = szs[i];
        long max_deg=16384;
        if (szs[i]==128)
            max_deg=2048;
        if (szs[i]==256)
            max_deg=256;
        if (szs[i]==512)
            max_deg=64;

        for (long d=1024/szs[i]; d<=max_deg; d=2*d)
        {
            // random matrix, degree deg
            Mat<zz_pX> pmat;
            random(pmat, dim, dim, d+1);

            // random unit lower triangular, degree deg
            Mat<zz_pX> trans;
            random(trans, dim, dim, d);
            for (long i = 0; i < dim; ++i)
            {
                set(trans[i][i]);
                for (long j = i+1; j < dim; ++j)
                    clear(trans[i][j]);
            }

            // random unit upper triangular, degree deg
            Mat<zz_pX> trans2;
            random(trans2, dim, dim, d);
            for (long i = 0; i < dim; ++i)
            {
                set(trans2[i][i]);
                for (long j = 0; j < i; ++j)
                    clear(trans2[i][j]);
            }

            // somehow random unimodular matrix of degree 2d
            multiply(trans, trans2, trans);

            // artificially non-reduced pmat
            Mat<zz_pX> trans_pmat;
            multiply(trans_pmat, trans, pmat);
            std::cout << "\ndegree input:\n" << degree_matrix(trans_pmat) << std::endl;

            double t=GetWallTime();
            Mat<zz_pX> reduced;
            reduced_form_gjv(reduced, trans_pmat);
            std::cout << GetWallTime()-t << "\t";
            std::cout << "\ndegree output:\n" << degree_matrix(reduced) << std::endl;

            // tests form
            std::cout << is_row_ordered_weak_popov(reduced) << "\t";
            std::cout << (deg(reduced) == d) << std::endl;
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
