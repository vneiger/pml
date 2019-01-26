#include <NTL/lzz_pX.h>
#include <NTL/matrix.h>
#include <NTL/vector.h>
#include <iomanip>
#include <NTL/BasicThreadPool.h>

#include "util.h"
#include "mat_lzz_pX_utils.h"
#include "mat_lzz_pX_arith.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* checks a product (sz,sz+1) x (sz+1,sz+2) in degree < deg   */
/* lhs constant                                               */
/*------------------------------------------------------------*/
void one_check(long sz, long d)
{
    Mat<zz_pX> a, a2, a3, c1, c2, c3;
    Mat<zz_p> b;

    random(a, sz, sz+1, d);
    random(b, sz+1, sz+2);

    double tt = GetWallTime();
    mul(c1, a, b);
    tt = GetWallTime()-tt;

    std::cout << sz << "\t" << d << "\t" << tt << "\t";

    a2 = a; // copya is just used to try to make the timings more reliable
    tt = GetWallTime();
    Vec<Mat<zz_p>> coeffs;
    conv(coeffs, a2);
    for (long k = 0; k < coeffs.length(); ++k)
        mul(coeffs[k], coeffs[k], b);
    conv(c2, coeffs);
    tt = GetWallTime()-tt;
    std::cout << tt << "\t";

    a3 = a; // copya is just used to try to make the timings more reliable
    tt = GetWallTime();
    Mat<zz_p> coeffmat;
    c3.SetDims(sz, sz+2);
    coeffmat.SetDims(sz*d,sz+1);
    for (long k = 0; k < d; ++k)
        for (long i = 0; i < sz; ++i)
            for (long j = 0; j < sz+1; ++j)
                coeffmat[i+k*sz][j] = coeff(a3[i][j], k);
    mul(coeffmat, coeffmat, b);
    for (long k = 0; k < d; ++k)
        for (long i = 0; i < sz; ++i)
            for (long j = 0; j < sz+2; ++j)
                SetCoeff(c3[i][j], k, coeffmat[i+k*sz][j]);
    tt = GetWallTime()-tt;
    std::cout << tt << std::endl;

    if (c1 != c2 || c1 != c3)
        LogicError("multiply mismatch");
}


/*------------------------------------------------------------*/
/* for a given prime, checks some (size, degree)               */
/*------------------------------------------------------------*/
void all_checks()
{
    std::cout << std::fixed;
    std::cout << std::setprecision(5);

    VecLong szs =
    {
        1, 2, 3, 5, 10, 20, 30, 40, 50, 70, 100
    };

    VecLong degs =
    {
        1, 2, 3, 4, 5, 10, 15, 20, 25, 50, 60, 70, 100, 150, 200, 250, 300, 400, 500, 750, 1000
    };

    std::cout << "sz\tdeg\tmul\tviaconv\tbigmat" << std::endl;
    for (size_t si = 0; si < szs.size(); si++)
        for (size_t di = 0; di < degs.size(); di++)
            one_check(szs[si], degs[di]);
}


/*------------------------------------------------------------*/
/* main calls check                                           */
/*------------------------------------------------------------*/
int main()
{
    SetNumThreads(1);
    zz_p::FFTInit(0);
    all_checks();
    zz_p::UserFFTInit(786433);
    all_checks();
    zz_p::init(288230376151711813);
    all_checks();
    zz_p::init(786433);
    all_checks();

    return 0;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
