#include <NTL/lzz_pX.h>
#include <NTL/matrix.h>
#include <NTL/vector.h>
#include <iomanip>

#include "util.h"
#include "mat_lzz_pX_extra.h"

PML_CLIENT

/*------------------------------------------------------------*/
/* checks a product (s,s) x (s,1) in degree < deg             */
/* checks a product (s,s) x (s,s) in degree < deg             */
/*------------------------------------------------------------*/
void one_check(long sz, long deg)
{
    Mat<zz_pX> a, b, c1, c2;
    random(a, sz, sz, deg);

    mat_lzz_pX_lmultiplier_geometric mula(a, deg-1);
    mat_lzz_pX_lmultiplier_3_primes mul3(a, deg-1);
    std::unique_ptr<mat_lzz_pX_lmultiplier> mul = get_lmultiplier(a, deg-1);
    
    random(b, sz, 1, deg);
    multiply_evaluate_geometric(c1, a, b);

    mula.multiply(c2, b);
    if (c1 != c2)
        LogicError("Error in geometric lmultiplier");
    
    if (deg < 100)
    {
        mat_lzz_pX_lmultiplier_dense muld(a, deg-1);
        muld.multiply(c2, b);
        if (c1 != c2)
            LogicError("Error in dense lmultiplier");
    }

    mul3.multiply(c2, b);
    if (c1 != c2)
        LogicError("Error in 3 primes lmultiplier");

    mul->multiply(c2, b);
    if (c1 != c2)
        LogicError("Error in get_lmultiplier");

    if (is_FFT_prime())
    {
        mat_lzz_pX_lmultiplier_FFT_direct mulFD(a, deg-1);
        mat_lzz_pX_lmultiplier_FFT_matmul mulFM(a, deg-1);

        mulFD.multiply(c2, b);
        if (c1 != c2)
            LogicError("Error in FFT direct lmultiplier");

        mulFM.multiply(c2, b);
        if (c1 != c2)
            LogicError("Error in FFT matmul lmultiplier");
    }

    random(b, sz, sz, deg);
    multiply_evaluate_geometric(c1, a, b);

    mula.multiply(c2, b);
    if (c1 != c2)
        LogicError("Error in geometric lmultiplier");

    if (deg < 100)
    {
        mat_lzz_pX_lmultiplier_dense muld(a, deg-1);
        muld.multiply(c2, b);
        if (c1 != c2)
            LogicError("Error in dense lmultiplier");
    }

    mul3.multiply(c2, b);
    if (c1 != c2)
        LogicError("Error in 3 primes lmultiplier");

    mul->multiply(c2, b);
    if (c1 != c2)
        LogicError("Error in get_lmultiplier");

    if (is_FFT_prime())
    {
        mat_lzz_pX_lmultiplier_FFT_direct mulFD(a, deg-1);
        mat_lzz_pX_lmultiplier_FFT_matmul mulFM(a, deg-1);

        mulFD.multiply(c2, b);
        if (c1 != c2)
            LogicError("Error in FFT direct lmultiplier");

        mulFM.multiply(c2, b);
        if (c1 != c2)
            LogicError("Error in FFT matmul lmultiplier");
    }
}

/*------------------------------------------------------------*/
/* for a give prime, checks some (size, degree)               */
/*------------------------------------------------------------*/
void all_checks()
{
    VecLong szs =
        {
            1, 2, 3, 5, 10, 20, 30, 80, 100, 150
        };

    VecLong degs =
        {
            1, 2, 3, 4, 5, 10, 15, 20, 25, 50, 60, 70, 100, 150, 200, 250, 300, 400
        };

    for (size_t si = 0; si < szs.size(); si++)
        for (size_t di = 0; di < degs.size(); di++)
            one_check(szs[si], degs[di]);

}

/*------------------------------------------------------------*/
/* checks some primes                                         */
/*------------------------------------------------------------*/
void check()
{
    zz_p::FFTInit(0);
    all_checks();
    zz_p::UserFFTInit(786433);
    all_checks();
    zz_p::init(288230376151711813);
    all_checks();
    zz_p::init(786433);
    all_checks();
}  

/*------------------------------------------------------------*/
/* main calls check                                           */
/*------------------------------------------------------------*/
int main(int argc, char ** argv)
{
    check();
    return 0;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
