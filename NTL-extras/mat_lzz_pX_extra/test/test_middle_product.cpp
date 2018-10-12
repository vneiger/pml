#include <NTL/lzz_pX.h>
#include <NTL/matrix.h>
#include <NTL/vector.h>
#include <iomanip>

#include "util.h"
#include "mat_lzz_pX_extra.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* checks one product / middle product                        */
/*------------------------------------------------------------*/
void one_check(long sz, long dg)
{
    Mat<zz_pX> a, b1, b2, c;

    for (long dA = dg - 1; dA < dg + 2; dA++)
        for (long dB = dg - 1; dB < dg + 2; dB++)
        {
            random(a, sz, sz+1, dA + 1);
            random(c, sz+1, sz+2, dA + dB + 1);

            // testing all possible call sequences            
            middle_product(b1, a, c, dA, dB);

            multiply(b2, a, c);
            b2 >>= dA;
            trunc(b2, b2, dB + 1);
            if (b1 != b2)
            {
                LogicError("Error in middle product");
            }

            middle_product_3_primes(b2, a, c, dA, dB);
            if (b1 != b2)
            {
                LogicError("Error in 3 primes middle product");
            }

            if (is_FFT_prime())
            {
                middle_product_FFT(b2, a, c, dA, dB);
                if (b1 != b2)
                {
                    LogicError("Error in FFT middle product");
                }
            }

            if (sz <= 10 || max(dA, dB) <= 8)
            {
                middle_product_naive(b2, a, c, dA, dB);
                if (b1 != b2)
                {
                    LogicError("Error in FFT middle product");
                }
            }
        }
}

/*------------------------------------------------------------*/
/* for a give prime, checks some (size, dgree)               */
/*------------------------------------------------------------*/
void all_checks()
{
    std::vector<long> szs =
    {
        1, 2, 3, 5, 10, 20, 30
    };

    std::vector<long> dgs =
    {
        1, 2, 3, 4, 5, 10, 15, 20, 25, 50, 60, 70, 100, 150, 200, 250, 300, 400
    };

    for (size_t si = 0; si < szs.size(); si++)
        for (size_t di = 0; di < dgs.size(); di++)
            one_check(szs[si], dgs[di]);
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
