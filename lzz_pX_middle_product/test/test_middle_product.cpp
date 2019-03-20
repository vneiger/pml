#include <NTL/lzz_pX.h>
#include <assert.h>

#include "lzz_pX_middle_product.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* checks middle products                                     */
/*------------------------------------------------------------*/
void check(long p)
{
    if (p == 0)
        zz_p::FFTInit(0);
    else
        zz_p::init(p);

    for (long dA = 0; dA < 300; dA++)
    {
        for (long dB = 0; dB < 300; dB++)
        {
            zz_pX a, b, c;
            a = random_zz_pX(dA + 1);
            c = random_zz_pX(dA + dB + 1);

            b = middle_product(a, c, dA, dB);
            if (b != trunc(RightShift(a*c, dA), dB + 1))
            {
                LogicError("Error in middle product.");
            }
        }
        {
            zz_pX a, b, c;
            a = random_zz_pX(dA + 1);
            c = random_zz_pX(2*dA + 1);
            b = middle_product(a, c, dA + 1);
            if (b != trunc(RightShift(a*c, dA), dA + 1))
            {
                LogicError("Error in middle product.");
            }
        }
        if (p==0 && dA > 3)
        {
            std::cout << dA << "\t" << "entering variant..." << std::endl;
            for (long dB = 3; dB < 300; dB++)
            {
                std::cout << dB << std::endl;
                zz_pX a, b, c;
                a = random_zz_pX(dA + 1);
                //c = random_zz_pX(dA + dB + 1);
                c = random_zz_pX(2*dA + 1);

                middle_FFT_variant(b, a, c, dA, dA);
                if (b != trunc(RightShift(a*c, dA), dA + 1))
                {
                    std::cout << b << std::endl;
                    std::cout << trunc(RightShift(a*c, dA), dA + 1) << std::endl;
                    LogicError("Error in middle product.");
                }
            }
        }
    }
}


/*------------------------------------------------------------*/
/* main just calls check()                                    */
/*------------------------------------------------------------*/
int main()
{
    check(0);
    check(23068673);
    check(288230376151711813);
    return 0;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
