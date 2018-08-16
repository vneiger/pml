#include <NTL/lzz_pX.h>
#include <assert.h>

#include "lzz_pX_middle_product.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* checks middle products                                     */
/*------------------------------------------------------------*/
void check()
{
    zz_p::init(1125899906842624);
  
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
    }
}


/*------------------------------------------------------------*/
/* main just calls check()                                    */
/*------------------------------------------------------------*/
int main(int argc, char** argv)
{
    check();
    return 0;
}
