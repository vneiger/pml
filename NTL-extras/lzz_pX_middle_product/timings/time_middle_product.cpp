#include <NTL/lzz_pX.h>
#include <iomanip>

#include "util.h"
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

    for (long dA = 0; dA < 300; dA += 5)
    {
        for (long dB = 0; dB < 300; dB += 5)
        {
            zz_pX a, b, c;
            
            cout << p << " ";
            cout << dA << " " << dB << " ";
            double t;
            long nb;
            const double thresh = 0.02;

            a = random_zz_pX(dA + 1);
            c = random_zz_pX(dA + dB + 1);
            

            nb = 0;
            t = get_time();
            do
            {
                b = middle_product(a, c, dA, dB);
                nb++;
            }
            while ( (get_time() - t) < thresh);
            t = (get_time() - t) / nb;
            cout << t << " ";
            
            nb = 0;
            t = get_time();
            do
            {
                c = a * b;
                nb++;
            }
            while ( (get_time() - t) < thresh);
            t = (get_time() - t) / nb;
            cout << t << " ";

            cout << endl;
        }
    }
}


/*------------------------------------------------------------*/
/* main just calls check()                                    */
/*------------------------------------------------------------*/
int main(int argc, char** argv)
{
    std::cout << std::fixed;
    std::cout << std::setprecision(8);
    warmup(); 
    check(0);
    check(23068673);
    check(288230376151711813);
    return 0;
}
