#include <NTL/lzz_pX.h>
#include <NTL/matrix.h>
#include <NTL/vector.h>
#include <iomanip>

#include "util.h"
#include "magma_output.h"
#include "mat_lzz_pX_extra.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* checks one product / middle product                        */
/*------------------------------------------------------------*/
void one_check(long sz, long deg)
{
    Mat<zz_pX> a, b, bb, c;
    const double thres = 0.01;

    for (long dA = deg - 1; dA < deg + 2; dA++)
        for (long dB = deg - 1; dB < deg + 2; dB++)
        {
            double t_mid, t_direct, t_naive;
            long nb;

            random(a, sz, sz, dA + 1);
            random(c, sz, sz, dA + dB + 1);
            
            cout << sz << " " << dA << " " << dB << " ";

            t_mid = get_time();
            nb = 0;
            do
            {
                middle_product(b, a, c, dA, dB);
                nb++;
            }
            while ((get_time()-t_mid) <= thres);
            t_mid = (get_time()-t_mid) / nb;

            t_naive = get_time();
            nb = 0;
            do
            {
                multiply(bb, a, c);
                nb++;
            }
            while ((get_time()-t_naive) <= thres);
            t_naive = (get_time()-t_naive) / nb;
            cout << (t_mid / t_naive) << " ";

            random(b, sz, sz, dB + 1);
            t_direct = get_time();
            nb = 0;
            do
            {
                multiply(c, a, b);
                nb++;
            }
            while ((get_time()-t_direct) <= thres);
            t_direct = (get_time()-t_direct) / nb;

            cout << (t_mid / t_direct) << " ";

            cout << endl;
        }
}

/*------------------------------------------------------------*/
/* checks some products                                       */
/*------------------------------------------------------------*/
void check(long p)
{
    vector<long> sizes = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 150, 200, 250};
    vector<long> degrees = {15, 20, 25, 30, 40, 50, 60, 70, 80, 90, 100, 150, 200, 250};

    if (p == 0)
        zz_p::FFTInit(0);
    else
        zz_p::init(p);

    for (size_t i = 0; i < sizes.size(); i++)
        for (size_t j = 0; j < degrees.size(); j++)
            one_check(sizes[i], degrees[j]);

    
}

/*------------------------------------------------------------*/
/* main calls check                                           */
/*------------------------------------------------------------*/
int main(int argc, char ** argv)
{
    std::cout << std::fixed;
    std::cout << std::setprecision(8);
    warmup();
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
