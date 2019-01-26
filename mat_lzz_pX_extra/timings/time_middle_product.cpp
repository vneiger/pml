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

    for (long dA = deg; dA < deg + 1; dA++)
        for (long dB = deg; dB < deg + 1; dB++)
        {
            double t_geom, t_FFT_direct, t_FFT_matmul, t_dense, t_3_primes, t_middle, t_direct, t_naive;
            long nb;

            random(a, sz, sz, dA + 1);
            random(b, sz, sz, dB + 1);
            random(c, sz, sz, dA + dB + 1);
            
            cout << sz << " " << dA << " " << dB << " ";

            if (is_FFT_prime())
            {
                t_FFT_direct = get_time();
                nb = 0;
                do
                {
                    middle_product_evaluate_FFT_direct(b, a, c, dA, dB);
                    nb++;
                }
                while ((get_time()-t_FFT_direct) <= thres);
                t_FFT_direct = (get_time()-t_FFT_direct) / nb;

                t_FFT_matmul = get_time();
                nb = 0;
                do
                {
                    middle_product_evaluate_FFT_matmul(b, a, c, dA, dB);
                    nb++;
                }
                while ((get_time()-t_FFT_matmul) <= thres);
                t_FFT_matmul = (get_time()-t_FFT_matmul) / nb;


                t_middle = get_time();
                nb = 0;
                do
                {
                    middle_product(b, a, c, dA, dB);
                    nb++;
                }
                while ((get_time()-t_middle) <= thres);
                t_middle = (get_time()-t_middle) / nb;


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


                random(c, sz, sz, dA + dB + 1);
                t_naive = get_time();
                nb = 0;
                do
                {
                    multiply(b, a, c);
                    nb++;
                }
                while ((get_time()-t_naive) <= thres);
                t_naive = (get_time()-t_naive) / nb;

                cout << t_FFT_direct << " " << t_FFT_matmul << "   " << t_middle << "   " << t_direct << " " << t_naive;
            }
            else
            {
                t_geom = get_time();
                nb = 0;
                do
                {
                    t_multiply_evaluate_geometric(b, a, c, dA, dB);
                    nb++;
                }
                while ((get_time()-t_geom) <= thres);
                t_geom = (get_time()-t_geom) / nb;
                

                t_dense = get_time();
                nb = 0;
                do
                {
                    middle_product_evaluate_dense(b, a, c, dA, dB);
                    nb++;
                }
                while ((get_time()-t_dense) <= thres);
                t_dense = (get_time()-t_dense) / nb;
                
                
                t_3_primes = get_time();
                nb = 0;
                do
                {
                    middle_product_3_primes(b, a, c, dA, dB);
                    nb++;
                }
                while ((get_time()-t_3_primes) <= thres);
                t_3_primes = (get_time()-t_3_primes) / nb;


                t_middle = get_time();
                nb = 0;
                do
                {
                    middle_product(b, a, c, dA, dB);
                    nb++;
                }
                while ((get_time()-t_middle) <= thres);
                t_middle = (get_time()-t_middle) / nb;


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


                random(c, sz, sz, dA + dB + 1);
                t_naive = get_time();
                nb = 0;
                do
                {
                    multiply(b, a, c);
                    nb++;
                }
                while ((get_time()-t_naive) <= thres);
                t_naive = (get_time()-t_naive) / nb;

                cout << t_geom << " " << t_dense << " " << t_3_primes << "   " << t_middle << "   " << t_direct << " " << t_naive;
            }
            cout << endl;
        }
}

/*------------------------------------------------------------*/
/* checks some products                                       */
/*------------------------------------------------------------*/
void check(long p)
{
    VecLong sizes =
    {
        20, 30, 50, 100, 150, 200, 300
    };

    VecLong degrees =
    {
        1, 2, 3, 4, 5, 10, 15, 20, 25, 50, 60, 70, 100, 150, 200, 250, 300
    };

    if (p == 0)
        zz_p::FFTInit(0);
    else
        zz_p::init(p);

    cout << p << endl;
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
    // check(0);
    // check(23068673);
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
