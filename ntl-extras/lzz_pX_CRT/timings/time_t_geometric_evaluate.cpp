#include <NTL/lzz_pX.h>
#include <NTL/vector.h>
#include <iomanip>

#include "util.h"
#include "lzz_pX_CRT.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* compares FFT / non-FFT and full / half degree              */
/*------------------------------------------------------------*/
void check(long p)
{
    if (p == 0)
        zz_p::FFTInit(0);
    else
        zz_p::init(p);

    for (long j = 1; j < 300; j++)
    {
        zz_pX f;
        zz_pX_Multipoint_Geometric ev;
        Vec<zz_p> val;
        zz_p a;

        a = random_zz_p();
        ev = zz_pX_Multipoint_Geometric(a, j);

        f = random_zz_pX(j);

        cout << p << " ";
        cout << j << " ";
        double t;
        long nb;
        const double thresh = 0.02;
        
        // reference FFT full
        ev.set_FFT_evaluate();
        nb = 0;
        t = get_time();
        do
        {
            ev.evaluate(val, f);
            nb++;
        }
        while ( (get_time() - t) < thresh);
        t = (get_time() - t) / nb;
        cout << t << " ";

        // FFT full
        ev.set_FFT_evaluate();
        nb = 0;
        t = get_time();
        do
        {
            ev.t_evaluate(f, val);
            nb++;
        }
        while ( (get_time() - t) < thresh);
        t = (get_time() - t) / nb;
        cout << t << " ";

        // non-FFT, full
        ev.unset_FFT_evaluate();
        nb = 0;
        t = get_time();
        do
        {
            ev.t_evaluate(f, val);
            nb++;
        }
        while ( (get_time() - t) < thresh);
        t = (get_time() - t) / nb;
        cout << t << " ";

        ev.prepare_degree((j / 2) - 1);
       
        // FFT, half
        ev.set_FFT_evaluate();
        nb = 0;
        t = get_time();
        do
        {
            ev.t_evaluate(f, val, j / 2);
            nb++;
        }
        while ( (get_time() - t) < thresh);
        t = (get_time() - t) / nb;
        cout << t << " ";

        // non-FFT, half
        ev.unset_FFT_evaluate();
        nb = 0;
        t = get_time();
        do
        {
            ev.t_evaluate(f, val, j / 2);
            nb++;
        }
        while ( (get_time() - t) < thresh);
        t = (get_time() - t) / nb;
        cout << t << " ";

        cout << endl;
    }
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
