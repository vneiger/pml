#include <NTL/lzz_pX.h>
#include <NTL/vector.h>
#include <iomanip>

#include "util.h"
#include "lzz_pX_CRT.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* warm-up the CPU                                            */
/*------------------------------------------------------------*/
void warmup()
{
    double t = get_time();
    do
    {
        ;
    } 
    while (get_time() - t < 1);
}

/*------------------------------------------------------------*/
/* compares geometric evaluation to subproduct tree one       */
/*------------------------------------------------------------*/
void check(long p)
{
    if (p == 0)
        zz_p::FFTInit(0);
    else
        zz_p::init(p);

    for (long j = 10; j < 300; j++)
    {
        zz_pX f;
        zz_pX_Multipoint_Geometric ev;
        Vec<zz_p> q, val1, val2, valG;
        zz_p a;

        a = random_zz_p();
        q.SetLength(j);
        for (long i = 0; i < j; i++)
        {
            q[i] = power(a, 2*i);
        }
        
        ev = zz_pX_Multipoint_Geometric(a, j);

        f = random_zz_pX(j);

        cout << p << " ";
        cout << j << " ";
        double t;
        long nb;
        const double thresh = 0.05;
        
        ev.set_FFT_evaluate();
        nb = 0;
        t = get_time();
        do
        {
            ev.evaluate(val1, f);
            nb++;
        }
        while ( (get_time() - t) < thresh);
        t = (get_time() - t) / nb;
        cout << t << " ";

        ev.unset_FFT_evaluate();
        nb = 0;
        t = get_time();
        do
        {
            ev.evaluate(val2, f);
            nb++;
        }
        while ( (get_time() - t) < thresh);
        t = (get_time() - t) / nb;
        cout << t << " ";

        f = random_zz_pX(j / 2);

        ev.set_FFT_evaluate();
        nb = 0;
        t = get_time();
        do
        {
            ev.evaluate(val1, f);
            nb++;
        }
        while ( (get_time() - t) < thresh);
        t = (get_time() - t) / nb;
        cout << t << " ";

        ev.unset_FFT_evaluate();
        nb = 0;
        t = get_time();
        do
        {
            ev.evaluate(val2, f);
            nb++;
        }
        while ( (get_time() - t) < thresh);
        t = (get_time() - t) / nb;
        cout << t << " ";

        ev.prepare_degree((j / 2) - 1);
       
        ev.set_FFT_evaluate();
        nb = 0;
        t = get_time();
        do
        {
            ev.evaluate(val1, f);
            nb++;
        }
        while ( (get_time() - t) < thresh);
        t = (get_time() - t) / nb;
        cout << t << " ";

        ev.unset_FFT_evaluate();
        nb = 0;
        t = get_time();
        do
        {
            ev.evaluate(val2, f);
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
