#include <NTL/lzz_pX.h>
#include <NTL/vector.h>
#include <iomanip>

#include "util.h"
#include "lzz_pX_CRT.h"
#include "vec_lzz_p_extra.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* compares FFT / non-FFT                                     */
/* FFT non-FFT thresholded                                    */
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
        Vec<zz_p> q, val;
        zz_p a;

        a = random_zz_p();
        ev = zz_pX_Multipoint_Geometric(a, j);

        cout << p << " ";
        cout << j << " ";
        double t;
        long nb;
        const double thresh = 0.02;

        random_vec_zz_p(val, j);

        // FFT        
        ev.set_FFT_interpolate();
        nb = 0;
        t = get_time();
        do
        {
            ev.interpolate(f, val);
            nb++;
        }
        while ( (get_time() - t) < thresh);
        t = (get_time() - t) / nb;
        cout << t << " ";

        // non-FFT
        ev.unset_FFT_interpolate();
        nb = 0;
        t = get_time();
        do
        {
            ev.interpolate(f, val);
            nb++;
        }
        while ( (get_time() - t) < thresh);
        t = (get_time() - t) / nb;
        cout << t << " ";

        // re-builds ev, using built-in thresholds
        ev = zz_pX_Multipoint_Geometric(a, j);
        nb = 0;
        t = get_time();
        do
        {
            ev.interpolate(f, val);
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
