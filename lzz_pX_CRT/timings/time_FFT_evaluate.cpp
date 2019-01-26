#include <NTL/lzz_pX.h>
#include <NTL/vector.h>
#include <iomanip>

#include "util.h"
#include "lzz_pX_CRT.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* FFT evaluation and interpolation                           */
/*------------------------------------------------------------*/
void check()
{
    zz_p::FFTInit(0);

    for (long j = 1; j < 300; j++)
    {
        zz_pX f;
        zz_pX_Multipoint_FFT ev;
        Vec<zz_p> val;

        ev = get_FFT_points(j);

        f = random_zz_pX(j);

        cout << j << " ";
        double t;
        long nb;
        const double thresh = 0.02;

        // FFT evaluate
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

        // FFT interpolate
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
    check();
    return 0;
}
