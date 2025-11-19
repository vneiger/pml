#include <NTL/FFT.h>
#include <NTL/lzz_p.h>
#include <NTL/vec_lzz_p.h>
#include <iomanip>
#include <ios>
#include <limits.h>

#include <NTL/FFT_impl.h>
#include <NTL/lzz_pX.h>
#include <NTL/tools.h>

#include "util.h"

// WARNING: for this to run up to 26 bits,
// NTL was compiled with NTL_FFTMaxRoot = 26
// (instead of default 25)

PML_CLIENT

void time_fft(long nbits, long depth)
{
    zz_pInfoT *info = zz_pInfo;
    FFTPrimeInfo *p_info = info->p_info;

    const double thresh = 0.2;

    const long len = (1<<depth);

    const long inc = len/16;

    long nb;
    double t_ntl_fft, t_ntl_ifft;

    long * a = (long *) malloc(len * sizeof(long));
    long * x = (long *) malloc(len * sizeof(long));
    for (long i = 0; i < len; i++)
    {
        a[i] = random_zz_p()._zz_p__rep;
        x[i] = random_zz_p()._zz_p__rep;
    }


    cout << nbits << "\t" << depth << "\t";
    for (long xn = inc + len/2; xn <= len; xn += inc)
    {
        nb = 0;
        t_ntl_fft = GetWallTime();
        do
        {
            new_fft(x, a, depth, *p_info, xn, xn);
            nb++;
        }
        while ((GetWallTime()-t_ntl_fft) <= thresh);
        t_ntl_fft = (GetWallTime()-t_ntl_fft) / nb;
        cout << t_ntl_fft << "\t";
    }

    for (long xn = inc + len/2; xn <= len; xn += inc)
    {
        t_ntl_ifft = GetWallTime();
        nb = 0;
        do
        {
            new_ifft(x, a, depth, *p_info, xn);
            nb++;
        }
        while ((GetWallTime()-t_ntl_ifft) <= thresh);
        t_ntl_ifft = (GetWallTime()-t_ntl_ifft) / nb;
        cout << t_ntl_ifft << "\t";
    }
    std::cout << std::endl;
}

int main(int argc, char ** argv)
{
    std::cout << std::scientific;
    std::cout << std::setprecision(1);

    if (argc != 2 && argc != 3)
    {
        std::cout << "Usage: " << argv[0] << " nbits (depth)" << std::endl;
        std::cout << "exiting..." << std::endl;
        return 0;
    }

    long nbits;

    if (argc >= 2)
    {
        nbits = atoi(argv[1]);

        if (nbits==20)
            zz_p::UserFFTInit(786433); // max 2**18
        else if (nbits==24)
            zz_p::UserFFTInit(7340033); // max 2**20
        else if (nbits==27)
            zz_p::UserFFTInit(113246209);  // max 2**22
        else if (nbits==31)
            zz_p::UserFFTInit(469762049);  // max 2**26
        else if (nbits==40)
            zz_p::UserFFTInit(890333298689); // max used 2**26
        else if (nbits==50)
            zz_p::UserFFTInit(969071659057153); // max used 2**26
        else if (nbits==60)
            zz_p::UserFFTInit(826040633134153729); // max used 2**26
                                                   //std::cout << "MAXROOT" << NTL_FFTMaxRoot << std::endl;
                                                   //std::cout << zz_pInfo->MaxRoot << std::endl;
                                                   //std::cout << CalcMaxRoot(123856518830358529) << std::endl;
    }

    if (argc == 3)
    {
        const long depth = atoi(argv[2]);

        cout << "nbits\tdepth\tfft -->\t\t\t\t\t\t\t<--fft |ifft -->\t\t\t\t\t\t<--ifft" << endl;
        warmup();
        time_fft(nbits, depth);
    }

    if (argc == 2)
    {
        cout << "nbits\tdepth\tfft -->\t\t\t\t\t\t\t<--fft |ifft -->\t\t\t\t\t\t<--ifft" << endl;
        warmup();
        for (long depth = 4; depth < 25; depth++)
            time_fft(nbits, depth);
    }
    return 0;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
