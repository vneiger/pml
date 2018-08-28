#include <NTL/lzz_pX.h>
#include <NTL/vector.h>
#include <iomanip>

#include "lzz_pX_CRT.h"
#include "vec_lzz_p_extra.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* computes rnd . FFT(f) using direct and tranposed FFT       */
/*------------------------------------------------------------*/
void check(){

    zz_p::FFTInit(0);

    for (long j = 1; j < 200; j+=1)
    {
        zz_p res1, res2;
        zz_pX f, g;
        zz_pX_Multipoint_FFT ev;
        Vec<zz_p> valF, rnd;

        ev = get_FFT_points(j);
        f = random_zz_pX(j);
        random_vec_zz_p(rnd, ev.length());

        ev.evaluate(valF, f);
        ev.t_evaluate(g, rnd);

        res1 = 0;
        for (long i = 0; i < ev.length(); i++)
            res1 += rnd[i] * valF[i];

        res2 = 0;
        for (long i = 0; i < ev.length(); i++)
            res2 += coeff(f, i) * coeff(g, i);

        if (res1 != res2) 
        {
            cerr << "error for j=" << j << endl;
            exit (-1);
        }
    }
}


/*------------------------------------------------------------*/
/* main calls check                                           */
/*------------------------------------------------------------*/
int main(int argc, char ** argv){
    check();
    return 0;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
