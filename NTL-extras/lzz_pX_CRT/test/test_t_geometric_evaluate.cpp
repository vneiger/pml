#include <NTL/lzz_pX.h>
#include <NTL/vector.h>
#include <iomanip>

#include "lzz_pX_CRT.h"
#include "vec_lzz_p_extra.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* computes rnd . eval(f) using direct / tranposed geom. eval */
/*------------------------------------------------------------*/
void check(long p)
{
    if (p == 0)
        zz_p::FFTInit(0);
    else
        zz_p::init(p);

    for (long j = 1; j < 500; j++)
    {
        zz_p a, b, res1, res2;
        zz_pX f, g, h, g2, h2;
        zz_pX_Multipoint_Geometric ev;
        Vec<zz_p> valF, rnd;

        a = random_zz_p();  // yuk --- TODO: fix!!!
        b = random_zz_p();

        // with b = 1
        ev = zz_pX_Multipoint_Geometric(a, j);
        f = random_zz_pX(j);
        random_vec_zz_p(rnd, j);

        ev.evaluate(valF, f);

        long degPrep = max(0, j/2 - 1); 
        long degH = max(0, j/2 - 0); // can be up to j / 2

        // FFT
        ev.set_FFT_evaluate();
        ev.t_evaluate(g, rnd);
        ev.t_evaluate(h, rnd, degH);

        res1 = 0;
        for (long i = 0; i < j; i++)
            res1 += rnd[i] * valF[i];
        res2 = 0;
        for (long i = 0; i < j; i++)
            res2 += coeff(f, i) * coeff(g, i);
        if (res1 != res2) 
        {
            cerr << "error for j=" << j << endl;
            exit (-1);
        }

        if (trunc(g, degH) != h)
        {
            cerr << "short t_evaluate error for j=" << j << endl;
            exit (-1);
        }

        // not FFT
        ev.unset_FFT_evaluate();
        ev.t_evaluate(g2, rnd);
        ev.t_evaluate(h2, rnd, degH);
        if (g != g2 || h != h2)
        {
            cerr << "error for j=" << j << endl;
            exit (-1);
        }

        // with prepared degree
        // need: degH - 1 <= degPrep
        ev.prepare_degree(degPrep);
        ev.set_FFT_evaluate();
        ev.t_evaluate(g2, rnd);
        ev.t_evaluate(h2, rnd, degH);
        if (g != g2 || h != h2)
        {
            cerr << "error for j=" << j << endl;
            exit (-1);
        }


        // same with b = random
        ev = zz_pX_Multipoint_Geometric(a, b, j);
        f = random_zz_pX(j);
        random_vec_zz_p(rnd, j);

        ev.evaluate(valF, f);

        degPrep = max(0, j/2 - 1); 
        degH = max(0, j/2 - 0); // can be up to j / 2

        // FFT
        ev.set_FFT_evaluate();
        ev.t_evaluate(g, rnd);
        ev.t_evaluate(h, rnd, degH);

        res1 = 0;
        for (long i = 0; i < j; i++)
            res1 += rnd[i] * valF[i];
        res2 = 0;
        for (long i = 0; i < j; i++)
            res2 += coeff(f, i) * coeff(g, i);
        if (res1 != res2) 
        {
            cerr << "error for j=" << j << endl;
            exit (-1);
        }

        if (trunc(g, degH) != h)
        {
            cerr << "short t_evaluate error for j=" << j << endl;
            exit (-1);
        }

        // not FFT
        ev.unset_FFT_evaluate();
        ev.t_evaluate(g2, rnd);
        ev.t_evaluate(h2, rnd, degH);
        if (g != g2 || h != h2)
        {
            cerr << "error for j=" << j << endl;
            exit (-1);
        }

        // with prepared degree
        // need: degH - 1 <= degPrep
        ev.prepare_degree(degPrep);
        ev.set_FFT_evaluate();
        ev.t_evaluate(g2, rnd);
        ev.t_evaluate(h2, rnd, degH);
        if (g != g2 || h != h2)
        {
            cerr << "error for j=" << j << endl;
            exit (-1);
        }

    }
}  

/*------------------------------------------------------------*/
/* main calls check                                           */
/*------------------------------------------------------------*/
int main(int argc, char ** argv)
{
    check(0);
    check(1125899906842679);
    return 0;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
