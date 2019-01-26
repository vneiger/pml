#include <NTL/lzz_pX.h>
#include <NTL/vector.h>

#include "lzz_pX_CRT.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* compares geometric evaluation to subproduct tree one       */
/*------------------------------------------------------------*/
void check(long p)
{
    if (p == 0)
        zz_p::FFTInit(0);
    else
        zz_p::init(p);

    for (long j = 1; j < 500; j++)
    {
        zz_pX f;
        zz_pX_Multipoint_General evG;
        zz_pX_Multipoint_Geometric ev;
        Vec<zz_p> q, val1, val2, valG;
        zz_p a, b;

        a = random_zz_p();
        b = random_zz_p();
        q.SetLength(j);

        // with b=random
        for (long i = 0; i < j; i++)
        {
            q[i] = b * power(a, 2*i);
        }

        evG = zz_pX_Multipoint_General(q);
        ev = zz_pX_Multipoint_Geometric(a, b, j);

        f = random_zz_pX(j);     

        evG.evaluate(valG, f);
        ev.set_FFT_evaluate();
        ev.evaluate(val1, f);
        ev.unset_FFT_evaluate();
        ev.evaluate(val2, f);

        if (valG != val1 || valG != val2) 
        {
            cerr << "error for j=" << j << endl;
            exit (-1);
        }

        f = random_zz_pX(j / 2);
        evG.evaluate(valG, f);

        ev.set_FFT_evaluate();
        ev.evaluate(val1, f);
        ev.unset_FFT_evaluate();
        ev.evaluate(val2, f);

        if (valG != val1 || valG != val2) 
        {
            cerr << "error for j=" << j << endl;
            exit (-1);
        }

        ev.prepare_degree((j / 2) - 1);
        ev.set_FFT_evaluate();
        ev.evaluate(val1, f);
        ev.unset_FFT_evaluate();
        ev.evaluate(val2, f);

        if (valG != val1 || valG != val2) 
        {
            cerr << "error for j=" << j << endl;
            exit (-1);
        }


        // with b=1
        for (long i = 0; i < j; i++)
        {
            q[i] = power(a, 2*i);
        }

        evG = zz_pX_Multipoint_General(q);
        ev = zz_pX_Multipoint_Geometric(a, j);

        f = random_zz_pX(j);     

        evG.evaluate(valG, f);
        ev.set_FFT_evaluate();
        ev.evaluate(val1, f);
        ev.unset_FFT_evaluate();
        ev.evaluate(val2, f);

        if (valG != val1 || valG != val2) 
        {
            cerr << "error for j=" << j << endl;
            exit (-1);
        }

        f = random_zz_pX(j / 2);
        evG.evaluate(valG, f);

        ev.set_FFT_evaluate();
        ev.evaluate(val1, f);
        ev.unset_FFT_evaluate();
        ev.evaluate(val2, f);

        if (valG != val1 || valG != val2) 
        {
            cerr << "error for j=" << j << endl;
            exit (-1);
        }

        ev.prepare_degree((j / 2) - 1);
        ev.set_FFT_evaluate();
        ev.evaluate(val1, f);
        ev.unset_FFT_evaluate();
        ev.evaluate(val2, f);

        if (valG != val1 || valG != val2) 
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
