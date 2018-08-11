#include <NTL/lzz_pX.h>
#include <NTL/vector.h>
#include <iomanip>

#include "lzz_pX_CRT.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* compares geometric evaluation to subproduct tree one       */
/* prints timings and aborts if results differ                */
/*------------------------------------------------------------*/
void check()
{
    long p = 1125899906842679;
    zz_p::init(p);

    for (long j = 1; j < 1100; j++)
    {
	double t;
	zz_p a, res1, res2;
        zz_pX f, g;
	zz_pX_Multipoint_General ev;
        Vec<zz_p> valF, rnd;

        ev = get_general_points(j);
        f = random_zz_pX(j);
        random(rnd, j);

        cout << j << " ";
        
	t = GetWallTime();
	for (long i = 0; i < 1000; i++)
	{
	    ev.evaluate(valF, f);
	}
	cout << GetWallTime()-t << " ";

	t = GetWallTime();
	for (long i = 0; i < 1000; i++)
	{
	    ev.t_evaluate(g, rnd);
	}
	cout << GetWallTime()-t << " ";

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

        cout << endl;
    }
}  

/*------------------------------------------------------------*/
/* main calls check                                           */
/*------------------------------------------------------------*/
int main(int argc, char ** argv)
{
    check();
    return 0;
}
