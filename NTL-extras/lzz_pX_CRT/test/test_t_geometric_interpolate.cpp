#include <NTL/lzz_pX.h>
#include <NTL/vector.h>
#include <iomanip>

#include "lzz_pX_CRT.h"
#include "vec_lzz_p_extra.h"

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
	// double t;
	zz_p a;
        zz_pX g;
	zz_pX_Multipoint_Geometric ev;
        Vec<zz_p> rnd1, rnd2;

        cout << j << " ";
        cout << "(" << ev.FFT_evaluate() << ") ";
        
        a = random_zz_p();
        ev = zz_pX_Multipoint_Geometric(a, j);
        random_vec_zz_p(rnd1, j);

	// for (long i = 0; i < 1000; i++)
	// {
	    ev.t_evaluate(g, rnd1);
	// }
	// cout << GetWallTime()-t << " ";

	// t = GetWallTime();
	// for (long i = 0; i < 1000; i++)
	// {
	    ev.t_interpolate(rnd2, g);
	// }
	// cout << GetWallTime()-t << " ";

	if (rnd1 != rnd2) 
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
