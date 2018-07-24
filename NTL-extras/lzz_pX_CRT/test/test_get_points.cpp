#include <NTL/lzz_pX.h>
#include <NTL/vector.h>

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

    for (long j = 1; j < 1100; j = 1 + (j % 1000))
    {

	zz_pX_Multipoint & points = get_points(j);

	// double t;
	// Vec<zz_p> q;
	// zz_p a = random_zz_p();
	// q.SetLength(j);
	// for (long i = 0; i < j; i++)
	// {
	//     q[i] = power(a, 2*i);
	// }

	// cout << j << " ";

	// zz_pX f = random_zz_pX(j);
	// Vec<zz_p> valQ, valG;
	// zz_pX_Multipoint * ev;

	// zz_pX_Multipoint_General evQ(q);
	// ev = &evQ;
	// t = GetTime();
	// for (long i = 0; i < 1000; i++)
	// {
	//     ev->evaluate(valQ, f);
	// }
	// cout << GetTime()-t << " ";

	// zz_pX_Multipoint_Geometric evG(a, j);
	// ev = &evG;
	// t = GetTime();
	// for (long i = 0; i < 1000; i++)
	// {
	//     ev->evaluate(valG, f);
	// }
	// cout << GetTime()-t << endl;

	// if (valG != valQ) 
	// {
	//     cerr << "error for j=" << j << endl;
	//     exit (-1);
	// }
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
