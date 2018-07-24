#include <NTL/lzz_pX.h>
#include <NTL/vector.h>

#include "lzz_pX_CRT.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* does a multipoint evaluation                               */
/*------------------------------------------------------------*/
void check(){

    long p = 1125899906842679;
    zz_p::init(p);

    for (long j = 1; j < 10; j++)
    {
	Vec<zz_p> q;
	q.SetLength(j);
	for (long i = 0; i < j; i++)
	{
	    q[i] = random_zz_p();
	}
      
	zz_pX_Multipoint * ev;
	zz_pX_Multipoint_General evQ(q);
	ev = &evQ;
	zz_pX f = random_zz_pX(2*j);
	Vec<zz_p> val;
	val.SetLength(j);
	ev->evaluate(val, f);
	for (long i = 0; i < j; i++)
	{
	    if (val[i] != eval(f, q[i]))
	    {
		cerr << "error Multipoint_General evaluate with j=" << j << endl;
		exit(-1);
	    }
	}
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
