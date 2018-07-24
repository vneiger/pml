#include <NTL/lzz_pX.h>
#include <NTL/vector.h>

#include "lzz_pX_CRT.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* initializes a zz_pX_Multipoint                             */
/*------------------------------------------------------------*/
void check(){

    zz_p::FFTInit(0);

    for (long j = 98; j < 200; j+=1)
    {
	long n = 1L << (NextPowerOfTwo(j));

	Vec<zz_p> val, pts;
	zz_pX f, X;

	zz_pX_Multipoint * ev;
	zz_pX_Multipoint_FFT evQ(n);
	ev = &evQ;

	X = 0;
	SetCoeff(X, 1, 1);
	ev->evaluate(pts, X);

	f = random_zz_pX(j);
	ev->evaluate(val, f);
    
	for (long i = 0; i < n; i++)
	{
	    if (val[i] != eval(f, pts[i]))
	    {
		cerr << "error evaluate at j=" << j << endl;
		exit(-1);
	    }
	}

	zz_pX g;
	ev->interpolate(g, val);
	if (f != g)
	{
	    cerr << "error interpolate at j=" << j << endl;
	    exit(-1);
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
