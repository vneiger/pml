#include <NTL/lzz_pX.h>
#include <NTL/matrix.h>
#include <NTL/vector.h>

#include "magma_output.h"
#include "mat_lzz_pX_extra.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* checks some products                                       */
/*------------------------------------------------------------*/
void one_check(long sz, long deg, long p)
{
    Mat<zz_pX> a, b, c1;
    double t;

    if (p == 0) // init zz_p with FFTInit()
    {
	zz_p::FFTInit(0);
    }
    else
    {
	zz_p::init(p);
    }

    random_mat_zz_pX(a, sz, sz, deg);
    random_mat_zz_pX(b, sz, sz, deg);
    cout << sz << " " << deg << " ";

    // naive algorithm, if the size is reasonable
    long do_naive = (sz <= 400)  && (deg <= 40);
    if (do_naive)
    {
	t = GetTime();
	multiply_naive(c1, a, b);
	cout << GetTime()-t << " ";
    }
    else 
    { 
	cout << "------ ";
    }

    // geometric evaluation
    t = GetTime();
    Mat<zz_pX> c2;
    multiply_evaluate_geometric(c2, a, b);
    cout << GetTime()-t << " ";

    if (do_naive)
    {
	if (c1 != c2)
	{
	    cout << "evaluate mismatch with p=" << p << ", sz=" << sz << ", deg=" << deg << endl;
	}
    }

    // naive transform
    t = GetTime();
    Mat<zz_pX> c3;
    multiply_transform_naive(c3, a, b);
    cout << GetTime()-t << " ";

    if (c3 != c2)
    {
	cout << "transform naive mismatch with p=" << p << ", sz=" << sz << ", deg=" << deg << endl;
    }

    if (deg == 2) // TODO: write is_karatsuba..
    {
	t = GetTime();
	Mat<zz_pX> c4;
	multiply_transform_karatsuba(c4, a, b);
	cout << GetTime()-t << " ";
	
	if (c4 != c2)
	{
	    cout << "transform karatsuba mismatch with p=" << p << ", sz=" << sz << ", deg=" << deg << endl;
	}
    }
    


    cout << endl;
}

/*------------------------------------------------------------*/
/* checks some products                                       */
/*------------------------------------------------------------*/
void check()
{
    
// TODO: detect small Fourier primes
    one_check(1000, 2, 1125899906842679);
    one_check(1000, 2, 23068673);
    one_check(1000, 2, 0);

}  

/*------------------------------------------------------------*/
/* main calls check                                           */
/*------------------------------------------------------------*/
int main(int argc, char ** argv)
{
    check();
    return 0;
}
