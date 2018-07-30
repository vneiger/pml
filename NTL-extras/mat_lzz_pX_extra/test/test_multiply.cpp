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

    cout << "size=" << sz << ", length=" << deg << " ";

    random_mat_zz_pX(a, sz, sz, deg);
    random_mat_zz_pX(b, sz, sz, deg);

    // naive algorithm, if the size is reasonable
    long do_naive = ((sz <= 400)  && (deg <= 40)) || ((sz <= 50) && (deg <= 200)) ||  ((sz <= 10) && (deg <= 500));
    if (do_naive)
    {
	t = GetTime();
	multiply_naive(c1, a, b);
	cout << GetTime()-t << " ";

	Mat<zz_pX> c0;
	t = GetTime();
	multiply_waksman(c0, a, b);
	if (c1 != c0)
	{
	    cout << "(waksman mismatch) ";
	}
	cout << GetTime()-t << " ";
    }
    else 
    { 
	cout << "------ ------ ";
    }

    // geometric evaluation -- should be done only if feasible
    t = GetTime();
    Mat<zz_pX> c2;
    multiply_evaluate_geometric(c2, a, b);
    cout << GetTime()-t << " ";

    if (do_naive)
    {
	if (c1 != c2)
	{
	    cout << "(evaluate mismatch) ";
	}
    }

    // naive transform -- always works
    t = GetTime();
    Mat<zz_pX> c3;
    multiply_transform_naive(c3, a, b);
    cout << GetTime()-t << " ";
    if (c3 != c2)
    {
	cout << "(transform naive mismatch) ";
    }

    if (deg == 2) // TODO: write is_karatsuba.
    {
	t = GetTime();
	Mat<zz_pX> c4;
	multiply_transform_karatsuba(c4, a, b);
	cout << GetTime()-t << " ";
	if (c4 != c2)
	{
	    cout << "(transform karatsuba mismatch) ";
	}
    }
    else 
    { 
	cout << "------ ";
    }


    if (deg == 3) // TODO: write is_montgomery
    {
	t = GetTime();
	Mat<zz_pX> c4;
	multiply_transform_montgomery3(c4, a, b);
	cout << GetTime()-t << " ";
	
	if (c4 != c2)
	{
	    cout << "(transform montgomery mismatch) ";
	}
    }
    else 
    { 
	cout << "------ ";
    }

    if (deg == 4) // TODO: write is_karatsuba4
    {
	t = GetTime();
	Mat<zz_pX> c4;
	multiply_transform_karatsuba4(c4, a, b);
	cout << GetTime()-t << " ";
	
	if (c4 != c2)
	{
	    cout << "(transform karatsuba4 mismatch) ";
	}
    }
    else 
    { 
	cout << "------ ";
    }


    {
	t = GetTime();
	Mat<zz_pX> c4;
	multiply_3_primes(c4, a, b);
	cout << GetTime()-t << " ";
	if (c4 != c2)
	{
	    cout << "(transform 3 primes mismatch) ";
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


    one_check(2, 500, 288230376151711813);
    one_check(2, 500, 23068673);
    one_check(2, 500, 0);

    one_check(10, 200, 288230376151711813);
    one_check(10, 200, 23068673);
    one_check(10, 200, 0);

    // one_check(100, 100, 288230376151711813);
    // one_check(100, 100, 23068673);
    // one_check(100, 100, 0);

    // one_check(1000, 2, 288230376151711813);
    // one_check(1000, 2, 23068673);
    // one_check(1000, 2, 0);

    // one_check(1000, 3, 288230376151711813);
    // one_check(1000, 3, 23068673);
    // one_check(1000, 3, 0);

    // one_check(1000, 4, 288230376151711813);
    // one_check(1000, 4, 23068673);
    // one_check(1000, 4, 0);
}  

/*------------------------------------------------------------*/
/* main calls check                                           */
/*------------------------------------------------------------*/
int main(int argc, char ** argv)
{
    check();
    return 0;
}
