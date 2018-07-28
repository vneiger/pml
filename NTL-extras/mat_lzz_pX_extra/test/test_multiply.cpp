#include <NTL/lzz_pX.h>
#include <NTL/matrix.h>
#include <NTL/vector.h>
#include <iomanip>

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

	// naive algorithms: if the size is reasonable
	bool do_naive = (sz <= 200)  && (deg <= 40);
	bool do_naive_transform = (sz*deg<20000);

	if (do_naive)
	{
		t = GetTime();
		multiply_naive(c1, a, b);
		cout << GetTime()-t << "\t";
	}
	else 
	{ 
		cout << "-----\t\t";
	}

	// geometric evaluation
	t = GetTime();
	Mat<zz_pX> c2;
	multiply_evaluate_geometric(c2, a, b);
	cout << GetTime()-t << "\t";

	if (do_naive && (c1 != c2))
			cout << "(evaluate mismatch) ";

	// naive transform
	if (do_naive_transform)
	{
		t = GetTime();
		Mat<zz_pX> c3;
		multiply_transform_naive(c3, a, b);
		cout << GetTime()-t << "\t";
		if (c3 != c2)
			cout << "(transform naive mismatch) ";
	}
	else
		cout << "-----------\t";

	// evaluate FFT
	if ( p==0 )
	{
		t = GetTime();
		Mat<zz_pX> c2bis;
		multiply_evaluate_FFT(c2bis, a, b);
		cout << GetTime()-t << "\t";
		if (c2bis != c2)
			cout << "(evaluate FFT mismatch) ";
	}
	else
		std::cout << "--------";

	cout << endl;
}

void one_check_smalldeg(long sz, long deg, long p)
{
	Mat<zz_pX> a, b, c1;
	double t;

	if (p == 0) // init zz_p with FFTInit()
		zz_p::FFTInit(0);
	else
		zz_p::init(p);

	random_mat_zz_pX(a, sz, sz, deg);
	random_mat_zz_pX(b, sz, sz, deg);
	cout << "size=" << sz << ", length=" << deg << " ";

	if (sz*sz*deg > 10000000000)
		throw std::invalid_argument("Calling one_check_smalldeg with large size and/or degree");

	// naive algorithms: if the size is reasonable
	bool do_naive = (sz*sz*sz*deg*deg <= 1000000000);
	if (do_naive)
	{
		t = GetTime();
		multiply_naive(c1, a, b);
		cout << GetTime()-t << " ";
	}

	Mat<zz_pX> c2;
	switch (deg) {
		case 2: // TODO: write is_karatsuba.
			t = GetTime();
			multiply_transform_karatsuba(c2, a, b);
			cout << GetTime()-t << " ";
			break;
		case 3: // TODO: write is_montgomery
			t = GetTime();
			multiply_transform_montgomery3(c2, a, b);
			cout << GetTime()-t << " ";
			break;
		case 4: // TODO: write is_karatsuba4
			t = GetTime();
			multiply_transform_karatsuba4(c2, a, b);
			cout << GetTime()-t << " ";
			break;
	}

	if (do_naive && c2 != c1)
		cout << "(transform karatsuba4 mismatch) ";

	cout << endl;
}

/*------------------------------------------------------------*/
/* checks some products                                       */
/*------------------------------------------------------------*/
void check(long sz=200, long deg=4)
{
	// TODO: detect small Fourier primes

	// small degree checks
	one_check_smalldeg(sz, 2, 1125899906842679);
	one_check_smalldeg(sz, 2, 23068673);
	one_check_smalldeg(sz, 2, 0);
	one_check_smalldeg(sz, 3, 1125899906842679);
	one_check_smalldeg(sz, 3, 23068673);
	one_check_smalldeg(sz, 3, 0);
	one_check_smalldeg(sz, 4, 1125899906842679);
	one_check_smalldeg(sz, 4, 23068673);
	one_check_smalldeg(sz, 4, 0);

	// big degrees checks
	cout << "size=" << sz << ", length=" << deg << " " << endl;
	cout << "naive\t\t" << "eval geom\t" << "naive trans\t" << "eval FFT" << endl;
	one_check(sz, deg, 1125899906842679);
	one_check(sz, deg, 23068673);
	one_check(sz, deg, 0);
}  

/*------------------------------------------------------------*/
/* main calls check                                           */
/*------------------------------------------------------------*/
int main(int argc, char ** argv)
{

	std::cout << std::fixed;
	std::cout << std::setprecision(8);

	if (argc==1)
		check();
	else if (argc==3)
		check(atoi(argv[1]),atoi(argv[2]));
	else
		throw std::invalid_argument("Usage: ./test_multiply OR ./test_multiply size degree");

	return 0;
}
