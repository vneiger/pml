#include <NTL/lzz_pX.h>
#include <NTL/matrix.h>
#include <NTL/vector.h>
#include <iomanip>
#include <vector>
#include <numeric>
#include <algorithm>
#include <random>
#include <NTL/BasicThreadPool.h>

#include "mat_lzz_pX_extra.h"

NTL_CLIENT


/********************************************
*  tests the approximant basis algorithms  *
********************************************/

int main(int argc, char *argv[])
{
	SetNumThreads(4);

	bool verify=false;

	if (argc!=5 && argc!=6)
		throw std::invalid_argument("Usage: ./test_appbas_varmbasis rdim cdim order nbits (verify)");

	long rdim = atoi(argv[1]);
	long cdim = atoi(argv[2]);
	long order = atoi(argv[3]);
	long nbits = atoi(argv[4]);
	if (argc==6)
		verify = (atoi(argv[5])==1);

	std::vector<long> shift(rdim,0);
	//std::vector<long> shift {0,1,0,1};
	//std::vector<long> shift {4,1,0,1};
	//std::iota(shift.begin(), shift.end(),0);
	//std::shuffle(shift.begin(), shift.end(), std::mt19937{std::random_device{}()});

	if (nbits==0)
		zz_p::FFTInit(0);
	else
		zz_p::init(NTL::GenPrime_long(nbits));

	std::cout << "Testing approximant basis (mbasis variant) with random input matrix" << std::endl;
	std::cout << "--prime =\t" << zz_p::modulus();
	if (nbits==0) std::cout << "  (FFT prime)";
	std::cout << std::endl;
	std::cout << "--rdim =\t" << rdim << std::endl;
	std::cout << "--cdim =\t" << cdim << std::endl;
	std::cout << "--order =\t" << order << std::endl;
	std::cout << "--shift =\t";
	if (shift.size()<50)
		std::cout << shift << std::endl; 
	else
		std::cout << "length " << shift.size() << std::endl;

	double t1,t2,t1w,t2w;

	// build random matrix
  Mat<zz_pX> pmat;
	Vec<Mat<zz_p>> matp;
	matp.SetLength(order);
	t1w = GetWallTime(); t1 = GetTime();
  random_mat_zz_pX(pmat, rdim, cdim, order);
	for (long d = 0; d < order; ++d)
		matp[d] = coeff(pmat,d);
	t2w = GetWallTime(); t2 = GetTime();

	std::cout << "Time(random mat creation; transfer): " << (t2w-t1w) << "s,  " << (t2-t1) << "s\n";

	Mat<zz_p> kerbas;
	std::vector<long> pivdeg;
	// to warm up
	long warm_time = ceil(500000000/(rdim*cdim*cdim));
	std::cout << "warming up..." << std::endl;
	for (long i = 0; i < warm_time; ++i)
		pivdeg = popov_mbasis1(kerbas,coeff(pmat,0),shift);

	// mbasis Mat<zz_pX> version
	{
		std::cout << "~~~Testing mbasis - Mat<zz_pX> input only ~~~" << std::endl;
		t1w = GetWallTime(); t1 = GetTime();
		Mat<zz_pX> appbas;
		pivdeg = mbasis(appbas,pmat,order,shift);
		t2w = GetWallTime(); t2 = GetTime();

		std::cout << "Time(mbasis computation): " << (t2w-t1w) << "s,  " << (t2-t1) << "s\n";

		if (verify)
		{
			std::cout << "Verifying ordered weak Popov approximant basis..." << std::endl;
			t1w = GetWallTime(); t1 = GetTime();
			bool verif = is_approximant_basis(appbas,pmat,order,shift,ORD_WEAK_POPOV,true,false);
			t2w = GetWallTime(); t2 = GetTime();
			std::cout << (verif?"correct":"wrong") << std::endl;
			std::cout << "Time(verification): " << (t2w-t1w) << "s,  " << (t2-t1) << "s\n";

			if (rdim*cdim*order < 100)
			{
				std::cout << "Print output approx basis..." << std::endl;
				std::cout << appbas << std::endl;
				std::cout << "Print final residual..." << std::endl;
				Mat<zz_pX> residual;
				multiply_naive(residual,appbas,pmat);
				std::cout << residual << std::endl;
			}

			if (std::max(rdim,cdim)<33) {
				Mat<long> degmat;
				degree_matrix(degmat,appbas,shift,true);
				std::cout << "Print degree matrix of approx basis..." << std::endl;
				std::cout << degmat << std::endl;
			}
		}
	}

	// mbasis Vec<Mat<zz_p>> version input
	{
		std::cout << "~~~Testing mbasis - Vec<Mat<zz_p>> ~~~" << std::endl;
		t1w = GetWallTime(); t1 = GetTime();
		Mat<zz_pX> appbas;
		pivdeg = mbasis(appbas,matp,order,shift);
		t2w = GetWallTime(); t2 = GetTime();

		std::cout << "Time(mbasis computation): " << (t2w-t1w) << "s,  " << (t2-t1) << "s\n";

		if (verify)
		{
			std::cout << "Verifying ordered weak Popov approximant basis..." << std::endl;
			t1w = GetWallTime(); t1 = GetTime();
			bool verif = is_approximant_basis(appbas,pmat,order,shift,ORD_WEAK_POPOV,true,false);
			t2w = GetWallTime(); t2 = GetTime();
			std::cout << (verif?"correct":"wrong") << std::endl;
			std::cout << "Time(verification): " << (t2w-t1w) << "s,  " << (t2-t1) << "s\n";

			if (rdim*cdim*order < 100)
			{
				std::cout << "Print output approx basis..." << std::endl;
				std::cout << appbas << std::endl;
				std::cout << "Print final residual..." << std::endl;
				Mat<zz_pX> residual;
				multiply_naive(residual,appbas,pmat);
				std::cout << residual << std::endl;
			}

			if (std::max(rdim,cdim)<33) {
				Mat<long> degmat;
				degree_matrix(degmat,appbas,shift,true);
				std::cout << "Print degree matrix of approx basis..." << std::endl;
				std::cout << degmat << std::endl;
			}
		}
	}

	// mbasis Vec<Mat<zz_p>> version input&output
	{
		std::cout << "~~~Testing mbasis - Vec<Mat<zz_p>> input+output ~~~" << std::endl;
		t1w = GetWallTime(); t1 = GetTime();
		Vec<Mat<zz_p>> appbas;
		pivdeg = mbasis(appbas,matp,order,shift);
		t2w = GetWallTime(); t2 = GetTime();

		std::cout << "Time(mbasis computation): " << (t2w-t1w) << "s,  " << (t2-t1) << "s\n";

		if (verify)
		{
			std::cout << "Converting basis to polmat format..." << std::endl;
			Mat<zz_pX> appbas2;
			appbas2.SetDims(appbas[0].NumRows(),appbas[0].NumCols());
			for (long d = 0; d < appbas.length(); ++d) {
				for (long i = 0; i < appbas[0].NumRows(); ++i) {
					for (long j = 0; j < appbas[0].NumCols(); ++j) {
						SetCoeff(appbas2[i][j],d,appbas[d][i][j]);
					}
				}
			}
			std::cout << "Verifying ordered weak Popov approximant basis..." << std::endl;
			t1w = GetWallTime(); t1 = GetTime();
			bool verif = is_approximant_basis(appbas2,pmat,order,shift,ORD_WEAK_POPOV,true,false);
			t2w = GetWallTime(); t2 = GetTime();
			std::cout << (verif?"correct":"wrong") << std::endl;
			std::cout << "Time(verification): " << (t2w-t1w) << "s,  " << (t2-t1) << "s\n";

			if (rdim*cdim*order < 100)
			{
				std::cout << "Print output approx basis..." << std::endl;
				std::cout << appbas2 << std::endl;
				std::cout << "Print final residual..." << std::endl;
				Mat<zz_pX> residual;
				multiply_naive(residual,appbas2,pmat);
				std::cout << residual << std::endl;
			}

			if (std::max(rdim,cdim)<33) {
				Mat<long> degmat;
				degree_matrix(degmat,appbas2,shift,true);
				std::cout << "Print degree matrix of approx basis..." << std::endl;
				std::cout << degmat << std::endl;
			}
		}
	}


	// mbasis_vector Eric's version
	{
		std::cout << "~~~Testing mbasis_vector (Eric's version) ~~~" << std::endl;
		t1w = GetWallTime(); t1 = GetTime();
		Mat<zz_pX> appbas;
		pivdeg = mbasis_vector(appbas,pmat,order,shift);
		t2w = GetWallTime(); t2 = GetTime();

		std::cout << "Time(mbasisvector computation): " << (t2w-t1w) << "s,	" << (t2-t1) << "s\n";

		if (verify)
		{
			std::cout << "Converting basis to polmat format..." << std::endl;
			std::cout << "Verifying ordered weak Popov approximant basis..." << std::endl;
			t1w = GetWallTime(); t1 = GetTime();
			bool verif = is_approximant_basis(appbas,pmat,order,shift,ORD_WEAK_POPOV,true,false);
			t2w = GetWallTime(); t2 = GetTime();
			std::cout << (verif?"correct":"wrong") << std::endl;
			std::cout << "Time(verification): " << (t2w-t1w) << "s,  " << (t2-t1) << "s\n";

			if (rdim*cdim*order < 100)
			{
				std::cout << "Print output approx basis..." << std::endl;
				std::cout << appbas << std::endl;
				std::cout << "Print final residual..." << std::endl;
				Mat<zz_pX> residual;
				multiply_naive(residual,appbas,pmat);
				std::cout << residual << std::endl;
			}

			if (std::max(rdim,cdim)<33) {
				Mat<long> degmat;
				degree_matrix(degmat,appbas,shift,true);
				std::cout << "Print degree matrix of approx basis..." << std::endl;
				std::cout << degmat << std::endl;
			}
		}
	}

	return 0;
}
