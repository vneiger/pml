#include <NTL/lzz_pX.h>
#include <NTL/matrix.h>
#include <NTL/vector.h>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <numeric>
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

	bool order_wise=true;

	if (argc!=5 && argc!=6)
		throw std::invalid_argument("Usage: ./test_appbas rdim cdim deg nbits (order_wise)");
	if (argc == 6)
		order_wise = atoi(argv[5]) ? true : false;

	long rdim   = atoi(argv[1]);
	long cdim   = atoi(argv[2]);
	long degree = atoi(argv[3]);
	long nbits = atoi(argv[4]);
	std::vector<long> shift(rdim,0);
	std::iota(shift.begin(), shift.end(),0);
	std::shuffle(shift.begin(), shift.end(), std::mt19937{std::random_device{}()});

	long prime = NTL::GenPrime_long(nbits);
  zz_p::init(prime);

	std::cout << "Testing approximant basis computation with random input matrix" << std::endl;
	std::cout << "--prime =\t" << prime << std::endl;
	std::cout << "--rdim =\t" << rdim << std::endl;
	std::cout << "--cdim =\t" << cdim << std::endl;
	std::cout << "--degree <\t" << degree << std::endl;
	if (rdim<40)
		std::cout << "--shift =\t" << shift << std::endl;
	else
		std::cout << "--shift: shuffled iota" << std::endl;

	double t1,t2,t1w,t2w;

	// build random matrix
  Mat<zz_pX> pmat;
	t1w = GetWallTime(); t1 = GetTime();
  random_mat_zz_pX(pmat, rdim, cdim, degree);
	t2w =  GetWallTime(); t2 = GetTime();

	std::cout << "Time(random mat creation): " << (t2w-t1w) <<  "s,  " << (t2-t1) << "s\n";

	// to warm up, have more accurate timings for small instances
	if (rdim*cdim*cdim*degree*degree < 10000000)
	{
		std::cout << "warming up..." << std::endl;
		for (long i=0; i<4; ++i)
		{
			Mat<zz_pX> appbas;
			std::vector<long> order(cdim,degree);
			std::vector<long> pivdeg = appbas_iterative(appbas,pmat,order,shift,order_wise);
		}
	}

	{
		std::cout << "~~~Iterative approximant basis~~~" << std::endl;
		t1w = GetWallTime(); t1 = GetTime();
		Mat<zz_pX> appbas;
		std::vector<long> order(cdim,degree);
		std::vector<long> pivdeg = appbas_iterative(appbas,pmat,order,shift,order_wise);
		t2w =	GetWallTime(); t2 = GetTime();

		std::cout << "Time(appbas computation): " << (t2w-t1w) << "s,	" << (t2-t1) <<	"s\n";

		std::cout << "Verifying ordered weak Popov approximant basis..." << std::endl;
		t1w = GetWallTime(); t1 = GetTime();
		bool verif = is_approximant_basis(appbas,pmat,order,shift,ORD_WEAK_POPOV,true,false);
		t2w =	GetWallTime(); t2 = GetTime();
		std::cout << (verif?"correct":"wrong") << std::endl;
		std::cout << "Time(verification): " << (t2w-t1w) << "s,	" << (t2-t1) << "s\n";

		Mat<zz_pX> residual;
		multiply_naive(residual,appbas,pmat);

		if (rdim*cdim*degree < 100)
		{
			std::cout << "Print output approx basis..." << std::endl;
			std::cout << appbas << std::endl;
			std::cout << "Print final residual..." << std::endl;
			std::cout << residual << std::endl;
		}

		if (std::max(rdim,cdim)<33) {
			Mat<long> degmat;
			degree_matrix(degmat,appbas,shift,true);
			std::cout << "Print degree matrix of approx basis..." << std::endl;
			std::cout << degmat << std::endl;
		}
	}
	{
		std::cout << "~~~Iterative Popov approximant basis~~~" << std::endl;
		t1w = GetWallTime(); t1 = GetTime();
		Mat<zz_pX> appbas;
		std::vector<long> order(cdim,degree);
		std::vector<long> pivdeg = popov_appbas_iterative(appbas,pmat,order,shift,order_wise);
		t2w =	GetWallTime(); t2 = GetTime();

		std::cout << "Time(appbas computation): " << (t2w-t1w) << "s,	" << (t2-t1) <<	"s\n";

		std::cout << "Verifying Popov approximant basis..." << std::endl;
		t1w = GetWallTime(); t1 = GetTime();
		bool verif = is_approximant_basis(appbas,pmat,order,shift,POPOV,true,false);
		t2w =	GetWallTime(); t2 = GetTime();
		std::cout << (verif?"correct":"wrong") << std::endl;
		std::cout << "Time(verification): " << (t2w-t1w) << "s,	" << (t2-t1) << "s\n";

		Mat<zz_pX> residual;
		multiply_naive(residual,appbas,pmat);

		if (rdim*cdim*degree < 100)
		{
			std::cout << "Print output approx basis..." << std::endl;
			std::cout << appbas << std::endl;
			std::cout << "Print final residual..." << std::endl;
			std::cout << residual << std::endl;
		}

		if (std::max(rdim,cdim)<33) {
			Mat<long> degmat;
			degree_matrix(degmat,appbas,shift,true);
			std::cout << "Print degree matrix of approx basis..." << std::endl;
			std::cout << degmat << std::endl;
		}
	}

	return 0;
}
