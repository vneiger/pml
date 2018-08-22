#include <NTL/lzz_pX.h>
#include <NTL/matrix.h>
#include <NTL/vector.h>
#include <iomanip>
#include <vector>
#include <numeric>
#include <algorithm>
#include <random>
#include <NTL/BasicThreadPool.h>

#include "util.h"
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
		throw std::invalid_argument("Usage: ./test_appbas_pmbasis rdim cdim order nbits (verify)");

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
		//zz_p::FFTInit(0);
		zz_p::UserFFTInit(1769473); // --> small FFT prime like in LinBox
	else
		zz_p::init(NTL::GenPrime_long(nbits));

	std::cout << "Testing approximant basis (pmbasis) with random input matrix" << std::endl;
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
	t1w = GetWallTime(); t1 = GetTime();
  random_mat_zz_pX(pmat, rdim, cdim, order);
	t2w = GetWallTime(); t2 = GetTime();

	std::cout << "Time(random mat creation): " << (t2w-t1w) << "s,  " << (t2-t1) << "s\n";

	Mat<zz_p> kerbas;
	std::vector<long> pivdeg;

	std::cout << "warming up..." << std::endl;
	warmup();

	// pmbasis
	{
		std::cout << "~~~Testing pmbasis~~~" << std::endl;
		t1w = GetWallTime(); t1 = GetTime();
		Mat<zz_pX> appbas;
		pivdeg = pmbasis(appbas,pmat,order,shift);
		t2w = GetWallTime(); t2 = GetTime();

		std::cout << "Time(pmbasis computation): " << (t2w-t1w) << "s,  " << (t2-t1) << "s\n";

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

	// popov_pmbasis
	{
		std::cout << "~~~Testing popov_pmbasis~~~" << std::endl;
		t1w = GetWallTime(); t1 = GetTime();
		Mat<zz_pX> appbas;
		pivdeg = popov_pmbasis(appbas,pmat,order,shift);
		t2w = GetWallTime(); t2 = GetTime();

		std::cout << "Time(popov_pmbasis computation): " << (t2w-t1w) << "s,  " << (t2-t1) << "s\n";

		if (verify)
		{
			std::cout << "Verifying Popov approximant basis..." << std::endl;
			t1w = GetWallTime(); t1 = GetTime();
			bool verif = is_approximant_basis(appbas,pmat,order,shift,POPOV,true,false);
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
