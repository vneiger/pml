#include <NTL/lzz_pX.h>
#include <NTL/matrix.h>
#include <NTL/vector.h>
#include <iomanip>
#include <vector>
#include <chrono>

#include "mat_lzz_pX_extra.h"

NTL_CLIENT

/********************************************
*  tests the approximant basis algorithms  *
********************************************/

int main(int argc, char *argv[])
{

	bool order_wise=true;

	if (argc!=5 && argc!=6)
		throw std::invalid_argument("Usage: ./test_appbas rdim cdim deg nbits (order_wise)");
	if (argc == 6)
		order_wise = atoi(argv[5]) ? true : false;

	long rdim   = atoi(argv[1]);
	long cdim   = atoi(argv[2]);
	long degree = atoi(argv[3]);
	long nbits = atoi(argv[4]);

	long prime = NTL::GenPrime_long(nbits);
  zz_p::init(prime);

	std::cout << "Testing approximant basis computation with random input matrix" << std::endl;
	std::cout << "--prime =\t" << prime << std::endl;
	std::cout << "--rdim =\t" << rdim << std::endl;
	std::cout << "--cdim =\t" << cdim << std::endl;
	std::cout << "--degree <\t" << degree << std::endl;
	std::cout << "--shift : uniform" << std::endl;

	// build random matrix
  Mat<zz_pX> pmat;
	auto start = std::chrono::system_clock::now();
  random_mat_zz_pX(pmat, rdim, cdim, degree);
	auto end = std::chrono::system_clock::now();

	std::cout << "Time(random mat creation): " <<
		(std::chrono::duration<double> (end-start)).count() << "s\n";

	start = std::chrono::system_clock::now();
	Mat<zz_pX> appbas;
	std::vector<long> order(cdim,degree);
	std::vector<long> shift(rdim,0);
	std::vector<long> rdeg = appbas_iterative(appbas,pmat,order,shift,order_wise);
	end = std::chrono::system_clock::now();

	std::cout << "Time(appbas computation): " <<
		(std::chrono::duration<double> (end-start)).count() << "s\n";

	Mat<zz_pX> prod;
	multiply_naive(prod,appbas,pmat);

	if (rdim*cdim*degree < 100)
	{
		std::cout << appbas << std::endl;
		std::cout << prod << std::endl;
	}

	if (std::max(rdim,cdim)<33) {
		Mat<long> degmat;
		degree_matrix(degmat,appbas);
		std::cout << degmat << std::endl;
	}

	return 0;
}
