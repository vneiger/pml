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

int main(int argc, char *argv[]) {

	if (argc!=5) {
		throw "Usage: ./test_appbas prime rdim cdim deg";
	}

	long prime = atoi(argv[1]);
	long rdim = atoi(argv[2]);
	long cdim = atoi(argv[3]);
	long degree = atoi(argv[4]);

  zz_p::init(prime);

	// build random matrix
  Mat<zz_pX> pmat;
  random_mat_zz_pX(pmat, rdim, cdim, degree);

	return 0;
}
