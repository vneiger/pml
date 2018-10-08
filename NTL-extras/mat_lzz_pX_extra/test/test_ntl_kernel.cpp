#include <NTL/BasicThreadPool.h>
#include <NTL/lzz_pX.h>
#include <NTL/matrix.h>
#include <NTL/vector.h>
#include <iomanip>

#include "util.h"
#include "mat_lzz_pX_extra.h"

NTL_CLIENT

void test(long r, long c){
	cout << "testing: " << r << " " << c << endl;
	Mat<zz_p> M;
	M.SetDims(r,c);
	Mat<zz_p> tmp;
	kernel(tmp, M);
	cout << "success" << endl;
}

int main(){
	zz_p::init(13);
	Mat<zz_p> M;

	test(4096,1);
	test(4096,1000);
	test(4096,2048);
	test(4096,3000);
	test(4096,4000);
}
