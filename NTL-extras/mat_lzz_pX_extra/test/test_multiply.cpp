#include <NTL/lzz_pX.h>
#include <NTL/matrix.h>
#include <NTL/vector.h>

#include "magma_output.h"
#include "mat_lzz_pX_extra.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* initializes a zz_pX_Multipoint                             */
/*------------------------------------------------------------*/
void check(){

  long i = 1000;
  long j = 2;
  Mat<zz_pX> a, b, c1, c2;

  cout << i << " " << j << " ";

  double t;
  zz_p::init(1125899906842679);
  zz_p::init(23068673);
  zz_p::UserFFTInit(23068673); // TODO: detect this

  random_mat_zz_pX(a, i, i, j);
  random_mat_zz_pX(b, i, i, j);

  t = GetTime();
  multiply_naive(c1, a, b);
  cout << GetTime()-t << endl;

  t = GetTime();
  multiply_evaluate_geometric(c2, a, b);
  cout << GetTime()-t << endl;

  if (c1 != c2)
    cout << "mismatch 1\n";

  zz_p::FFTInit(1);

  random_mat_zz_pX(a, i, i, j);
  random_mat_zz_pX(b, i, i, j);

  multiply_naive(c1, a, b);

  t = GetTime();
  multiply_evaluate_FFT(c2, a, b);
  cout << GetTime()-t << endl;

  if (c1 != c2)
    cout << "mismatch 2\n";
}  

int main(int argc, char ** argv){
  check();
  return 0;
}
