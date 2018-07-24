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

  long i = 2000;
  Mat<zz_p> a, b, c1;

  cout << i << " ";

  double t;
  zz_p::init(90011);

  a = random_mat_zz_p(i, i);
  b = random_mat_zz_p(i, i);

  t = GetTime();
  c1 = a*b;
  cout << GetTime()-t << endl;
}  

int main(int argc, char ** argv){
  check();
  return 0;
}
