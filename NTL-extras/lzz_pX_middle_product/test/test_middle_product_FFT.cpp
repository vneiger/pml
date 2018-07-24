#include <NTL/lzz_pX.h>
#include <assert.h>

#include "lzz_pX_middle_product.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* if opt = 1, runs a check                                   */
/* else, runs timings                                         */
/*------------------------------------------------------------*/
void check(int opt){
  zz_p::FFTInit(0);
  
  for (long i = 1; i < 10000; i += 1){

    zz_pX a, b, c;
    a = random_zz_pX(2*i-1);
    b = random_zz_pX(i);

    if (opt == 1){
      c = middle_product(a, b, i);
      assert (c == trunc(RightShift(a*b, i-1), i));
    }
    else{
      cout << i << " ";
      
      double t;
      
      t = GetTime();
      for (long j = 0; j < 1000; j++)
	c = middle_product(a, b, i);
      t = GetTime() - t;
      cout << t << " ";

      t = GetTime();
      for (long j = 0; j < 1000; j++)
	a = b*c;
      t = GetTime() - t;
      cout << t << " ";
      
      cout << endl;
    }
  }
}

/*------------------------------------------------------------*/
/* main just calls check()                                    */
/* if not argument is given, runs timings                     */
/* if the argument 1 is given, runs check                     */
/*------------------------------------------------------------*/
int main(int argc, char** argv){
  int opt = 0;
  if (argc > 1)
    opt = atoi(argv[1]);
  check(opt);

  return 0;
}
