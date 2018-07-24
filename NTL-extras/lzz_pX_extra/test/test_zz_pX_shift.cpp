#include <NTL/lzz_pX.h>
#include <assert.h>

#include "lzz_pX_extra.h"
#include "magma_output.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* if opt = 1, runs a check                                   */
/* else, runs timings                                         */
/*------------------------------------------------------------*/
void check(int opt){
  zz_p::FFTInit(0);
  
  long i = 29;
  zz_pX a, b, d;
  zz_p c = random_zz_p();
  zz_pX_shift s(c, i);
  a = random_zz_pX(i-2);
  s.shift(b, a);
  shift(d, a, c);

  assert (b == d);

  magma_init();
  magma_init_X();
  magma_assign(a, "a");
  magma_assign(b, "b");
  cout << "c:=k!" << c << ";\n";
  cout << "print b eq Evaluate(a, x+c);\n";
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
