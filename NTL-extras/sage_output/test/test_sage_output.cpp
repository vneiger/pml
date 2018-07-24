#include <NTL/lzz_pX.h>
#include <NTL/vector.h>

#include "sage_output.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* initializes output, print some vectors and polys           */
/* check takes an extra argument, not used here               */
/*------------------------------------------------------------*/
void check(int opt){

  long p = 1125899906842679;
  zz_p::init(p);

  sage_init();
  sage_init_X();

  Vec<zz_p> v;
  v.SetLength(10);
  for (long i = 0; i < v.length(); i++)
    v[i] = random_zz_p();
  sage_assign(v, "v");

  zz_pX f = random_zz_pX(10);
  sage_assign(f, "f");
  sage_output(f, "x");

}  

int main(int argc, char ** argv){
  int opt = 0;
  if (argc > 1)
    opt = atoi(argv[1]);
  check(opt);
  return 0;
}
