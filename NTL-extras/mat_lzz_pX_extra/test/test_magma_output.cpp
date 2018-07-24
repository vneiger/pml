#include <NTL/lzz_pX.h>
#include <NTL/matrix.h>
#include <NTL/vector.h>

#include "magma_output.h"
#include "mat_lzz_pX_extra.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* initializes a zz_pX_Multipoint                             */
/* check takes an extra argument, not used here               */
/*------------------------------------------------------------*/
void check(int opt){

  zz_p::init(1125899906842679);
  Mat<zz_pX> a;
  random_mat_zz_pX(a, 3, 4, 2);
  magma_init();
  magma_init_X();
  magma_assign(a, "x", "a");

}  

int main(int argc, char ** argv){
  int opt = 0;
  if (argc > 1)
    opt = atoi(argv[1]);
  check(opt);
  return 0;
}
