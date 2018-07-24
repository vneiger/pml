#include <NTL/lzz_pX.h>
#include <NTL/vector.h>

#include "lzz_p_extra.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* does a multipoint evaluation                               */
/*------------------------------------------------------------*/
void check(int opt){

  long p = 101;
  zz_p::init(p);
  zz_p a;

  long ord = 100;
  element_of_order(a, ord);

  for (long i = 1; i <= ord; i++)
    cout << (power(a, i) == 1) << " ";

  cout << endl;

}  

int main(int argc, char ** argv){
  int opt = 0;
  if (argc > 1)
    opt = atoi(argv[1]);
  check(opt);
  return 0;
}
