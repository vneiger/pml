#ifndef __LZZ_PX_EXTRA__H
#define __LZZ_PX_EXTRA__H

#include <NTL/lzz_pX.h>

NTL_CLIENT

/*------------------------------------------------------------*/
/* a class that does Taylor shift                             */
/*------------------------------------------------------------*/

class zz_pX_shift {
 private:
  Vec<zz_p> fact, ifact;
  zz_pX v;
  long d;
 public:  
  /*------------------------------------------------------------*/
  /* constructor inits a few arrays                             */
  /* d is an upper bound on the degrees of the inputs           */
  /*------------------------------------------------------------*/
  zz_pX_shift(const zz_p& c, long d);

  /*------------------------------------------------------------*/
  /* g = f(x+c)                                                 */
  /*------------------------------------------------------------*/
  void shift(zz_pX& g, const zz_pX& f);

};

/*------------------------------------------------------------*/
/* g = f(x+c)                                                 */
/*------------------------------------------------------------*/
void shift(zz_pX& g, const zz_pX& f, const zz_p& c);

#endif
