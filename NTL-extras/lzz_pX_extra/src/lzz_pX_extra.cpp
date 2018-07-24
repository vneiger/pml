#include "vec_lzz_p_extra.h"
#include "lzz_pX_extra.h"

NTL_CLIENT


/*------------------------------------------------------------*/
/* constructor inits a few arrays                             */
/* d is an upper bound on the degrees of the inputs           */
/*------------------------------------------------------------*/
zz_pX_shift::zz_pX_shift(const zz_p& c, long d){
  // TODO: assert that d+1 is a unit in zz_p

  this->d = d;
  fact.SetLength(d+1);
  fact[0] = to_zz_p(1);
  for (long i = 1; i <= d; i++)
    fact[i] = i * fact[i-1];
  inv(ifact, fact);
  
  v.rep.SetLength(d+1);  
  v.rep[0] = to_zz_p(1);
  for (long i = 1; i <= d; i++)
    v.rep[i] = c * v.rep[i-1];
  for (long i = 0; i <= d; i++)
    v.rep[i] *= ifact[i];
}

/*------------------------------------------------------------*/
/* g = f(x+c)                                                 */
/*------------------------------------------------------------*/
void zz_pX_shift::shift(zz_pX& g, const zz_pX& f){

  zz_pX u, w;
  u.rep.SetLength(d+1);
  for (long i = 0; i < f.rep.length(); i++)
    u.rep[d-i] = f.rep[i] * fact[i];
  u.normalize();
  
  w = trunc(u * v, d + 1);

  g.rep.SetLength(d+1);
  
  for (long i = 0; i < w.rep.length(); i++)
    g.rep[d-i] = w.rep[i] * ifact[d-i];
  for (long i = w.rep.length(); i <= d; i++)
    g.rep[d-i] = 0;
  g.normalize();
}


/*------------------------------------------------------------*/
/* g = f(x+c)                                                 */
/*------------------------------------------------------------*/
void shift(zz_pX& g, const zz_pX& f, const zz_p& c){
  zz_pX_shift s(c, deg(f));
  s.shift(g, f);
}
