#include <NTL/mat_lzz_p.h>
#include <NTL/lzz_pX.h>

#include "vec_lzz_p_extra.h"
#include "cauchy_geometric_special.h"

NTL_CLIENT

/*----------------------------------------------------*/
/*----------------------------------------------------*/
/* build and constructor                              */
/*----------------------------------------------------*/
/*----------------------------------------------------*/

void build(cauchy_geometric_special& C, const zz_p& a, const zz_p& b, const zz_p& q, 
	   const long n, const long m){

  C.n = n;
  C.m = m;
  C.q = q;
  C.a = a;
  C.b = b;

  if (m == 0 || n == 0)
    return;

  vec_zz_p P;
  P.SetLength(n); // C.P[i] = 1/q^i
  zz_p iq = 1/q, pow_iq = to_zz_p(1);
  for (long i = 0; i < n; i++){
    P[i] = pow_iq;
    pow_iq *= iq;
  }  
  // so now pow_iq = 1/q^n
  vec_zz_p tmp, tmp_n, tmp_m, itmp_n, itmp_m;
  
  tmp.SetLength(m+n-1);
  // 1/(a-b q^(i-n+1))
  zz_p pow_iq_tmp = pow_iq;
  for (long i = 0; i < m+n-1; i++){
    pow_iq_tmp *= q;
    tmp[i] = a - b*pow_iq_tmp;
  }
  vec_zz_p Qab;
  inv(Qab, tmp);

  tmp_n.SetLength(n-1);
  for (long i = 0; i < n-1; i++){
    pow_iq *= q;
    tmp_n[i] = 1 - pow_iq;
  }
  inv(itmp_n, tmp_n);

  tmp_m.SetLength(m-1);
  pow_iq *= q;
  for (long i = 0; i < m-1; i++){
    pow_iq *= q;
    tmp_m[i] = 1 - pow_iq;
  }
  inv(itmp_m, tmp_m);

  zz_p ia = 1/a, ib = 1/b;
  vec_zz_p Qa, Qb;
  Qa.SetLength(n+m-1);
  Qb.SetLength(n+m-1);
  
  for (long i = 0; i < n-1; i++){
    Qa[i] = itmp_n[i]*ia;
    Qb[i] = itmp_n[i]*ib;
  }
  Qa[n-1] = to_zz_p(0);
  Qb[n-1] = to_zz_p(0);
  for (long i = n; i < n+m-1; i++){
    Qa[i] = itmp_m[i-n]*ia;
    Qb[i] = itmp_m[i-n]*ib;
  }

  // prepare multipliers 
  long p = zz_p::modulus();
  const mulmod_t pinv = zz_p::ModulusInverse();

  C.P_pre.SetLength(n);
  for (long i = 0; i < n; i++)
    C.P_pre[i] = PrepMulModPrecon(rep(P[i]), p, pinv);

  C.Qab_pre.SetLength(n+m-1);
  C.Qa_pre.SetLength(n+m-1);
  C.Qb_pre.SetLength(n+m-1);

  for (long i = 0; i < n+m-1; i++){
    C.Qab_pre[i] = PrepMulModPrecon(rep(Qab[i]), p, pinv);
    C.Qa_pre[i] = PrepMulModPrecon(rep(Qa[i]), p, pinv);
    C.Qb_pre[i] = PrepMulModPrecon(rep(Qb[i]), p, pinv);
  }

  // prepare reps
  C.P_r.SetLength(n);
  for (long i = 0; i < n; i++)
    C.P_r[i] = rep(P[i]);

  C.Qab_r.SetLength(n+m-1);
  C.Qa_r.SetLength(n+m-1);
  C.Qb_r.SetLength(n+m-1);

  for (long i = 0; i < n+m-1; i++){
    C.Qab_r[i] = rep(Qab[i]);
    C.Qa_r[i] = rep(Qa[i]);
    C.Qb_r[i] = rep(Qb[i]);
  }
}

cauchy_geometric_special::cauchy_geometric_special(const zz_p& a, const zz_p& b, const zz_p& q, long n, long m){
  build(*this, a, b, q, n, m);
}


