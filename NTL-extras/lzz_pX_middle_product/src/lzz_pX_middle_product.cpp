#include <NTL/lzz_pX.h>

#include "lzz_pX_middle_product.h"

NTL_CLIENT


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* middle product stuff                                       */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* reverse in degree < d                                      */
/*------------------------------------------------------------*/
zz_pX rev(const zz_pX& a, long d){
  zz_pX rA;
  if (deg(a) >= d)
    Error("degree too large to reverse");
  rA.rep.SetLength(d);
  for (long i = 0; i < d; i++)
    rA.rep[i] = coeff(a, d-1-i);
  rA.normalize();
  return rA;
}

/*------------------------------------------*/
/* naive transposed product of (a,c)        */
/*------------------------------------------*/
void tPlainMul2(zz_p *b, long sb, const zz_p *a, long sa, const zz_p *c, long sc){
  long p = zz_p::modulus();
  mulmod_t pinv = PrepMulMod(p);

  for (long i = 0; i < sb; i++)
    b[i].LoopHole() = 0;

  for (long i = 0; i < sa; i++) {
    long ai = rep(a[i]);
    mulmod_precon_t apinv = PrepMulModPrecon(ai, p, pinv); 
    for (long j = 0; j < sb; j++) 
      b[j].LoopHole() = AddMod(MulModPrecon(rep(c[i+j]), ai, p, apinv), rep(b[j]), p);
  }
}

/*------------------------------------------------------------*/
/* size of the extra storage for Karatsuba                    */
/*------------------------------------------------------------*/
long Kar_stk_size(long N){
  long n, hn, sp;
  n = N;
  sp = 0;
  do {
    hn = (n+1) >> 1;
    sp += (hn << 2) - 1;
    n = hn;
  } while (n >= KARX);
  return sp;
}

/*------------------------------------------------------------*/
/* helper routine for Karatsuba                               */
/*------------------------------------------------------------*/
static inline void KarSub(zz_p *T, const zz_p *b, const long sb){
   long p = zz_p::modulus();

   for (long i = 0; i < sb; i++)
     T[i].LoopHole() = SubMod(rep(T[i]), rep(b[i]), p);
}

/*------------------------------------------------------------*/
/* helper routine for Karatsuba                               */
/*------------------------------------------------------------*/
static inline void KarFoldSub(zz_p *T, const zz_p *b, const long sb, const long hsa){
   long p = zz_p::modulus();
   long m = sb - hsa;
   long i;

   for (i = 0; i < m; i++)
     T[i].LoopHole() = SubMod(rep(b[hsa+i]), rep(b[i]), p);
   for (i = m; i < hsa; i++)
     T[i].LoopHole() = SubMod(0, rep(b[i]), p);
}

/*------------------------------------------------------------*/
/* helper routine for Karatsuba                               */
/*------------------------------------------------------------*/
static inline void KarAdd2(zz_p *T, const zz_p *a, const zz_p *b, const long s){
   long p = zz_p::modulus();
   for (long i = 0; i < s; i++)
     T[i].LoopHole() = AddMod(rep(a[i]), rep(b[i]), p);
}

/*------------------------------------------------------------*/
/* helper routine for Karatsuba                               */
/*------------------------------------------------------------*/
static inline void KarAdd(zz_p *T, const zz_p *b, const long sb){
   long p = zz_p::modulus();
   for (long i = 0; i < sb; i++)
     T[i].LoopHole() = AddMod(rep(T[i]), rep(b[i]), p);
}

/*------------------------------------------------------------*/
/* Karatsuba transposed product, switches to naive            */
/*------------------------------------------------------------*/
void tKarMul_aux(zz_p *b, const long sb, const zz_p *a, const long sa, const zz_p *c, const long sc, zz_p *stk){

  if (sa < KARX || sb < KARX){
    tPlainMul2(b, sb, a, sa, c, sc);
    return;
  }

  long hsa = (sa + 1) >> 1;
  long hsb = (sb + 1) >> 1;

  // Degenerate case I
  if (sa >= sb && sb <= hsa) {
    zz_p *T = stk; 
    stk += sb;
    tKarMul_aux(b, sb, a+hsa, sa-hsa, c+hsa, sc-hsa, stk);
    tKarMul_aux(T, sb, a, hsa, c, min(sb+hsa-1,sc), stk);
    KarAdd(b, T, sb);
    return;
  }

  // Degenerate case II
  if (sa < sb && sa <= hsb) {
    tKarMul_aux(b, hsb, a, sa, c, min(sa+hsb-1,sc), stk);
    tKarMul_aux(b+hsb, sb - hsb, a, sa, c+hsb, sc-hsb, stk);
    return;
  }

  long hs;
  if (sa < sb) 
    hs = hsb; 
  else 
    hs = hsa;
  long hs2 = hs << 1;
  long hs3 = hs2 + hs;

  zz_p *T = stk;
  stk += hs3;
  zz_p *T1 = T+hs;
  zz_p *T2 = T1+hs;

  KarFoldSub(T, a, sa, hs);
  tKarMul_aux(b, hs, T, hs, c+hs, min(hs2-1,sc-hs), stk);

  long i = sc - hs;
  KarAdd2(T, c + hs, c, i);
  for(; i < min(sc, sa+sb-1-hs); ++i) 
    T[i] = c[i];
  for(; i < sa+sb-1-hs; ++i) 
    T[i] = 0;
  tKarMul_aux(b+hs, sb - hs, a+hs, sa-hs, T1, sa+sb-1-hs2, stk);

  KarSub(b+hs, b, sb-hs);
  tKarMul_aux(T2, hs, a, hs, T, hs2-1, stk);
  KarAdd(b, T2, hs);
}



/*------------------------------------------------------------*/
/* middle product via FFT                                     */
/* returns x=(a*b div x^(N-1)) mod x^N                        */
/*------------------------------------------------------------*/
void middle_FFT(zz_pX& x, const zz_pX& a, const zz_pX& b, long N){
  long k = NextPowerOfTwo(2*N-1);
  fftRep R1(INIT_SIZE, k), R2(INIT_SIZE, k);
  TofftRep(R1, a, k);
  TofftRep(R2, b, k);
  mul(R1, R1, R2);
  FromfftRep(x, R1, N-1, 2*N-2);
}



/*------------------------------------------------------------*/
/* middle product of (a,b)                                    */
/*   c has length <= 2*N-1                                    */
/*   a has length <= N                                        */
/*   x has length <= N                                        */
/* returns x=(a*c div x^(N-1)) mod x^N                        */
/*------------------------------------------------------------*/
zz_pX middle_product(const zz_pX& c, const zz_pX& a, long N){
  
  zz_pX b;

  if (c == 0 || a == 0){
    clear(b);
    return b;
  }

  if (deg(c) >= 2*N-1 || deg(a) >= N)
    Error("degree mismatch for middle product");

  if (N > NTL_zz_pX_MUL_CROSSOVER){
    middle_FFT(b, a, c, N);
    return b;
  }

  zz_p *ap = new zz_p[N];
  const long ind_a = N-(deg(a)+1);
  const zz_p *cf_a = a.rep.elts();
  for (long i = 0; i < ind_a; i++)
    ap[i] = 0;
  for (long i = ind_a; i < N; i++)
    ap[i] = cf_a[N-1-i];

  zz_p *bp = new zz_p[N];

  zz_p *cp = new zz_p[2*N-1];
  const long ind_c = deg(c)+1;
  const zz_p *cf_c = c.rep.elts();
  for (long i = 0; i < ind_c; i++)
    cp[i] = cf_c[i];
  for (long i = ind_c; i < 2*N-1; i++)
    cp[i] = 0;

  if (N < KARX){
    tPlainMul2(bp, N, ap, N, cp, 2*N-1);
  }
  else {
    long sp = Kar_stk_size(N);
    zz_p *stk = new zz_p[sp];
    tKarMul_aux(bp, N, ap, N, cp, 2*N-1, stk);
    delete[] stk;
  }

  b.rep.SetLength(N);
  zz_p *cf_b = b.rep.elts();
  for (long i = 0; i < N; i++)
    cf_b[i] = bp[i];
  
  b.normalize();

  delete[] ap;
  delete[] bp;
  delete[] cp;

  return b;
}

