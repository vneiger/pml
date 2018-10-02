#include <NTL/lzz_pX.h>

#include "lzz_pX_middle_product.h"

NTL_CLIENT


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* middle product stuff                                       */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* DIRT: uses undocumented MulMod feature (see sp_arith.h)    */
/* (taken from NTL zz_pX.cpp                                  */
/*------------------------------------------------------------*/
static inline 
void reduce(zz_p& r, long a, long p, mulmod_t pinv)
{
   r.LoopHole() = MulMod(a, 1L, p, pinv);
}

/*------------------------------------------------------------*/
/* naive transposed product of (a,c)                          */
/* use_long = 1 -> delays reduction                           */
/*------------------------------------------------------------*/
void tPlainMul2(zz_p *b, long sb, const zz_p *a, long sa, const zz_p *c, long sc, long use_long)
{
    long p = zz_p::modulus();
    mulmod_t pinv = PrepMulMod(p);

    if (use_long == 0)
    {
        for (long i = 0; i < sb; i++)
            b[i].LoopHole() = 0;
        for (long i = 0; i < sa; i++) 
        {
            long ai = rep(a[i]);
            mulmod_precon_t apinv = PrepMulModPrecon(ai, p, pinv); 
            for (long j = 0; j < sb; j++) 
                b[j].LoopHole() = AddMod(MulModPrecon(rep(c[i+j]), ai, p, apinv), rep(b[j]), p);
        }
    }
    else
    {
        for (long j = 0; j < sb; j++) 
        {
            long accum = 0;
            for (long i = 0; i < sa; i++) 
                accum += rep(c[i+j]) * rep(a[i]);
            reduce(b[j], accum, p, pinv);
        }
    }
}

/*------------------------------------------------------------*/
/* size of the extra storage for Karatsuba                    */
/*------------------------------------------------------------*/
long Kar_stk_size(long N)
{
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
static inline void KarSub(zz_p *T, const zz_p *b, const long sb)
{
    long p = zz_p::modulus();

    for (long i = 0; i < sb; i++)
        T[i].LoopHole() = SubMod(rep(T[i]), rep(b[i]), p);
}

/*------------------------------------------------------------*/
/* helper routine for Karatsuba                               */
/*------------------------------------------------------------*/
static inline void KarFoldSub(zz_p *T, const zz_p *b, const long sb, const long hsa)
{
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
static inline void KarAdd2(zz_p *T, const zz_p *a, const zz_p *b, const long s)
{
    long p = zz_p::modulus();
    for (long i = 0; i < s; i++)
        T[i].LoopHole() = AddMod(rep(a[i]), rep(b[i]), p);
}

/*------------------------------------------------------------*/
/* helper routine for Karatsuba                               */
/*------------------------------------------------------------*/
static inline void KarAdd(zz_p *T, const zz_p *b, const long sb)
{
    long p = zz_p::modulus();
    for (long i = 0; i < sb; i++)
        T[i].LoopHole() = AddMod(rep(T[i]), rep(b[i]), p);
}

/*------------------------------------------------------------*/
/* Karatsuba transposed product, switches to naive            */
/* value of use_long forward to tPlainMul2                    */
/*------------------------------------------------------------*/
static void tKarMul_aux(zz_p *b, const long sb, const zz_p *a, const long sa, const zz_p *c, const long sc, zz_p *stk, long use_long)
{
    if (sa < KARX || sb < KARX)
    {
        tPlainMul2(b, sb, a, sa, c, sc, use_long);
        return;
    }

    long hsa = (sa + 1) >> 1;
    long hsb = (sb + 1) >> 1;

    // Degenerate case I
    if (sa >= sb && sb <= hsa) 
    {
        zz_p *T = stk; 
        stk += sb;
        tKarMul_aux(b, sb, a+hsa, sa-hsa, c+hsa, sc-hsa, stk, use_long);
        tKarMul_aux(T, sb, a, hsa, c, min(sb+hsa-1,sc), stk, use_long);
        KarAdd(b, T, sb);
        return;
    }

    // Degenerate case II
    if (sa < sb && sa <= hsb) 
    {
        tKarMul_aux(b, hsb, a, sa, c, min(sa+hsb-1,sc), stk, use_long);
        tKarMul_aux(b+hsb, sb - hsb, a, sa, c+hsb, sc-hsb, stk, use_long);
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
    tKarMul_aux(b, hs, T, hs, c+hs, min(hs2-1,sc-hs), stk, use_long);

    long i = sc - hs;
    KarAdd2(T, c + hs, c, i);
    for(; i < min(sc, sa+sb-1-hs); ++i) 
        T[i] = c[i];
    for(; i < sa+sb-1-hs; ++i) 
        T[i] = 0;
    tKarMul_aux(b+hs, sb - hs, a+hs, sa-hs, T1, sa+sb-1-hs2, stk, use_long);

    KarSub(b+hs, b, sb-hs);
    tKarMul_aux(T2, hs, a, hs, T, hs2-1, stk, use_long);
    KarAdd(b, T2, hs);
}

/*------------------------------------------------------------*/
/* Karatsuba transposed product, switches to naive            */
/*------------------------------------------------------------*/
void tKarMul_aux(zz_p *b, const long sb, const zz_p *a, const long sa, const zz_p *c, const long sc, zz_p *stk)
{
    long p = zz_p::modulus();
    long use_long = (p < NTL_SP_BOUND/KARX && p*KARX < NTL_SP_BOUND/p);
    tKarMul_aux(b, sb, a, sa, c, sc, stk, use_long);
}

/*------------------------------------------------------------*/
/* middle product via FFT                                     */
/* returns trunc( trunc(a, dA+1)*c div x^dA, dB+1 )           */
/* output can alias input                                     */
/*------------------------------------------------------------*/
void middle_FFT(zz_pX& x, const zz_pX& a, const zz_pX& c, long dA, long dB)
{
    long k = NextPowerOfTwo(dA + dB + 1);
    fftRep R1(INIT_SIZE, k), R2(INIT_SIZE, k);
    TofftRep(R1, a, k, 0, dA);
    TofftRep(R2, c, k, 0, dA + dB);
    mul(R1, R1, R2);
    FromfftRep(x, R1, dA, dA + dB);
}

/*------------------------------------------------------------*/
/* middle product of (a,c)                                    */
/* returns trunc( trunc(a, dA+1)*c div x^dA, dB+1 )           */
/* output can alias input                                     */
/*------------------------------------------------------------*/
void middle_product(zz_pX& b, const zz_pX& a, const zz_pX& c, long dA, long dB)
{

    if (c == 0 || a == 0)
    {
        clear(b);
        return;
    }

    // not sure what's best here
    if (dA > NTL_zz_pX_MUL_CROSSOVER && dB > NTL_zz_pX_MUL_CROSSOVER)
    {
        middle_FFT(b, a, c, dA, dB);
        return;
    }

    Vec<zz_p> ap;
    ap.SetLength(dA + 1);
    const long ind_a = dA - deg(a);
    const zz_p *cf_a = a.rep.elts();
    for (long i = 0; i < ind_a; i++)
        ap[i] = 0;
    for (long i = max(0, ind_a); i <= dA; i++)
        ap[i] = cf_a[dA - i];


    Vec<zz_p> bp;
    bp.SetLength(dB + 1);

    Vec<zz_p> cp;
    cp.SetLength(dA + dB + 1);
    const long ind_c = deg(c)+1;
    const zz_p *cf_c = c.rep.elts();
    for (long i = 0; i < min(ind_c, dA + dB + 1); i++)
        cp[i] = cf_c[i];
    for (long i = ind_c; i <= dA + dB; i++)
        cp[i] = 0;

    if (min(dA, dB) < KARX)
    {
        long p = zz_p::modulus();
        long use_long = (p < NTL_SP_BOUND/KARX && p*KARX < NTL_SP_BOUND/p);
        tPlainMul2(bp.elts(), dB + 1, ap.elts(), dA + 1, cp.elts(), dA + dB + 1, use_long);
    }
    else 
    {
        long sp = Kar_stk_size(1 + max(dA, dB));
        Vec<zz_p> stk;
        stk.SetLength(sp);
        tKarMul_aux(bp.elts(), dB + 1, ap.elts(), dA + 1, cp.elts(), dA + dB + 1, stk.elts());
    }

    b.rep.SetLength(dB + 1);
    zz_p *cf_b = b.rep.elts();
    for (long i = 0; i <= dB; i++)
        cf_b[i] = bp[i];

    b.normalize();
}


// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
