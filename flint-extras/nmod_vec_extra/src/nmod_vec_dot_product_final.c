#include <flint/nmod.h>

#include "nmod_extra.h"
#include "nmod_vec_extra.h"

ulong nmod_vec_dot_product_v2(nn_srcptr v1, nn_srcptr v2, ulong len, nmod_t mod, ulong n_limbs)
{
    if (n_limbs == 1)
    {
        ulong res = UWORD(0);

        for (ulong i = 0; i < len; i++)
            res += v1[i] * v2[i];

        NMOD_RED(res, res, mod);
        return res;
    }

    if (n_limbs == 2)
    {
        ulong s0, s1;
        ulong u0 = UWORD(0);
        ulong u1 = UWORD(0);

        ulong i = 0;
        for (; i+7 < len; i += 8)
        {
            umul_ppmm(s1, s0, v1[i+0], v2[i+0]);
            add_ssaaaa(u1, u0, u1, u0, s1, s0);
            umul_ppmm(s1, s0, v1[i+1], v2[i+1]);
            add_ssaaaa(u1, u0, u1, u0, s1, s0);
            umul_ppmm(s1, s0, v1[i+2], v2[i+2]);
            add_ssaaaa(u1, u0, u1, u0, s1, s0);
            umul_ppmm(s1, s0, v1[i+3], v2[i+3]);
            add_ssaaaa(u1, u0, u1, u0, s1, s0);
            umul_ppmm(s1, s0, v1[i+4], v2[i+4]);
            add_ssaaaa(u1, u0, u1, u0, s1, s0);
            umul_ppmm(s1, s0, v1[i+5], v2[i+5]);
            add_ssaaaa(u1, u0, u1, u0, s1, s0);
            umul_ppmm(s1, s0, v1[i+6], v2[i+6]);
            add_ssaaaa(u1, u0, u1, u0, s1, s0);
            umul_ppmm(s1, s0, v1[i+7], v2[i+7]);
            add_ssaaaa(u1, u0, u1, u0, s1, s0);
        }

        for (; i < len; i++)
        {
            umul_ppmm(s1, s0, v1[i], v2[i]);
            add_ssaaaa(u1, u0, u1, u0, s1, s0);
        }

        ulong res;
        NMOD2_RED2(res, u1, u0, mod);
        return res;
    }

    if (n_limbs == 3)
    {
        /* number of products we can do before overflow */
        const ulong maxbits = FLINT_BIT_COUNT(mod.n);
        const ulong log_nt = 2*FLINT_BITS - 2*maxbits;
        const ulong num_terms = (log_nt < FLINT_BITS) ? (UWORD(1) << log_nt) : (UWORD_MAX);

        ulong s0, s1, u0, u1;
        ulong t2 = UWORD(0);
        ulong t1 = UWORD(0);
        ulong t0 = UWORD(0);

        ulong i = 0;
        if (num_terms >= 8)
            for (; i+7 < len; i += 8) // FIXME 8 vs 4 ? (Vincent: not tested)
            {
                umul_ppmm(u1, u0, v1[i+0], v2[i+0]);
                umul_ppmm(s1, s0, v1[i+1], v2[i+1]);
                add_ssaaaa(u1, u0, u1, u0, s1, s0);
                umul_ppmm(s1, s0, v1[i+2], v2[i+2]);
                add_ssaaaa(u1, u0, u1, u0, s1, s0);
                umul_ppmm(s1, s0, v1[i+3], v2[i+3]);
                add_ssaaaa(u1, u0, u1, u0, s1, s0);
                umul_ppmm(s1, s0, v1[i+4], v2[i+4]);
                add_ssaaaa(u1, u0, u1, u0, s1, s0);
                umul_ppmm(s1, s0, v1[i+5], v2[i+5]);
                add_ssaaaa(u1, u0, u1, u0, s1, s0);
                umul_ppmm(s1, s0, v1[i+6], v2[i+6]);
                add_ssaaaa(u1, u0, u1, u0, s1, s0);
                umul_ppmm(s1, s0, v1[i+7], v2[i+7]);
                add_ssaaaa(u1, u0, u1, u0, s1, s0);
                add_sssaaaaaa(t2, t1, t0, t2, t1, t0, UWORD(0), u1, u0);
            }
        else
            for (; i+num_terms < len; i += num_terms)
            {
                umul_ppmm(u1, u0, v1[i], v2[i]);
                for (ulong ii = 1; ii < num_terms; ii++)
                {
                    umul_ppmm(s1, s0, v1[i+ii], v2[i+ii]);
                    add_ssaaaa(u1, u0, u1, u0, s1, s0);
                }
                add_sssaaaaaa(t2, t1, t0, t2, t1, t0, UWORD(0), u1, u0);
            }

        u0 = UWORD(0);
        u1 = UWORD(0);
        for (; i < len; i++)
        {
            umul_ppmm(s1, s0, v1[i], v2[i]);
            add_ssaaaa(u1, u0, u1, u0, s1, s0);
        }

        add_sssaaaaaa(t2, t1, t0, t2, t1, t0, UWORD(0), u1, u0);
        NMOD_RED(t2, t2, mod);

        ulong res;
        NMOD_RED3(res, t2, t1, t0, mod);
        return res;
    }

    return UWORD(0);
}



// NOTE: same constraints on len and mod as for avx version
ulong _nmod_vec_dot_mod32(nn_ptr a, nn_ptr b, ulong len, nmod_t mod, uint power_two)
{
    // the dot product without modular reduction is
    //       dp  =  dp_lo + 2**sp_nb * dp_hi
    // where sp_nb is a fixed number of bits where we split
    const ulong low_bits = (1L << DOT_SP_NB) - 1;
    ulong dp_lo0 = 0;
    ulong dp_lo1 = 0;
    ulong dp_lo2 = 0;
    ulong dp_lo3 = 0;
    ulong dp_hi0 = 0;
    ulong dp_hi1 = 0;
    ulong dp_hi2 = 0;
    ulong dp_hi3 = 0;

    // each pass in the loop handles 8*4 = 32 coefficients
    ulong k = 0;
    for (; k+31 < len; k += 32)
    {
        // no overflow in these lines if (2**sp_nb - 1) + 8*(p-1)**2 < 2**64    (*)
        dp_lo0 += a[k+ 0] * b[k+ 0];  dp_lo1 += a[k+ 1] * b[k+ 1];  dp_lo2 += a[k+ 2] * b[k+ 2];  dp_lo3 += a[k+ 3] * b[k+ 3];
        dp_lo0 += a[k+ 4] * b[k+ 4];  dp_lo1 += a[k+ 5] * b[k+ 5];  dp_lo2 += a[k+ 6] * b[k+ 6];  dp_lo3 += a[k+ 7] * b[k+ 7];
        dp_lo0 += a[k+ 8] * b[k+ 8];  dp_lo1 += a[k+ 9] * b[k+ 9];  dp_lo2 += a[k+10] * b[k+10];  dp_lo3 += a[k+11] * b[k+11];
        dp_lo0 += a[k+12] * b[k+12];  dp_lo1 += a[k+13] * b[k+13];  dp_lo2 += a[k+14] * b[k+14];  dp_lo3 += a[k+15] * b[k+15];
        dp_lo0 += a[k+16] * b[k+16];  dp_lo1 += a[k+17] * b[k+17];  dp_lo2 += a[k+18] * b[k+18];  dp_lo3 += a[k+19] * b[k+19];
        dp_lo0 += a[k+20] * b[k+20];  dp_lo1 += a[k+21] * b[k+21];  dp_lo2 += a[k+22] * b[k+22];  dp_lo3 += a[k+23] * b[k+23];
        dp_lo0 += a[k+24] * b[k+24];  dp_lo1 += a[k+25] * b[k+25];  dp_lo2 += a[k+26] * b[k+26];  dp_lo3 += a[k+27] * b[k+27];
        dp_lo0 += a[k+28] * b[k+28];  dp_lo1 += a[k+29] * b[k+29];  dp_lo2 += a[k+30] * b[k+30];  dp_lo3 += a[k+31] * b[k+31];

        // add bits sp_nb...63 to dp_hi; keep sp_nb lower bits in dp_lo
        dp_hi0 += (dp_lo0 >> DOT_SP_NB); dp_hi1 += (dp_lo1 >> DOT_SP_NB); dp_hi2 += (dp_lo2 >> DOT_SP_NB); dp_hi3 += (dp_lo3 >> DOT_SP_NB);
        dp_lo0 &= low_bits; dp_lo0 &= low_bits; dp_lo0 &= low_bits; dp_lo0 &= low_bits;
    }

    // left with len-k < 32 coefficients,
    // this handles (len-k)//4 < 8 blocks of 4 coefficients, no overflow in this loop if (*) satisfied
    for (; k + 3 < len; k += 4)
    {
        dp_lo0 += a[k+ 0] * b[k+ 0];  dp_lo1 += a[k+ 1] * b[k+ 1];  dp_lo2 += a[k+ 2] * b[k+ 2];  dp_lo3 += a[k+ 3] * b[k+ 3];
    }

    dp_hi0 = dp_hi0 + dp_hi1 + dp_hi2 + dp_hi3 + (dp_lo0 >> DOT_SP_NB) + (dp_lo1 >> DOT_SP_NB) + (dp_lo2 >> DOT_SP_NB) + (dp_lo3 >> DOT_SP_NB);
    dp_lo0 = (dp_lo0 & low_bits) + (dp_lo1 & low_bits) + (dp_lo2 & low_bits) + (dp_lo3 & low_bits);

    // left with < 4 coefficients
    // no overflow in dp_last if 3*(p-1)**2 < 2**64, ok if (*) satisfied
    ulong dp_last = 0;
    for (; k < len; k++)
        dp_last += a[k] * b[k];

    dp_lo0 += dp_last & low_bits;
    dp_hi0 += dp_last >> DOT_SP_NB;
    ulong dp = ((ulong)power_two * dp_hi0) + dp_lo0;
    NMOD_RED(dp, dp, mod);
    return dp;
}

// NOTE: same constraints on len and mod as for avx version
ulong _nmod_vec_dot_mod32_v2(nn_ptr a, nn_ptr b, ulong len, nmod_t mod, uint power_two)
{
    // the dot product without modular reduction is
    //       dp  =  dp_lo + 2**sp_nb * dp_hi
    // where sp_nb is a fixed number of bits where we split
    const ulong low_bits = (1L << DOT_SP_NB) - 1;
    ulong dp_lo = 0;
    ulong dp_hi = 0;

    // each pass in the loop handles 8 coefficients
    ulong k = 0;
    for (; k < len; k += 8)
    {
        // no overflow in these lines if (2**sp_nb - 1) + 8*(p-1)**2 < 2**64    (*)
        dp_lo += a[k+0] * b[k+0];
        dp_lo += a[k+1] * b[k+1];
        dp_lo += a[k+2] * b[k+2];
        dp_lo += a[k+3] * b[k+3];
        dp_lo += a[k+4] * b[k+4];
        dp_lo += a[k+5] * b[k+5];
        dp_lo += a[k+6] * b[k+6];
        dp_lo += a[k+7] * b[k+7];

        // add bits sp_nb...63 to dp_hi; keep sp_nb lower bits in dp_lo
        dp_hi += (dp_lo >> DOT_SP_NB);
        dp_lo &= low_bits;
    }

    // left with len-k < 8 coefficients
    // no overflow in dp_last if 7*(p-1)**2 < 2**64, ok FIXME if (*) satisfied
    ulong dp_last = 0;
    for (; k < len; k++)
        dp_last += a[k] * b[k];

    dp_lo += dp_last & low_bits;
    dp_hi += dp_last >> DOT_SP_NB;
    ulong dp = ((ulong)power_two * dp_hi) + dp_lo;
    NMOD_RED(dp, dp, mod);
    return dp;
}

/*-------------------------------*/
/* dot product for small modulus */
/*-------------------------------*/

// in short: with current DOT_SP_NB value 55,
// -> modulus up to about 2**30.5 (more precisely, up to 1517016615)
// -> length of dot product up to 


// NOTES:
// -> constraint 0:
// as can be seen below (*), we will have the condition
//      8 * (p-1)**2 <= 2**64 - 2**sp_nb
// i.e. (p-1)**2  <=  2**61 - 2**(sp_nb-3),
// one can take p up to 1 + floor(sqrt(2**61 - 2**(sp_nb-3)))
// in particular, p-1 < 2**30.5
//
// -> constraint 1:
// the above representation of dp requires len * (p-1)**2 <= 2**sp_nb * (2**64-1)
// since (p-1)**2 is < 2**61, a sufficient condition is len * 2**61 <= 2**sp_nb * 2**63,
// i.e. len <= 2**(sp_nb+2)   (which is not very restrictive)
//
// -> constraint 2:
// having dp_lo and dp_hi, we will actually compute dp_lo + power_two * dp_hi
// if computing power_two * dp_hi with single-word, this requires
//    (p-1) dp_hi < 2**64, 
// since (p-1) dp_hi <= (p-1) * len * (p-1)**2 / 2**sp_nb,
// it suffices to ensure
//     len * (p-1)**3 < 2**(64+sp_nb) 
// which is more restrictive than constraint 1 (as soon as p > 2; if p==2 both constraints are fine)
// since (p-1)**3 < 2**91.5, this accepts at least len <= 2**(sp_nb - 27.5)
//
//
// sage: for sp_nb in range(40,64):
// ....:     pmax = 1 + floor(sqrt(2**61 - 2**(sp_nb-3)))
// ....:     lenmax = ceil(2**(64+sp_nb) / (pmax-1)**3) - 1
// ....:     print(f"{sp_nb}\t{pmax.nbits()}\t{pmax}\t{lenmax.nbits()}\t{lenmax}")
// ....:
// 40      31      1518500205      13      5792
// 41      31      1518500160      14      11585
// 42      31      1518500069      15      23170
// 43      31      1518499888      16      46340
// 44      31      1518499526      17      92682
// 45      31      1518498802      18      185364
// 46      31      1518497354      19      370729
// 47      31      1518494458      20      741463
// 48      31      1518488665      21      1482944
// 49      31      1518477080      22      2965956
// 50      31      1518453909      23      5932184
// 51      31      1518407566      24      11865455
// 52      31      1518314875      25      23735258
// 53      31      1518129478      26      47487909
// 54      31      1517758614      27      95045458
// 55      31      1517016615      28      190369983
// 56      31      1515531528      29      381860339
// 57      31      1512556978      30      768235276
// 58      31      1506590261      31      1554798112
// 59      31      1494585366      32      3185130939
// 60      31      1470281545      33      6691414686
// 61      31      1420426920      34      14842011900
// 62      31      1315059793      36      37406145210
// 63      31      1073741825      37      137438953471

//
// NOTE: same constraints on len and mod as for avx version
ulong _nmod_vec_dot_mod32_v3(nn_ptr a, nn_ptr b, ulong len, nmod_t mod, uint power_two)
{
    // we compute dp_lo and dp_hi such that the dot product
    // without modular reduction is   dp  =  dp_lo + 2**DOT_SP_NB * dp_hi
    const ulong low_bits = (1L << DOT_SP_NB) - 1;
    ulong dp_lo = 0;
    uint dp_hi = 0;

    // handle batches of 8 coefficients
    // (experiments show benefit against 4; no real benefit for > 8)
    ulong k = 0;
    for (; k+7 < len; k+=8)
    {
        dp_lo += a[k+0] * b[k+0] + 
                 a[k+1] * b[k+1] + 
                 a[k+2] * b[k+2] + 
                 a[k+3] * b[k+3] +
                 a[k+4] * b[k+4] + 
                 a[k+5] * b[k+5] + 
                 a[k+6] * b[k+6] +
                 a[k+7] * b[k+7];

        // add bits sp_nb...63 to dp_hi; keep sp_nb lower bits in dp_lo
        dp_hi += dp_lo >> DOT_SP_NB;
        dp_lo &= low_bits;
    }

    // add the remaining < 8 coefficients
    for (; k < len; k++)
        dp_lo += a[k] * b[k];

    ulong dp;
    NMOD_RED(dp, ((ulong)power_two * dp_hi) + dp_lo, mod);
    return dp;
}



ulong nmod_vec_dot_mod32(nn_ptr a, nn_ptr b, ulong len, nmod_t mod)
{
    // FIXME add check for len and p?
    uint power_two;
    NMOD_RED(power_two, 1L<<DOT_SP_NB, mod);
    return _nmod_vec_dot_mod32(a, b, len, mod, power_two);
}



/*------------------------------------------------------------*/
/* AVX2 for small moduli (up to ~31 bits)                     */
/*------------------------------------------------------------*/


ulong _nmod_vec_dot_mod32_avx2(nn_ptr a, nn_ptr b, ulong len, nmod_t mod, uint power_two)
{
    // the dot product without modular reduction is
    //       dp  =  dp_lo + 2**sp_nb * dp_hi
    // where sp_nb is a fixed number of bits where we split
    const vec4n low_bits = vec4n_set_n((1L << DOT_SP_NB) - 1);
    vec4n dp_lo = vec4n_zero();
    vec4n dp_hi = vec4n_zero();

    // each pass in the loop handles 8x4 = 32 coefficients
    ulong k = 0;
    for (; k+31 < len; k += 32)
    {
        // no overflow in next 8 lines if (2**sp_nb - 1) + 8*(p-1)**2 < 2**64    (*)
        dp_lo = vec4n_add(dp_lo, vec4n_mul(vec4n_load_unaligned(a+k+ 0), vec4n_load_unaligned(b+k+ 0)));
        dp_lo = vec4n_add(dp_lo, vec4n_mul(vec4n_load_unaligned(a+k+ 4), vec4n_load_unaligned(b+k+ 4)));
        dp_lo = vec4n_add(dp_lo, vec4n_mul(vec4n_load_unaligned(a+k+ 8), vec4n_load_unaligned(b+k+ 8)));
        dp_lo = vec4n_add(dp_lo, vec4n_mul(vec4n_load_unaligned(a+k+12), vec4n_load_unaligned(b+k+12)));
        dp_lo = vec4n_add(dp_lo, vec4n_mul(vec4n_load_unaligned(a+k+16), vec4n_load_unaligned(b+k+16)));
        dp_lo = vec4n_add(dp_lo, vec4n_mul(vec4n_load_unaligned(a+k+20), vec4n_load_unaligned(b+k+20)));
        dp_lo = vec4n_add(dp_lo, vec4n_mul(vec4n_load_unaligned(a+k+24), vec4n_load_unaligned(b+k+24)));
        dp_lo = vec4n_add(dp_lo, vec4n_mul(vec4n_load_unaligned(a+k+28), vec4n_load_unaligned(b+k+28)));

        // add bits sp_nb...63 to dp_hi; keep sp_nb lower bits in dp_lo
        // FIXME vec4n_bit_shift_right  uses _mm256_srl_epi64, a priori slower?
        dp_hi = vec4n_add(dp_hi, _mm256_srli_epi64(dp_lo, DOT_SP_NB));
        dp_lo = vec4n_bit_and(dp_lo, low_bits);
    }

    // left with len-k < 32 coefficients,
    // this handles (len-k)//4 < 8 blocks of 4 coefficients, no overflow in this loop if (*) satisfied
    for (; k + 3 < len; k += 4)
        dp_lo = vec4n_add(dp_lo, vec4n_mul(vec4n_load_unaligned(a+k), vec4n_load_unaligned(b+k)));

    dp_hi = vec4n_add(dp_hi, _mm256_srli_epi64(dp_lo, DOT_SP_NB));
    dp_lo = vec4n_bit_and(dp_lo, low_bits);

    // left with < 4 coefficients
    // no overflow in dp_last if 3*(p-1)**2 < 2**64, ok if (*) satisfied
    ulong dp_last = 0;
    for (; k < len; k++)
        dp_last += a[k] * b[k];

    const ulong total_lo = dp_lo[0] + dp_lo[1] + dp_lo[2] + dp_lo[3] + (dp_last & ((1L << DOT_SP_NB) - 1));
    const uint total_hi = dp_hi[0] + dp_hi[1] + dp_hi[2] + dp_hi[3] + (dp_last >> DOT_SP_NB);
    // FIXME check uint is fine; may add additional requirement
    ulong dp = ((ulong)power_two * total_hi) + total_lo;
    NMOD_RED(dp, dp, mod);
    return dp;
}

ulong nmod_vec_dot_mod32_avx2(nn_ptr a, nn_ptr b, ulong len, nmod_t mod)
{
    // FIXME add check for len and p?
    uint power_two;
    NMOD_RED(power_two, 1L<<DOT_SP_NB, mod);
    return _nmod_vec_dot_mod32_avx2(a, b, len, mod, power_two);
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
