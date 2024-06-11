#include <flint/nmod.h>

#include "nmod_extra.h"
#include "nmod_vec_extra.h"

/*-------------------------------*/
/* dot product for small modulus */
/*-------------------------------*/

// in short: with current DOT_SP_NB value 56,
// -> modulus n up to about 2**30.5
//       (more precisely, n <= 1515531528)
// -> length of dot product up to at least 134744072
//       (more precisely, len <= (2**88 - 2**56) / (n-1)**2)

// APPROACH:
//
// Let n = mod.n, s = DOT_SP_NB
// As input, take power_two == 2**s % n
//
// -> avoiding modular reductions altogether, compute dp_lo and dp_hi such that
// the dot product without modular reduction is dp  =  dp_lo + 2**s * dp_hi
// -> finally, compute (dp_lo + power_two * dp_hi)  %  n
// -> done through repeating this: accumulate a few terms,
// save higher bits in dp_hi and lower ones in dp_lo

// PARAMETER CONSTRAINTS:
//
// -> constraint (C0):
// we will accumulate 8 terms (each a product of two integers reduced modulo n)
// on top of an s-bit integer, so we require
//     2**s - 1 + 8 * (p-1)**2  <  2**64
// so one can take any modulus with
//     n <= 1 + floor(sqrt(2**61 - 2**(s-3)))
// in particular, n-1 < 2**30.5, (n-1)**2 < 2**61, (n-1)**3 < 2**91.5
//
// -> constraint (C1):
// in the above representation of dp we will use an uint32 for dp_hi,
// so we require      len * (n-1)**2 <= 2**s * (2**32 - 1)
//
// -> constraint (C2):
// [this constraint is void since dp_hi fits in uint32, but these notes are left
// in case one would want to use ulong for dp_hi, which allows for larger lengths]
// if ones wishes power_two * dp_hi to fit in a single word, this requires
//    (n-1) dp_hi < 2**64
// and since (n-1) dp_hi <= (n-1) * len * (n-1)**2 / 2**s, it suffices to ensure
//     len * (n-1)**3 < 2**(64+s) 
//
// sage: for s in range(40,64):
// ....:     nmax = 1 + floor(sqrt(2**61 - 2**(s-3)))             # (C0)
// ....:     lenmax = floor(2**s * (2**32 - 1) / (nmax-1)**2)     # (C1)
// ....:     lenmax_bis = ceil(2**(64+s) / (nmax-1)**3) - 1       # (C2)
// ....:     print(f"{s}\t{nmax.nbits()}\t{nmax}\t{lenmax}\t{lenmax_bis}")
// ....:
// s       nbits   nmax            (C1) for nmax     (C2) for nmax
// 40      31      1518500205      2048            5792
// 41      31      1518500160      4096            11585
// 42      31      1518500069      8192            23170
// 43      31      1518499888      16384           46340
// 44      31      1518499526      32768           92682
// 45      31      1518498802      65536           185364
// 46      31      1518497354      131072          370729
// 47      31      1518494458      262146          741463
// 48      31      1518488665      524296          1482944
// 49      31      1518477080      1048608         2965956
// 50      31      1518453909      2097280         5932184
// 51      31      1518407566      4194816         11865455
// 52      31      1518314875      8390656         23735258
// 53      31      1518129478      16785412        47487909
// 54      31      1517758614      33587232        95045458
// 55      31      1517016615      67240192        190369983
// 56      31      1515531528      134744072       381860339
// 57      31      1512556978      270549121       768235276
// 58      31      1506590261      545392672       1554798112
// 59      31      1494585366      1108378657      3185130939
// 60      31      1470281545      2290649226      6691414686
// 61      31      1420426920      4908534053      14842011900
// 62      31      1315059793      11453246122     37406145210
// 63      31      1073741825      34359738360     137438953471


// scalar version (may be automatically vectorized to some extent)
ulong _nmod_vec_dot_mod32(nn_ptr a, nn_ptr b, ulong len, nmod_t mod, uint power_two)
{
    // we compute dp_lo and dp_hi such that the dot product
    // without modular reduction is   dp  =  dp_lo + 2**DOT_SP_NB * dp_hi
    const ulong low_bits = (1L << DOT_SP_NB) - 1;
    ulong dp_lo = 0;
    uint dp_hi = 0;

    // handle batches of 8 terms
    // no overflow thanks to condition (C0)
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

    // handle remaining terms (< 8, so again no overflow thank to (C0))
    for (; k < len; k++)
        dp_lo += a[k] * b[k];

    // compute and reduce dp_lo + power_two * dp_hi
    // no overflow in the multiplication thanks to (C1)
    ulong dp;
    NMOD_RED(dp, ((ulong)power_two * dp_hi) + dp_lo, mod);
    return dp;
}


// explicitly vectorized version, flag AVX2
ulong _nmod_vec_dot_mod32_avx2(nn_ptr a, nn_ptr b, ulong len, nmod_t mod, uint power_two)
{
    // we compute 4n-vectors dp_lo and dp_hi such that the dot product
    // without modular reduction is   dp  =  sum(dp_lo) + 2**DOT_SP_NB * sum(dp_hi)
    const vec4n low_bits = vec4n_set_n((1L << DOT_SP_NB) - 1);
    vec4n dp_lo = vec4n_zero();
    vec4n dp_hi = vec4n_zero();

    // handle batches of 8 terms (hence 32 coefficients per loop iteration)
    // no overflow thanks to condition (C0)
    ulong k = 0;
    for (; k+31 < len; k += 32)
    {
        dp_lo = vec4n_add(dp_lo, vec4n_mul(vec4n_load_unaligned(a+k+ 0), vec4n_load_unaligned(b+k+ 0)));
        dp_lo = vec4n_add(dp_lo, vec4n_mul(vec4n_load_unaligned(a+k+ 4), vec4n_load_unaligned(b+k+ 4)));
        dp_lo = vec4n_add(dp_lo, vec4n_mul(vec4n_load_unaligned(a+k+ 8), vec4n_load_unaligned(b+k+ 8)));
        dp_lo = vec4n_add(dp_lo, vec4n_mul(vec4n_load_unaligned(a+k+12), vec4n_load_unaligned(b+k+12)));
        dp_lo = vec4n_add(dp_lo, vec4n_mul(vec4n_load_unaligned(a+k+16), vec4n_load_unaligned(b+k+16)));
        dp_lo = vec4n_add(dp_lo, vec4n_mul(vec4n_load_unaligned(a+k+20), vec4n_load_unaligned(b+k+20)));
        dp_lo = vec4n_add(dp_lo, vec4n_mul(vec4n_load_unaligned(a+k+24), vec4n_load_unaligned(b+k+24)));
        dp_lo = vec4n_add(dp_lo, vec4n_mul(vec4n_load_unaligned(a+k+28), vec4n_load_unaligned(b+k+28)));

        // add bits sp_nb...63 to dp_hi; keep sp_nb lower bits in dp_lo
        dp_hi = vec4n_add(dp_hi, vec4n_bit_shift_right(dp_lo, DOT_SP_NB));
        dp_lo = vec4n_bit_and(dp_lo, low_bits);
    }

    // handle remaining terms 4 by 4
    // (at most 7 x 4 remaining terms, so again no overflow thanks to (C0))
    for (; k + 3 < len; k += 4)
        dp_lo = vec4n_add(dp_lo, vec4n_mul(vec4n_load_unaligned(a+k), vec4n_load_unaligned(b+k)));

    dp_hi = vec4n_add(dp_hi, _mm256_srli_epi64(dp_lo, DOT_SP_NB));
    dp_lo = vec4n_bit_and(dp_lo, low_bits);

    ulong total_lo = dp_lo[0] + dp_lo[1] + dp_lo[2] + dp_lo[3];
    const uint total_hi = dp_hi[0] + dp_hi[1] + dp_hi[2] + dp_hi[3] + (total_lo >> DOT_SP_NB);
    total_lo &= (1L << DOT_SP_NB) - 1;

    // left with < 4 terms (no overflow thanks to (C0))
    for (; k < len; k++)
        total_lo += a[k] * b[k];

    // compute and reduce dp_lo + power_two * dp_hi
    // no overflow in the multiplication thanks to (C1)
    ulong dp;
    NMOD_RED(dp, ((ulong)power_two * total_hi) + total_lo, mod);
    return dp;
}


ulong nmod_vec_dot_product(nn_srcptr v1, nn_srcptr v2, ulong len, nmod_t mod, ulong n_limbs)
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




/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
