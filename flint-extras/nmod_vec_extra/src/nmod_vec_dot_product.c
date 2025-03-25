#include <flint/longlong.h>
#include <flint/nmod.h>
#include <flint/nmod_vec.h>
#include <immintrin.h>

#include "nmod_vec_extra.h"

/************
*  hsum
*  https://stackoverflow.com/questions/60108658/fastest-method-to-calculate-sum-of-all-packed-32-bit-integers-using-avx512-or-av
************/

//static inline
//uint hsum_epi32_avx(__m128i x)
//{
//    __m128i hi64  = _mm_unpackhi_epi64(x, x);           // 3-operand non-destructive AVX lets us save a byte without needing a movdqa
//    __m128i sum64 = _mm_add_epi32(hi64, x);
//    __m128i hi32  = _mm_shuffle_epi32(sum64, _MM_SHUFFLE(2, 3, 0, 1));    // Swap the low two elements
//    __m128i sum32 = _mm_add_epi32(sum64, hi32);
//    return _mm_cvtsi128_si32(sum32);       // movd
//}
//
//// only needs AVX2
//uint hsum_8x32(__m256i v)
//{
//    __m128i sum128 = _mm_add_epi32(
//                 _mm256_castsi256_si128(v),
//                 _mm256_extracti128_si256(v, 1)); // silly GCC uses a longer AXV512VL instruction if AVX512 is enabled :/
//    return hsum_epi32_avx(sum128);
//}



/* ------------------------------------------------------------ */
/* number of limbs needed for a dot product of length len       */
/* all entries 1st vector have <= max_bits1 bits <= FLINT_BITS  */
/* all entries 2nd vector have <= max_bits2 bits <= FLINT_BITS  */
/* returns 0, 1, 2, 3                                           */
/* ------------------------------------------------------------ */
// TODO adapt from new flint dot_params
static inline
ulong _nmod_vec_dot_bound_limbs_unbalanced(ulong len, ulong max_bits1, ulong max_bits2)
{
    const ulong a1 = (max_bits1 == FLINT_BITS) ? (UWORD_MAX) : (UWORD(1) << max_bits1) - 1;
    const ulong a2 = (max_bits2 == FLINT_BITS) ? (UWORD_MAX) : (UWORD(1) << max_bits2) - 1;

    ulong t2, t1, t0, u1, u0;
    umul_ppmm(t1, t0, a1, a2);
    umul_ppmm(t2, t1, t1, len);
    umul_ppmm(u1, u0, t0, len);
    add_sssaaaaaa(t2, t1, t0,  t2, t1, UWORD(0),  UWORD(0), u1, u0);

    if (t2 != 0)
        return 3;
    if (t1 != 0)
        return 2;
    return (t0 != 0);
}

/*  ------------------------------------------------------------ */
/** v1 and v2 have length at least len, len < 2^FLINT_BITS       */
/** computes sum(v1[i]*v2[i], 0 <= i < len) modulo mod.n         */
/** uses 1 limb to store the result before reduction             */
/*  ------------------------------------------------------------ */
static inline
ulong _nmod_vec_dot_product_1(nn_srcptr v1, nn_srcptr v2, ulong len, nmod_t mod)
{
    ulong res; slong i;
    _NMOD_VEC_DOT1(res, i, len, v1[i], v2[i], mod);
    return res;
}

/*  ------------------------------------------------------------ */
/** v1 and v2 have length at least len, len < 2^FLINT_BITS       */
/** computes sum(v1[i]*v2[i], 0 <= i < len) modulo mod.n         */
/** uses 2 limbs to store the result before reduction            */
/*  ------------------------------------------------------------ */
static inline
ulong _nmod_vec_dot_product_2(nn_srcptr v1, nn_srcptr v2, ulong len, nmod_t mod)
{
    ulong res; slong i;
    _NMOD_VEC_DOT2(res, i, len, v1[i], v2[i], mod);
    return res;
}

/*  ------------------------------------------------------------ */
/** v1 and v2 have length at least len, len < 2^FLINT_BITS       */
/** all entries of v1 have <= max_bits1 bits <= FLINT_BITS       */
/** all entries of v2 have <= max_bits2 bits <= FLINT_BITS       */
/** computes sum(v1[i]*v2[i], 0 <= i < len) modulo mod.n         */
/** uses 3 limbs to store the result before reduction            */
/*  ------------------------------------------------------------ */
// TODO adapt with new flint dot_params
static inline
ulong _nmod_vec_dot_product_3(nn_srcptr v1, nn_srcptr v2, ulong len, ulong max_bits1, ulong max_bits2, nmod_t mod)
{
    /* number of products we can do before overflow */
    const ulong log_nt = 2*FLINT_BITS - (max_bits1 + max_bits2);
    const ulong num_terms = (log_nt < FLINT_BITS) ? (UWORD(1) << log_nt) : (UWORD_MAX);

    ulong s0, s1, u0, u1;
    ulong t2 = UWORD(0);
    ulong t1 = UWORD(0);
    ulong t0 = UWORD(0);

    ulong i = 0;
    if (num_terms >= 8)
        for (; i+7 < len; i += 8)
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

/*  ------------------------------------------------------------ */
/** v1 and v2 have length at least len, len <= 2^FLINT_BITS      */
/** all entries of v1 have <= max_bits1 bits <= FLINT_BITS       */
/** all entries of v2 have <= max_bits2 bits <= FLINT_BITS       */
/** computes sum(v1[i]*v2[i], 0 <= i < len) modulo mod.n         */
/** does not assume input is reduced modulo mod.n                */
/*  ------------------------------------------------------------ */
ulong nmod_vec_dot_product_unbalanced(nn_srcptr v1, nn_srcptr v2, ulong len, ulong max_bits1, ulong max_bits2, nmod_t mod)
{
    const ulong n_limbs = _nmod_vec_dot_bound_limbs_unbalanced(len, max_bits1, max_bits2);

    if (n_limbs == 2)
        return _nmod_vec_dot_product_2(v1, v2, len, mod);
    if (n_limbs == 3)
        return _nmod_vec_dot_product_3(v1, v2, len, max_bits1, max_bits2, mod);
    if (n_limbs == 1)
        return _nmod_vec_dot_product_1(v1, v2, len, mod);
    return UWORD(0);
}








/*------------------------------------------------------------*/
/* EXPERIMENTAL */
/*------------------------------------------------------------*/

#if HAVE_AVX512
// dot product using single limb, avx512
ulong _nmod_vec_dot_product_1_avx512(nn_srcptr vec1, nn_srcptr vec2, ulong len, nmod_t mod)
{
    // compute 4 vertical sub-dot products
    __m512i res_vec = _mm512_setzero_si512();
    ulong i;
    for (i = 0; i+7 < len; i += 8)
    {
        // load + multiplication + addition
        __m512i x = _mm512_loadu_si512((__m512i *) (vec1+i));
        __m512i y = _mm512_loadu_si512((__m512i *) (vec2+i));
        x = _mm512_mul_epu32(x, y);
        res_vec = _mm512_add_epi64(res_vec, x);
    }

    // horizontal add
    //ulong res = res_vec[0] + res_vec[1] + res_vec[2] + res_vec[3]
              //+ res_vec[4] + res_vec[5] + res_vec[6] + res_vec[7];
    ulong res = _mm512_reduce_add_epi64(res_vec);

    // scalar loop for leftover entries
    for (; i < len; ++i)
        res += vec1[i] * vec2[i];
    NMOD_RED(res, res, mod);

    return res;
}
#endif


#define DOT_SPLITXX_BITS 26
#define DOT_SPLITXX_MASK 0x3FFFFFF // (1L << DOT_SPLITXX_BITS) - 1

#define __ll_lowhi_partsXX(tlo,thi,t)     \
      thi = (unsigned int) ((t) >> DOT_SPLITXX_BITS);           \
      tlo = ((unsigned int)(t)) & DOT_SPLITXX_MASK;
#include <flint/machine_vectors.h>

ulong _nmod_vec_dot2_half_avx(nn_srcptr v1, nn_srcptr v2, ulong len, nmod_t mod)
{
    const vec4n low_bits = vec4n_set_n(DOT_SPLIT_MASK);
    vec4n dp_lo = vec4n_zero();
    vec4n dp_hi = vec4n_zero();

    ulong i = 0;
    // DOT_SPLIT_BITS == 56: we can accumulate up to 2**8 == 256 integers of <= DOT_SPLIT_BITS bits without overflow
    for ( ; i+255 < len; i += 256)
    {
        ulong j = 0;
        for ( ; j+3 < 256; j += 4)
        {
            __m256i prod = vec4n_mul(vec4n_load_unaligned(v1+i+j), vec4n_load_unaligned(v2+i+j));
            dp_hi = vec4n_add(dp_hi, vec4n_bit_shift_right(prod, DOT_SPLIT_BITS));
            dp_lo = vec4n_add(dp_lo, vec4n_bit_and(prod, low_bits));
        }
        // dp_lo might be very close to full 64 bits: move its bits 56..63 to dp_hi
        dp_hi = vec4n_add(dp_hi, vec4n_bit_shift_right(dp_lo, DOT_SPLIT_BITS));
        dp_lo = vec4n_bit_and(dp_lo, low_bits);
    }

    // less than 256 terms remaining
    // we can accumulate all of the next <= 252 ones
    for ( ; i+3 < len; i += 4)
    {
        __m256i prod = vec4n_mul(vec4n_load_unaligned(v1+i), vec4n_load_unaligned(v2+i));
        dp_hi = vec4n_add(dp_hi, vec4n_bit_shift_right(prod, DOT_SPLIT_BITS));
        dp_lo = vec4n_add(dp_lo, vec4n_bit_and(prod, low_bits));
    }

    // since only <= 252 were accumulated, we can safely sum 4 terms horizontally
    ulong hsum_lo = vec4n_horizontal_sum(dp_lo);
    ulong hsum_hi = vec4n_horizontal_sum(dp_hi) + (hsum_lo >> DOT_SPLIT_BITS);
    hsum_lo &= DOT_SPLIT_MASK;

    for ( ; i < len; i++)
    {
        ulong prod = v1[i] * v2[i];
        hsum_hi += (prod >> DOT_SPLIT_BITS);
        hsum_lo += (prod & DOT_SPLIT_MASK);
    }

    ulong res;
    // TODO replace this with some powmod2_precomp
    NMOD_RED(res, ((1L<<DOT_SPLIT_BITS) % mod.n) * hsum_hi + hsum_lo, mod);
    return res;
}

ulong _nmod_vec_dot2_half_avx_int(const uint * v1, const uint * v2, ulong len, nmod_t mod)
{
    const vec4n low_bits = vec4n_set_n(DOT_SPLIT_MASK);
    vec4n dp_lo = vec4n_zero();
    vec4n dp_hi = vec4n_zero();

    ulong i = 0;
    // DOT_SPLIT_BITS == 56: we can accumulate up to 2**8 == 256 integers of <= DOT_SPLIT_BITS bits without overflow
    for ( ; i+255 < len; i += 256)
    {
        ulong j = 0;
        for ( ; j+7 < 256; j += 8)
        {
            __m256i v1i = _mm256_loadu_si256((__m256i *) (v1+i+j));
            __m256i v2i = _mm256_loadu_si256((__m256i *) (v2+i+j));
            __m256i prod = vec4n_mul(v1i, v2i);
            dp_hi = vec4n_add(dp_hi, vec4n_bit_shift_right(prod, DOT_SPLIT_BITS));
            dp_lo = vec4n_add(dp_lo, vec4n_bit_and(prod, low_bits));
            prod = vec4n_mul(_mm256_srli_epi64(v1i, 32), _mm256_srli_epi64(v2i, 32));
            dp_hi = vec4n_add(dp_hi, vec4n_bit_shift_right(prod, DOT_SPLIT_BITS));
            dp_lo = vec4n_add(dp_lo, vec4n_bit_and(prod, low_bits));
        }
        // dp_lo might be very close to full 64 bits: move its bits 56..63 to dp_hi
        dp_hi = vec4n_add(dp_hi, vec4n_bit_shift_right(dp_lo, DOT_SPLIT_BITS));
        dp_lo = vec4n_bit_and(dp_lo, low_bits);
    }

    // less than 256 terms remaining
    // we can accumulate all of the next <= 248 ones
    for ( ; i+7 < len; i += 8)
    {
        __m256i v1i = _mm256_loadu_si256((__m256i *) (v1+i));
        __m256i v2i = _mm256_loadu_si256((__m256i *) (v2+i));
        __m256i prod = vec4n_mul(v1i, v2i);
        dp_hi = vec4n_add(dp_hi, vec4n_bit_shift_right(prod, DOT_SPLIT_BITS));
        dp_lo = vec4n_add(dp_lo, vec4n_bit_and(prod, low_bits));
        prod = vec4n_mul(_mm256_srli_epi64(v1i, 32), _mm256_srli_epi64(v2i, 32));
        dp_hi = vec4n_add(dp_hi, vec4n_bit_shift_right(prod, DOT_SPLIT_BITS));
        dp_lo = vec4n_add(dp_lo, vec4n_bit_and(prod, low_bits));
    }

    // since only <= 248 were accumulated, we can safely sum 4 terms horizontally
    ulong hsum_lo = vec4n_horizontal_sum(dp_lo);
    ulong hsum_hi = vec4n_horizontal_sum(dp_hi) + (hsum_lo >> DOT_SPLIT_BITS);
    hsum_lo &= DOT_SPLIT_MASK;

    for ( ; i < len; i++)
    {
        ulong prod = (ulong)v1[i] * v2[i];
        hsum_hi += (prod >> DOT_SPLIT_BITS);
        hsum_lo += (prod & DOT_SPLIT_MASK);
    }

    ulong res;
    // TODO replace this with some powmod2_precomp
    NMOD_RED(res, ((1L<<DOT_SPLIT_BITS) % mod.n) * hsum_hi + hsum_lo, mod);
    return res;
}

// TODO benchmark more, integrate, give precise conditions for when this works
// (or better, really do a hand-made avx512 version...)
// --> if splitting at 26, each product is 52, can allow at most 12 additional bits,
// i.e. not more than xxx terms (this depends on the size of the high part since
// they are not balanced... could make sense to balance them to allow more terms,
// but do this only if this really is interesting in terms of speed)
ulong _nmod_vec_dot_product_split26(nn_srcptr v1, nn_srcptr v2, ulong len, nmod_t mod)
{
    uint v1hi, v1lo, v2hi, v2lo;
    ulong ulo = UWORD(0);
    ulong umi = UWORD(0);
    ulong uhi = UWORD(0);
    for (ulong i = 0; i < len; i++)
    {
        __ll_lowhi_partsXX(v1lo, v1hi, v1[i]);
        __ll_lowhi_partsXX(v2lo, v2hi, v2[i]);
        ulo += (ulong)v1lo * v2lo;
        umi += (ulong)v1lo * v2hi + (ulong)v1hi * v2lo;
        uhi += (ulong)v1hi * v2hi;
    }

    // result: ulo + 2**26 umi + 2**52 uhi
    // hi = (umi >> 38) + (uhi >> 12)  ||  lo = (umi << 26) + (uhi << 52) + ulo
    add_ssaaaa(uhi, ulo, (umi>>38), (umi<<26), (uhi>>12), ((uhi<<52)+ulo));
    ulong res;
    NMOD2_RED2(res, uhi, ulo, mod);
    return res;
}

ulong _nmod_vec_dot_product_split26_avx(nn_srcptr v1, nn_srcptr v2, ulong len, nmod_t mod)
{
    const vec4n low_bits = vec4n_set_n(DOT_SPLITXX_MASK);
    ulong i = 0;
    vec4n dp_vlo = vec4n_zero();
    vec4n dp_vmi = vec4n_zero();
    vec4n dp_vhi = vec4n_zero();
    for ( ; i+3 < len; i+=4)
    {
        vec4n v1hi = vec4n_load_unaligned(v1+i);
        vec4n v1lo = vec4n_bit_and(v1hi, low_bits);
        v1hi = vec4n_bit_shift_right(v1hi, DOT_SPLITXX_BITS);
        vec4n v2hi = vec4n_load_unaligned(v2+i);
        vec4n v2lo = vec4n_bit_and(v2hi, low_bits);
        v2hi = vec4n_bit_shift_right(v2hi, DOT_SPLITXX_BITS);

        dp_vlo = vec4n_add(dp_vlo, vec4n_mul(v1lo, v2lo));
        dp_vmi = vec4n_add(dp_vmi, vec4n_mul(v1lo, v2hi));
        dp_vmi = vec4n_add(dp_vmi, vec4n_mul(v1hi, v2lo));
        dp_vhi = vec4n_add(dp_vhi, vec4n_mul(v1hi, v2hi));
    }

    ulong dp_lo = vec4n_horizontal_sum(dp_vlo);
    ulong dp_mi = vec4n_horizontal_sum(dp_vmi);
    ulong dp_hi = vec4n_horizontal_sum(dp_vhi);

    for ( ; i < len; i++)
    {
        unsigned int v1hi, v1lo, v2hi, v2lo;
        __ll_lowhi_partsXX(v1lo, v1hi, v1[i]);
        __ll_lowhi_partsXX(v2lo, v2hi, v2[i]);
        dp_lo += (ulong)v1lo * v2lo;
        dp_mi += (ulong)v1lo * v2hi + (ulong)v1hi * v2lo;
        dp_hi += (ulong)v1hi * v2hi;
    }

    // result: ulo + 2**26 umi + 2**52 uhi
    // hi = (umi >> 38) + (uhi >> 12)  ||  lo = (umi << 26) + (uhi << 52) + ulo
    add_ssaaaa(dp_hi, dp_lo, (dp_mi>>38), (dp_mi<<26), (dp_hi>>12), ((dp_hi<<52)+dp_lo));
    ulong res;
    NMOD2_RED2(res, dp_hi, dp_lo, mod);
    return res;
}

#if HAVE_AVX_IFMA
ulong _nmod_vec_dot_product_avx_ifma(nn_srcptr v1, nn_srcptr v2, ulong len, nmod_t mod)
{
    ulong i = 0;
    vec4n dp_vlo = vec4n_zero();
    vec4n dp_vhi = vec4n_zero();
    for ( ; i+3 < len; i+=4)
    {
        __m256i v1i = vec4n_load_unaligned(v1+i);
        __m256i v2i = vec4n_load_unaligned(v2+i);
        dp_vlo = _mm256_madd52lo_epu64(dp_vlo, v1i, v2i);
        dp_vhi = _mm256_madd52hi_epu64(dp_vhi, v1i, v2i);
    }

    ulong dp_lo = vec4n_horizontal_sum(dp_vlo);
    ulong dp_hi = vec4n_horizontal_sum(dp_vhi);

    add_ssaaaa(dp_hi, dp_lo, 0, (dp_hi<<52), (dp_hi>>12), dp_lo);
    for ( ; i < len; i++)
    {
        ulong s0, s1;
        umul_ppmm(s1, s0, v1[i], v2[i]);
        add_ssaaaa(dp_hi, dp_lo, dp_hi, dp_lo, s1, s0);
    }

    ulong res;
    NMOD2_RED2(res, dp_hi, dp_lo, mod);
    return res;
}
#endif

#if HAVE_AVX512
ulong _nmod_vec_dot_product_avx512_ifma(nn_srcptr v1, nn_srcptr v2, ulong len, nmod_t mod)
{
    ulong i = 0;
    __m512i dp_vlo = _mm512_setzero_si512();
    __m512i dp_vhi = _mm512_setzero_si512();
    for ( ; i+7 < len; i+=8)
    {
        __m512i v1i = _mm512_loadu_si512(v1+i);
        __m512i v2i = _mm512_loadu_si512(v2+i);
        dp_vlo = _mm512_madd52lo_epu64(dp_vlo, v1i, v2i);
        dp_vhi = _mm512_madd52hi_epu64(dp_vhi, v1i, v2i);
    }

    ulong dp_lo = _mm512_reduce_add_epi64(dp_vlo);
    ulong dp_hi = _mm512_reduce_add_epi64(dp_vhi);

    add_ssaaaa(dp_hi, dp_lo, 0, (dp_hi<<52), (dp_hi>>12), dp_lo);
    for ( ; i < len; i++)
    {
        ulong s0, s1;
        umul_ppmm(s1, s0, v1[i], v2[i]);
        add_ssaaaa(dp_hi, dp_lo, dp_hi, dp_lo, s1, s0);
    }

    ulong res;
    NMOD2_RED2(res, dp_hi, dp_lo, mod);
    return res;
}
#endif


/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
