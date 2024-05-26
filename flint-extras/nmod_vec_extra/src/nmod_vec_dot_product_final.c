#include <flint/nmod.h>

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

// dot product using single limb, avx2
ulong _nmod_vec_dot_product_1_avx2(nn_srcptr vec1, nn_srcptr vec2, ulong len, nmod_t mod)
{
    // compute 4 vertical sub-dot products
    __m256i res_vec = _mm256_setzero_si256();
    ulong i;
    for (i = 0; i+3 < len; i += 4)
    {
        // load + multiplication + addition
        __m256i x = _mm256_loadu_si256((__m256i *) (vec1+i));
        __m256i y = _mm256_loadu_si256((__m256i *) (vec2+i));
        x = _mm256_mul_epu32(x, y);
        res_vec = _mm256_add_epi64(res_vec, x);
    }

    // horizontal add
    ulong res = res_vec[0] + res_vec[1] + res_vec[2] + res_vec[3];

    // scalar loop for leftover entries
    for (; i < len; ++i)
        res += vec1[i] * vec2[i];
    NMOD_RED(res, res, mod);

    return res;
}

// dot product using single limb, avx512
ulong _nmod_vec_dot_product_1_avx512(nn_srcptr vec1, nn_srcptr vec2, ulong len, nmod_t mod) {
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

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
