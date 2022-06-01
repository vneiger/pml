#include "nmod_poly_extra.h"

#ifdef HAS_AVX2
#include <immintrin.h>

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* using AVX2 to do FFT nmod_poly with p < 2^32               */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* in-place FFT in size 2^3                                   */
/* powers_w_in are the roots of 1 needed                      */
/* input/output is a vector of 32 bit unsigned ints           */
/*------------------------------------------------------------*/
static void _fft_3(mp_hlimb_t * x,
                   mp_hlimb_t * powers_w, mp_hlimb_t * i_powers_w,
                   const nmod_t mod)
{
    slong r, i;
    mp_hlimb_t p, p2, w, iw;
    mp_hlimb_t  *x0, *x1;
        
    p = (mp_hlimb_t) mod.n;
    p2 = p + p;

    x0 = x;
    x1 = x + 4;

    for (i = 0; i < 4; i++)
    {
        mp_hlimb_t u0, u1, t1;
        mp_hlimb_signed_t t0;
        u0 = x0[i];
        u1 = x1[i];
        t0 = u0 + u1 - p2;
        x0[i] = (t0 < 0) ? t0+p2 : t0;
        t1 = (u0 + p2) - u1;
        x1[i] = mul_mod_precon_unreduced_32(t1, powers_w[i], p, i_powers_w[i]);
    }
    x0 += 8;
    x1 += 8;

    powers_w += 4;
    i_powers_w += 4;
    
    // last two layers
    w = powers_w[1];
    iw = i_powers_w[1];

    for (r = 0; r < 2; r++)
    {
        mp_hlimb_t u0, u1, u2, u3, v0, v1, v2, v3, z0, z1, z2, z3;
        
        u0 = x[0];
        u1 = x[1];
        u2 = x[2];
        u3 = x[3];

        v0 = u0 + u2;
        v0 -= (v0 >= p2) ? (p2) : 0;
        v2 = u0 - u2;
        v2 += ((mp_hlimb_signed_t) v2 < 0) ? (p2) : 0;
        v1 = u1 + u3;
        v1 -= (v1 >= p2) ? (p2) : 0;
        v3 = u1 - u3 + p2;
        v3 = mul_mod_precon_unreduced_32(v3, w, p, iw);
        z0 = v0 + v1;
        z0 -= (z0 >= p2) ? (p2) : 0;
        z0 -= (z0 >= p) ? (p) : 0;
        z1 = v0 - v1;
        z1 += ((mp_hlimb_signed_t) z1 < 0) ? (p2) : 0;
        z1 -= (z1 >= p) ? (p) : 0;
        z2 = v2 + v3;
        z2 -= (z2 >= p2) ? (p2) : 0;
        z2 -= (z2 >= p) ? (p) : 0;
        z3 = v2 - v3;
        z3 += ((mp_hlimb_signed_t) z3 < 0) ? (p2) : 0;
        z3 -= (z3 >= p) ? (p) : 0;
        
        x[0] = z0;
        x[1] = z1;
        x[2] = z2;
        x[3] = z3;
        x += 4;
    }

}



/*------------------------------------------------------------*/
/* in-place FFT in size 2^k, k >= 4                           */
/* powers_w_in are the roots of 1 needed                      */
/* input/output is a vector of 32 bit unsigned ints           */
/*------------------------------------------------------------*/
static void _fft_k(mp_hlimb_t * x,
                   mp_hlimb_t * powers_w, mp_hlimb_t * i_powers_w,
                   const nmod_t mod, const ulong k)
{
    /* N=size of block, M=number of blocks */
    ulong N, M;
    slong r, i, nb;
    mp_hlimb_t p, p2, w, iw;
    mp_hlimb_t  *x0, *x1;
    __m256i avx_p, avx_p2, avx_p2_minus_1;
        
    N = 1L << k;
    M = 1;
    p = (mp_hlimb_t) mod.n;
    p2 = p + p;

    avx_p = _mm256_set1_epi32(p);
    avx_p2 = _mm256_set1_epi32(p + p);
    avx_p2_minus_1 = _mm256_set1_epi32(p2 - 1);
    
    // the first layer is different: no need to reduce the sums
    x0 = x;
    x1 = x + N/2;
    nb = N/16;
    
    for (i = 0; i < nb; i++)
    {
        __m256i avx_u0, avx_u1, avx_i_w, avx_w, avx_t0, avx_t1, avx_q, avx_t2;
        avx_u0 = _mm256_load_si256((__m256i const*) x0);
        avx_u1 = _mm256_load_si256((__m256i const*) x1);
        avx_w = _mm256_load_si256((__m256i const*) (powers_w + 8*i));
        avx_i_w = _mm256_load_si256((__m256i const*) (i_powers_w + 8*i));
        
        avx_t0 = _mm256_add_epi32(avx_u0, avx_u1);
        
        avx_t1 = _mm256_sub_epi32(_mm256_add_epi32(avx_u0, avx_p2), avx_u1);
        avx_q = mm256_mulhi_epi32(avx_t1, avx_i_w);
        avx_t2 = _mm256_sub_epi32(_mm256_mullo_epi32(avx_t1, avx_w), _mm256_mullo_epi32(avx_q, avx_p));
        
        _mm256_store_si256((__m256i*) x0, avx_t0);
        _mm256_store_si256((__m256i*) x1, avx_t2);
        x0 += 8;
        x1 += 8;
    }
    
    powers_w += N / 2;
    i_powers_w += N / 2;
    N /= 2;
    M = 2;
    
    // all middle layers, N >= 16
    for (; N > 8; N /= 2, M *= 2)
    {
        x0 = x;
        x1 = x + N/2;
        
        for (r = 0; r < M; r++)
        {
            slong nb;
            nb = N/16;
            
            for (i = 0; i < nb; i++)
            {
                __m256i avx_u0, avx_u1, avx_i_w, avx_w, avx_t0, avx_t1, avx_q, avx_t2, avx_cmp;
                avx_u0 = _mm256_load_si256((__m256i const*) x0);
                avx_u1 = _mm256_load_si256((__m256i const*) x1);
                avx_w = _mm256_load_si256((__m256i const*) (powers_w + 8*i));
                avx_i_w = _mm256_load_si256((__m256i const*) (i_powers_w + 8*i));
                
                avx_t0 = _mm256_add_epi32(avx_u0, avx_u1);
                avx_cmp = _mm256_cmpgt_epi32(avx_t0, avx_p2_minus_1);
                avx_t0 = _mm256_sub_epi32(avx_t0, _mm256_and_si256(avx_cmp, avx_p2));
                
                avx_t1 = _mm256_sub_epi32(_mm256_add_epi32(avx_u0, avx_p2), avx_u1);
                avx_q = mm256_mulhi_epi32(avx_t1, avx_i_w);
                avx_t2 = _mm256_sub_epi32(_mm256_mullo_epi32(avx_t1, avx_w), _mm256_mullo_epi32(avx_q, avx_p));
                
                _mm256_store_si256((__m256i*) x0, avx_t0);
                _mm256_store_si256((__m256i*) x1, avx_t2);
                x0 += 8;
                x1 += 8;
            }
            x0 += N / 2;
            x1 += N / 2;
        }
        powers_w += N / 2;
        i_powers_w += N / 2;
    }
    
    // now N = 8
    x0 = x;
    x1 = x + 4;

    for (r = 0; r < M; r++)
    {
        for (i = 0; i < 4; i++)
        {
            mp_hlimb_t u0, u1, t1;
            mp_hlimb_signed_t t0;
            u0 = x0[i];
            u1 = x1[i];
            t0 = u0 + u1 - p2;
            x0[i] = (t0 < 0) ? t0+p2 : t0;
            t1 = (u0 + p2) - u1;
            x1[i] = mul_mod_precon_unreduced_32(t1, powers_w[i], p, i_powers_w[i]);
        }
        x0 += 8;
        x1 += 8;
    }
    powers_w += 4;
    i_powers_w += 4;
    M *= 2;
    
    // last two layers
    // not much SIMD here
    w = powers_w[1];
    iw = i_powers_w[1];

    for (r = 0; r < M; r++)
    {
        mp_hlimb_t u0, u1, u2, u3, v0, v1, v2, v3, z0, z1, z2, z3;
        
        u0 = x[0];
        u1 = x[1];
        u2 = x[2];
        u3 = x[3];

        v0 = u0 + u2;
        v0 -= (v0 >= p2) ? (p2) : 0;
        v2 = u0 - u2;
        v2 += ((mp_hlimb_signed_t) v2 < 0) ? (p2) : 0;
        v1 = u1 + u3;
        v1 -= (v1 >= p2) ? (p2) : 0;
        v3 = u1 - u3 + p2;
        v3 = mul_mod_precon_unreduced_32(v3, w, p, iw);
        z0 = v0 + v1;
        z0 -= (z0 >= p2) ? (p2) : 0;
        z0 -= (z0 >= p) ? (p) : 0;
        z1 = v0 - v1;
        z1 += ((mp_hlimb_signed_t) z1 < 0) ? (p2) : 0;
        z1 -= (z1 >= p) ? (p) : 0;
        z2 = v2 + v3;
        z2 -= (z2 >= p2) ? (p2) : 0;
        z2 -= (z2 >= p) ? (p) : 0;
        z3 = v2 - v3;
        z3 += ((mp_hlimb_signed_t) z3 < 0) ? (p2) : 0;
        z3 -= (z3 >= p) ? (p) : 0;
        
        x[0] = z0;
        x[1] = z1;
        x[2] = z2;
        x[3] = z3;
        x += 4;
    }
}


/*------------------------------------------------------------*/
/* fft evaluation                                             */
/* returns x[i] = poly(w^bitreverse_n(i)), n=2^k              */
/* x must have length >= n                                    */
/*------------------------------------------------------------*/
void nmod_avx2_32_fft_evaluate(mp_ptr x, const nmod_poly_t poly, const nmod_avx2_32_fft_t F, const ulong k)
{
    slong i, N;
    mp_hlimb_t * x_32;

    N = 1L << k;
    for (i = 0; i < N; i++)
    {
        x[i] = nmod_poly_get_coeff_ui(poly, i);
    }
    
    // special cases
    switch(k)
    {
        case 0:
            break;
        case 1:
            _fft_1(x, F->mod);
            break;
        case 2:
            _fft_2(x, F->powers_w[2][1], F->mod);
            break;
        case 3:
            // too many copies
            x_32 = (mp_hlimb_t *) aligned_alloc(32, 32);
            for (i = 0; i < 8; i++)
            {
                x_32[i] = (mp_hlimb_t) x[i];
            }
            
            _fft_3(x_32, F->powers_w[k], F->i_powers_w[k], F->mod);
            
            for (i = 0; i < 8; i++)
            {
                x[i] = (mp_limb_t) x_32[i];
            }
            
            free(x_32);
            break;
        default:
            x_32 = (mp_hlimb_t *) aligned_alloc(32, 4 * N);
            for (i = 0; i < N; i++)
            {
                x_32[i] = (mp_hlimb_t) x[i];
            }
            
            _fft_k(x_32, F->powers_w[k], F->i_powers_w[k], F->mod, k);
            
            for (i = 0; i < N; i++)
            {
                x[i] = (mp_limb_t) x_32[i];
            }
            
            free(x_32);
            break;
    }
}


#endif
