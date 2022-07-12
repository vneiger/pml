#include "nmod_poly_extra.h"
#include "nmod_vec_extra.h"

#ifdef HAS_AVX2
#include <immintrin.h>

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* using AVX2 to do inverse FFT nmod_poly with p < 2^32       */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/*  in-place size 1 inverse butterfly                         */
/*------------------------------------------------------------*/
static inline void _inv_fft_32_1(mp_hlimb_t * x, const nmod_t mod)
{
    mp_limb_t u0, u1, t0, t1;
    
    u0 = x[0];
    u1 = x[1];
    BUTTERFLY32(t0, t1, u0, u1, mod);
    x[0] = t0;
    x[1] = t1;
}

/*------------------------------------------------------------*/
/*  in-place size 2 inverse butterfly                         */
/*------------------------------------------------------------*/
static inline void _inv_fft_32_2(mp_hlimb_t * x, const mp_limb_t w, mp_hlimb_t iw, const nmod_t mod)
{
    mp_limb_t u0, u1, u2, u3, v0, v1, v2, v3, z0, z1, z2, z3;
    
    u0 = x[0];
    u1 = x[1];
    u2 = x[2];
    u3 = x[3];

    BUTTERFLY32(v0, v1, u0, u1, mod);
    BUTTERFLY32(v2, v3, u2, u3, mod);
    BUTTERFLY32(z0, z1, v0, v2, mod);
    v3 = mul_mod_precon_32(v3, w, mod.n, iw);
    BUTTERFLY32(z2, z3, v1, v3, mod);

    x[0] = z0;
    x[1] = z2;
    x[2] = z1;
    x[3] = z3;
}

/*------------------------------------------------------------*/
/* in-place inverse FFT in size 2^3                           */
/* powers_w_in are the roots of 1 needed                      */
/* input/output is a vector of 32 bit unsigned ints           */
/*------------------------------------------------------------*/
static void _inv_fft_32_3(mp_hlimb_t * x, mp_hlimb_t * powers_inv_w, mp_hlimb_t * i_powers_inv_w, const nmod_t mod)
{
    mp_hlimb_t w, iw, p, p2;
    ulong r, i;
    mp_hlimb_t * x0, * x1, * x2;

    // new 
    powers_inv_w = powers_inv_w + 4;
    i_powers_inv_w = i_powers_inv_w + 4;
    
    w = powers_inv_w[1];
    iw = i_powers_inv_w[1];

    p = mod.n;
    p2 = p+p;

    // first two layers: no reduction
    // not much SIMD here
    x2 = x;
    for (r = 0; r < 2; r++, x2 += 4)
    {
        mp_limb_t u0, u1, u2, u3, v0, v1, v2, v3;
        u0 = x2[0];
        u1 = x2[1];
        u2 = x2[2];
        u3 = x2[3];
        
        v0 = u0 + u1; // [0,2p)
        v2 = u2 + u3; // [0,2p)
        v1 = u0 - u1 + p; // [0,2p)
        v3 = u2 - u3 + p; // [0,2p)
        v3 = mul_mod_precon_unreduced_32(v3, w, p, iw); // [0,2p)
        
        x2[0] = v0 + v2; // [0,4p)
        x2[1] = v1 + v3; // [0,4p)
        x2[2] = (v0 + p2) - v2; // [0,4p)
        x2[3] = (v1 + p2) - v3; // [0, 4p) 
    }
    
    powers_inv_w -= 4;
    i_powers_inv_w -= 4;
    /* powers_inv_w += 2; */
    /* i_powers_inv_w += 2; */

    // last layer, with full reduction
    x0 = x;
    x1 = x + 4;
   
    for (i = 0; i < 4; i++)
    {
        mp_hlimb_t u0, u1, q;
        mp_hlimb_signed_t z;
        z = x0[i] - p2;
        u0 = (z < 0) ? (mp_hlimb_t)(z + p2) : (mp_hlimb_t)z;
        u1 = mul_mod_precon_unreduced_32(x1[i], powers_inv_w[i], p, i_powers_inv_w[i]);  // [0,2p)
        q = u0 + u1;
        q -= (q >= p2) ? p2 : 0;
        q -= (q >= p) ? p : 0;
        x0[i] = q;
        q = (u0 + p2) - u1;
        q -= (q >= p2) ? p2 : 0;
        q -= (q >= p) ? p : 0;
        x1[i] = q;
    }
}



/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* using AVX2 to do inverse FFT nmod_poly with p < 2^32       */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
static void _inv_fft_32_k(mp_hlimb_t * x, mp_hlimb_t * powers_inv_w, mp_hlimb_t * i_powers_inv_w,
                          const nmod_t mod, const ulong k)
{
    mp_hlimb_t w, iw, p, p2;
    ulong N, M, r, i, nb;
    mp_hlimb_t * x0, * x1, * x2;
    __m256i avx_p, avx_p_minus_1, avx_p2, avx_p2_minus_1;

    N = 4;
    M = 1L << (k-2);    // number of blocks of size 4 = n/4

    powers_inv_w = powers_inv_w + (1 << k) - 4;
    i_powers_inv_w = i_powers_inv_w + (1 << k) - 4;

    w = powers_inv_w[1];
    iw = i_powers_inv_w[1];

    p = mod.n;
    p2 = p+p;
    avx_p = _mm256_set1_epi32(p);
    avx_p_minus_1 = _mm256_set1_epi32(p - 1);
    avx_p2 = _mm256_set1_epi32(p + p);
    avx_p2_minus_1 = _mm256_set1_epi32(p + p - 1);

    // first two layers: no reduction
    // not much SIMD here
    x2 = x;
    for (r = 0; r < M; r++, x2 += 4)
    {
        mp_limb_t u0, u1, u2, u3, v0, v1, v2, v3;
        u0 = x2[0];
        u1 = x2[1];
        u2 = x2[2];
        u3 = x2[3];
        
        v0 = u0 + u1; // [0,2p)
        v2 = u2 + u3; // [0,2p)
        v1 = u0 - u1 + p; // [0,2p)
        v3 = u2 - u3 + p; // [0,2p)
        v3 = mul_mod_precon_unreduced_32(v3, w, p, iw); // [0,2p)
        
        x2[0] = v0 + v2; // [0,4p)
        x2[1] = v1 + v3; // [0,4p)
        x2[2] = (v0 + p2) - v2; // [0,4p)
        x2[3] = (v1 + p2) - v3; // [0, 4p) 
    }
    
    powers_inv_w -= 4;
    i_powers_inv_w -= 4;
    N = 8;
    M /= 2;

    // now N=8, M=2^k/8
    x0 = x;
    x1 = x + 4;
    for (r = 0; r < M; r++, x0 += 8, x1 += 8)
    {
        // N is still too small for AVX2
        for (i = 0; i < 4; i++)
        {
            mp_hlimb_t u0, u1;
            mp_hlimb_signed_t z;
            z = x0[i] - p2; // (-2p,2p)
            u0 = (z < 0) ? (mp_hlimb_t)(z + p2) : (mp_hlimb_t)z; // [0,2p)
            u1 = mul_mod_precon_unreduced_32(x1[i], powers_inv_w[i], p, i_powers_inv_w[i]);  // [0,2p)
            x0[i] = u0 + u1; // [0,4p)
            x1[i] = (u0 + p2) - u1; // [0,4p)
        }
    }
    powers_inv_w -= 8;
    i_powers_inv_w -= 8;
    N = 16;
    M /= 2;
    
    // now N=16, M=2^k/8
    // middle layers
    for (; M >= 2; N *= 2, M /= 2)
    {
        x0 = x;
        x1 = x + N/2;
        for (r = 0; r < M; r++, x0 += N/2, x1 += N/2)
        {
            // N/2 >= 8 so nb >= 1
            nb = N/16;
            for (i = 0; i < nb; i++)
            {
                __m256i avx_x0, avx_x1, avx_i_w, avx_w, avx_q, avx_cmp, avx_z;
                avx_x0 = _mm256_load_si256((__m256i const*) x0);
                avx_x1 = _mm256_load_si256((__m256i const*) x1);
                avx_w = _mm256_load_si256((__m256i const*) (powers_inv_w + 8*i));
                avx_i_w = _mm256_load_si256((__m256i const*) (i_powers_inv_w + 8*i));

                avx_cmp = _mm256_cmpgt_epi32(avx_x0, avx_p2_minus_1);
                avx_x0 = _mm256_sub_epi32(avx_x0, _mm256_and_si256(avx_cmp, avx_p2));

                avx_q = mm256_mulhi_epi32(avx_x1, avx_i_w);
                avx_z = _mm256_sub_epi32(_mm256_mullo_epi32(avx_x1, avx_w), _mm256_mullo_epi32(avx_q, avx_p));
                
                avx_x1 = _mm256_sub_epi32(_mm256_add_epi32(avx_x0, avx_p2), avx_z);
                avx_x0 = _mm256_add_epi32(avx_x0, avx_z);

                _mm256_store_si256((__m256i*) x0, avx_x0);
                _mm256_store_si256((__m256i*) x1, avx_x1);

                x0 += 8;
                x1 += 8;
            }
        }
        powers_inv_w -= N;
        i_powers_inv_w -= N;
    }

    // last layer, with full reduction
    x0 = x;
    x1 = x + N/2;
    nb = N/16;  // nb >= 1
    for (i = 0; i < nb; i++)
    {
        __m256i avx_x0, avx_x1, avx_i_w, avx_w, avx_q, avx_cmp, avx_z;
        avx_x0 = _mm256_load_si256((__m256i const*) x0);
        avx_x1 = _mm256_load_si256((__m256i const*) x1);
        avx_w = _mm256_load_si256((__m256i const*) (powers_inv_w + 8*i));
        avx_i_w = _mm256_load_si256((__m256i const*) (i_powers_inv_w + 8*i));
        
        avx_cmp = _mm256_cmpgt_epi32(avx_x0, avx_p2_minus_1);
        avx_x0 = _mm256_sub_epi32(avx_x0, _mm256_and_si256(avx_cmp, avx_p2));

        avx_q = mm256_mulhi_epi32(avx_x1, avx_i_w);  
        avx_z = _mm256_sub_epi32(_mm256_mullo_epi32(avx_x1, avx_w), _mm256_mullo_epi32(avx_q, avx_p));
        
        avx_x1 = _mm256_sub_epi32(_mm256_add_epi32(avx_x0, avx_p2), avx_z);
        avx_cmp = _mm256_cmpgt_epi32(avx_x1, avx_p2_minus_1);
        avx_x1 = _mm256_sub_epi32(avx_x1, _mm256_and_si256(avx_cmp, avx_p2));
        avx_cmp = _mm256_cmpgt_epi32(avx_x1, avx_p_minus_1);
        avx_x1 = _mm256_sub_epi32(avx_x1, _mm256_and_si256(avx_cmp, avx_p));

        avx_x0 = _mm256_add_epi32(avx_x0, avx_z);
        avx_cmp = _mm256_cmpgt_epi32(avx_x0, avx_p2_minus_1);
        avx_x0 = _mm256_sub_epi32(avx_x0, _mm256_and_si256(avx_cmp, avx_p2));
        avx_cmp = _mm256_cmpgt_epi32(avx_x0, avx_p_minus_1);
        avx_x0 = _mm256_sub_epi32(avx_x0, _mm256_and_si256(avx_cmp, avx_p));
        
        _mm256_store_si256((__m256i*) x0, avx_x0);
        _mm256_store_si256((__m256i*) x1, avx_x1);

        x0 += 8;
        x1 += 8;
    }
}

/*---------------------------------------------------------*/
/* reduces x mod (X^b2-1), assuming x has length a         */
/* assumes b2 divides a                                    */
/* all calculations are mod p                              */
/*---------------------------------------------------------*/
static inline
void fold_minus(mp_hlimb_t * tmp, mp_hlimb_t * x, ulong b2, ulong a, mp_hlimb_t p)
{
    ulong i, j, a8;
    mp_hlimb_t * acc, * tmp2;
    mp_hlimb_t acc0, acc1;
    __m256i avx_acc, avx_p, avx_p_minus_1;

    avx_p = _mm256_set1_epi32(p);
    avx_p_minus_1 = _mm256_set1_epi32(p - 1);
    
    switch(b2)
    {
        case 2:
            a8 = (a >> 3) << 3;
            acc = (mp_hlimb_t *) aligned_alloc(32, 32);
            avx_acc = _mm256_set1_epi32(0);
            for (i = 0; i < a8; i+= 8)
            {
                avx_acc = mm256_add_mod(avx_acc, _mm256_load_si256((__m256i const*) (x+i)), avx_p, avx_p_minus_1);
            }
            _mm256_store_si256((__m256i*) acc, avx_acc);
                        
            acc0 = n_addmod(n_addmod(n_addmod(acc[0], acc[2], p), acc[4], p), acc[6], p); // in [0..p)
            acc1 = n_addmod(n_addmod(n_addmod(acc[1], acc[3], p), acc[5], p), acc[7], p); // in [0..p)

            while (i < a)
            {
                acc0 = n_addmod(acc0, x[i], p);
                i++;
                acc1 = n_addmod(acc1, x[i], p);
                i++;
            }
            tmp[0] = acc0;
            tmp[1] = acc1;
            free(acc);
            break;

        case 4:
            if (a == 4)
            {
                tmp[0] = x[0];
                tmp[1] = x[1];
                tmp[2] = x[2];
                tmp[3] = x[3];
                break;
            }
            else
            {
                acc = (mp_hlimb_t *) aligned_alloc(32, 32);
                avx_acc = _mm256_set1_epi32(0);

                for (i = 0; i < a; i += 8)
                {
                    avx_acc = mm256_add_mod(avx_acc, _mm256_load_si256((__m256i const*) (x+i)), avx_p, avx_p_minus_1);
                }
                _mm256_store_si256((__m256i*) acc, avx_acc);
                
                tmp[0] = n_addmod(acc[0], acc[4], p);
                tmp[1] = n_addmod(acc[1], acc[5], p);
                tmp[2] = n_addmod(acc[2], acc[6], p);
                tmp[3] = n_addmod(acc[3], acc[7], p);

                free(acc);
            }
            break;

        default:
            // now b2 is at least 8
            tmp2 = tmp;
            for (i = 0; i < b2;  i += 8, x += 8, tmp2 += 8)
            {
                _mm256_store_si256((__m256i*) tmp2, _mm256_load_si256((__m256i const*) x));
                
            }
            for (i = b2; i < a; i += b2)
            {
                tmp2 = tmp;
                for (j = 0; j < b2; j += 8, x += 8, tmp2 += 8)
                {
                    _mm256_store_si256((__m256i*) tmp2,
                                       mm256_add_mod(_mm256_load_si256((__m256i const*) tmp2),
                                                     _mm256_load_si256((__m256i const*) x), avx_p, avx_p_minus_1));
                }
            }
            break;
    }
}

/*-----------------------------------------------------------*/
/* applies the CRT map to x in place                         */
/* x has length n                                            */
/* input is in [0,p), output is in [0,p)                     */
/* tmp is a temporary workspace, of length at least 2n       */
/*-----------------------------------------------------------*/
static inline
void CRT(mp_hlimb_t * x, mp_hlimb_t * tmp, ulong n, mp_hlimb_t p)
{
    ulong a, b, b2, b8, n2, n8, i, nn;
    mp_hlimb_t half;
    mp_hlimb_t * tmp2, * x2;
    __m256i avx_p, avx_p_minus_1, avx_half;

    nn = n;
    half = 1 + (p >> 1);

    avx_p = _mm256_set1_epi32(p);
    avx_p_minus_1 = _mm256_set1_epi32(p - 1);
    avx_half = _mm256_set1_epi32(half);
    
        
    if (nn <= 1)
    {
        return;
    }

    a = 1L;
    n2 = n >> 1;
    while (a <= n2)
    {
        a <<= 1;
    }
    
    while (n != a)
    {
        n = n - a;
        b = 1L;
        n2 = (n >> 1);
        while (b <= n2)
        {
            b <<= 1;
        }
        b2 = b << 1;

        fold_minus(tmp, x, b2, a, p);
        
        b8 = (b >> 3) << 3;
        tmp2 = tmp + b;
        x2 = x + a;
        for (i = 0; i < b8; i+= 8, tmp2 += 8, x2 += 8)
        {
            // here, b is at least 8
            __m256i avx_b;
            avx_b = _mm256_load_si256((__m256i const*) tmp2);
            _mm256_store_si256((__m256i*) x2,
                               mm256_add_mod(mm256_add_mod(avx_b, avx_b, avx_p, avx_p_minus_1),
                                             _mm256_load_si256((__m256i const*) x2), avx_p, avx_p_minus_1));
        }

        // less than 8 passes
        for (; i < b; i++)
        {
            x[a+i] = n_addmod(n_addmod(tmp[b+i], tmp[b+i], p), x[a+i], p);
        }

        tmp += b2;
        x += a;
        a = b;
    }

    while (n != nn)
    {
        b = a;
        a <<= 1;
        b2 = a;
        while (!(a & nn))
        {
            a <<= 1;
        }
        tmp -= b2;
        x -= a;

        b8 = (b >> 3) << 3;
        for (i = 0; i < b8; i+= 8)
        {
            __m256i u, u2;
            u = mm256_add_mod(_mm256_load_si256((__m256i const*) (tmp + i)),
                              _mm256_load_si256((__m256i const*) (tmp + b + i)),
                              avx_p, avx_p_minus_1);
            u = mm256_sub_mod(_mm256_load_si256((__m256i const*) (x + a + i)), u, avx_p, avx_p_minus_1);
            // u div 2
            u2 = _mm256_srli_epi32(u, 1);
            // u + (u mod 2 = 1 ? half : 0) gives u/2 mod p
            u = _mm256_mask_add_epi32(u2, _mm256_movemask_ps(_mm256_castsi256_ps(_mm256_slli_epi32(u, 31))), u2, avx_half);
            _mm256_store_si256((__m256i*) (x + i),
                               mm256_add_mod(_mm256_load_si256((__m256i const*) (x + i)), u, avx_p, avx_p_minus_1));
            _mm256_store_si256((__m256i*) (x + a + i), u);

        }

        // less than 8 passes
        for (; i < b; i++)
        {
            mp_hlimb_t u;
            u = n_addmod(tmp[i], tmp[b+i], p);
            u = n_submod(x[a+i], u, p);
            x[a+i] = (u >> 1) + (u & 1 ? half : 0);   // u/2 mod p, no multiplication
            /* x[a+i] = (u >> 1) + (u & 1)*half; */  // u/2 mod p 
            x[i] = n_addmod(x[i], x[a+i], p);
        }

        // now from b to n
        b8 = ((b + 7) >> 3) << 3;   // next multiple of 8
        n8 = (n >> 3) << 3;
        if (b8 > n8)
        {
            // less than 4 passes
            for (i = b; i < n; i++)
            {
                mp_hlimb_t u;
                u = x[a+i];
                x[a+i] = (u >> 1) + (u & 1 ? half : 0);   // u/2 mod p
                x[i] = n_addmod(x[i], x[a+i], p);
            }
        }
        else
        {
            for (i = b8; i < n8; i+=8)
            {
                __m256i u, u2;
                u = _mm256_load_si256((__m256i const*) (x + a + i));
                // u div 2
                u2 = _mm256_srli_epi32(u, 1);
                // u + (u mod 2 = 1 ? half : 0) gives u/2 mod p
                u = _mm256_mask_add_epi32(u2, _mm256_movemask_ps(_mm256_castsi256_ps(_mm256_slli_epi32(u, 31))), u2, avx_half);
                _mm256_store_si256((__m256i*) (x + i),
                                   mm256_add_mod(_mm256_load_si256((__m256i const*) (x + i)), u, avx_p, avx_p_minus_1));
                _mm256_store_si256((__m256i*) (x + a + i), u);
            }
            // less than 8 passes
            for (i = n8; i < n; i++)
            {
                mp_hlimb_t u;
                u = x[a+i];
                x[a+i] = (u >> 1) + (u & 1 ? half : 0);   // u/2 mod p
                x[i] = n_addmod(x[i], x[a+i], p);
            }
        }

        n = n + a;
    }
}


/*------------------------------------------------------------*/
/* inverse fft                                                */
/* given x[i] = poly(w^bitreverse(i,2^k)), returns poly       */
/*------------------------------------------------------------*/
void nmod_avx2_32_fft_interpolate(nmod_poly_t poly, mp_srcptr x, const nmod_32_fft_t F, const ulong k)
{
    ulong i, N;
    mp_hlimb_t * x_32;

    N = 1L << k;
    nmod_poly_fit_length(poly, N);
    poly->length = N;
    x_32 = (mp_hlimb_t *) aligned_alloc(32, 4*N > 32 ? 4*N : 32);

   
    for (i = 0; i < N; i++)
    {
        x_32[i] = (mp_hlimb_t) x[i];
    }

    switch(k)
    {
        case 0:
            break;
        case 1:
            _inv_fft_32_1(x_32, F->mod);
            break;
        case 2:
            _inv_fft_32_2(x_32, F->powers_inv_w_t[2][1], F->i_powers_inv_w_t[2][1], F->mod);
            break;
        case 3:
            _inv_fft_32_3(x_32, F->powers_inv_w_t[3], F->i_powers_inv_w_t[3], F->mod);
            break;
        default:
            _inv_fft_32_k(x_32, F->powers_inv_w_t[k], F->i_powers_inv_w_t[k], F->mod, k);
            break;
    }

    _nmod_vec_avx2_32_scalar_mul(x_32, x_32, N, F->powers_inv_2[k], F->i_powers_inv_2[k], F->mod);
    for (i = 0; i < N; i++)
    {
        poly->coeffs[i] = (mp_limb_t) x_32[i];
    }
    free(x_32);
    
    _nmod_poly_normalise(poly);
}


/*------------------------------------------------------------*/
/* inverse tft                                                */
/*------------------------------------------------------------*/
void nmod_avx2_32_tft_interpolate(nmod_poly_t poly, mp_srcptr x, const nmod_32_fft_t F, const ulong N)
{
    ulong i, nn, k, aa;
    mp_hlimb_t * wk, * wk2, * tmp, * powers, * i_powers;
    mp_hlimb_t p;
    
    nmod_poly_fit_length(poly, N);
    poly->length = N;

    wk = (mp_hlimb_t *) aligned_alloc(32, 4*N > 32 ? 4*N : 32);
    tmp = (mp_hlimb_t *) aligned_alloc(32, 8*N > 32 ? 8*N : 32);
    for (i = 0; i < N; i++)
    {
        wk[i] = x[i];
        tmp[i] = 0;
        tmp[i+N] = 0;
    }
    wk2 = wk;
    p = F->mod.n;
    
    nn = N;
    k = 0;
    aa = 1;
    
    while (aa <= nn/2)
    {
        k++;
        aa <<= 1;
    }

    do
    {
        switch(k)
        {
            case 0:
                break;
            case 1:
                _inv_fft_32_1(wk, F->mod);
                break;
            case 2:
                _inv_fft_32_2(wk, F->powers_inv_w_t[2][1], F->i_powers_inv_w_t[2][1], F->mod);
                break;
            case 3:
                _inv_fft_32_3(wk, F->powers_inv_w_t[3], F->i_powers_inv_w_t[3], F->mod);
                break;
            default:
                _inv_fft_32_k(wk, F->powers_inv_w_t[k], F->i_powers_inv_w_t[k], F->mod, k);
                break;
        }
        
        powers = F->powers_inv_w_over_2[k];
        i_powers = F->i_powers_inv_w_over_2[k];
        _nmod_vec_avx2_32_hadamard_mul(wk, wk, powers, i_powers, aa, F->mod);

        wk += aa;
        nn = nn-aa;
        
        aa = 1;
        k = 0;
        while (aa <= nn/2)
        {
            k++;
            aa <<= 1;
        }
    }
    while (nn != 0);

    wk = wk2;
    CRT(wk, tmp, N, p);

    for (i = 0; i < N; i++)
    {
        poly->coeffs[i] = wk[i];
    }
    
    free(wk);
    free(tmp);
    _nmod_poly_normalise(poly);
}


/*-----------------------------------------------------------*/
/* transpose tft in length N                                 */
/*-----------------------------------------------------------*/
void nmod_avx2_32_tft_evaluate_t(mp_ptr x, mp_srcptr A, const nmod_32_fft_t F, const ulong N)
{
    ulong n, a, b, b2, t, i, j, k, ell, lambda;
    mp_hlimb_t p;
    mp_hlimb_t * wk, * wk2, * powers_rho, * i_powers_rho;
    
    if (N == 0)
    {
        return;
    }
    
    if (N == 1)
    {
        x[0] = A[0];
        return;
    }

    wk = (mp_hlimb_t *) aligned_alloc(32, 8*N > 32 ? 8*N : 32);
    
    for (i = 0; i < N; i++)
    {
        wk[i] = A[i];
        wk[i + N] = 0;
    }

    p = F->mod.n;
    k = 0;
    a = 1;
    n = N;
    
    while (a <= n/2)
    {
        k++;
        a <<= 1;
    }

    wk2 = wk;
    do
    {
        // special cases
        switch(k)
        {
            case 0:
                break;
            case 1:
                _inv_fft_32_1(wk, F->mod);
                break;
            case 2:
                _inv_fft_32_2(wk, F->powers_w[2][1], F->i_powers_w[2][1], F->mod);
                break;
            case 3:
                _inv_fft_32_3(wk, F->powers_w[3], F->i_powers_w[3], F->mod);
                break;
            default:
                _inv_fft_32_k(wk, F->powers_w[k], F->i_powers_w[k], F->mod, k);
                break;
        }
        
        powers_rho = F->powers_w[k+1];
        i_powers_rho = F->i_powers_w[k+1];
        for (i = 0; i < a; i++)
        {
            wk[i] = mul_mod_precon_32(wk[i], powers_rho[i], p, i_powers_rho[i]);
        }
        wk += a;
        n = n-a;
        
        a = 1;
        k = 0;
        while (a <= n/2)
        {
            k++;
            a <<= 1;
        }
    }
    while (n != 0);
    wk = wk2;

    
    a = 1;
    k = 0;
    while (! (a & N))
    {
        k++;
        a <<= 1;
    }
    
    for (i = 0; i < a; i++)
    {
        if (wk[N-a+i] != 0) // neg_mod
        {
            wk[N+i] = p - wk[N-a+i];
        }
    }
    
    n = a;
    while (n != N)
    {
        b = a;
        ell = k;
        a <<= 1;
        k++;
        while (!(a & N))
        {
            a <<= 1;
            k++;
        }
        
        b2 = b << 1;
        lambda = 1L << (k-ell-1);

        n = n + a;
        t = b2;
        // we have to deal with i=0 separately
        // (since otherwise we would erase the source)
        switch(b)
        {
            case 1:
                for (i = 1; i < lambda; i++)
                {
                    mp_limb_t u, v;
                    u = wk[N-n+t];
                    v = wk[N-n+a];
                    wk[N-n+t] = n_addmod(v, u, p);
                    wk[N-n+a+t] = n_submod(v, u, p);
                    t++;
                    u = wk[N-n+t];
                    v = wk[N-n+a+1];
                    wk[N-n+t] = n_addmod(v, u, p);
                    wk[N-n+a+t] = n_submod(v, u, p);
                    t++;
                }
                
                mp_limb_t u, v;
                u = wk[N-n];
                v = wk[N-n+a];
                wk[N-n] = n_addmod(v, u, p);
                wk[N-n+a] = n_submod(v, u, p);
                u = wk[N-n+1];
                v = wk[N-n+a+1];
                wk[N-n+1] = n_addmod(v, u, p);
                wk[N-n+a+1] = n_submod(v, u, p);
                break;
            case 2:
            default:
                for (i = 1; i < lambda; i++)
                {
                    for (j = 0; j < b2; j++)
                    {
                        mp_limb_t u, v;
                        u = wk[N-n+t];
                        v = wk[N-n+a+j];
                        wk[N-n+t] = n_addmod(v, u, p);
                        wk[N-n+a+t] = n_submod(v, u, p);
                        t++;
                    }
                }
                for (j = 0; j < b2; j++)
                {
                    mp_limb_t u, v;
                    u = wk[N-n+j];
                    v = wk[N-n+a+j];
                    wk[N-n+j] = n_addmod(v, u, p);
                    wk[N-n+a+j] = n_submod(v, u, p);
                }
        }
    }
    
    for (i = 0; i < N; i++)
    {
        x[i] = wk[i];
    }

    free(wk);
}



#endif
