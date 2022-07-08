#include "nmod_poly_extra.h"

#ifdef HAS_AVX2
#include <immintrin.h>

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* using AVX2 to do FFT nmod_poly with p < 2^32               */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/*  in-place size 1 butterfly                                 */
/*  o1, o2 = i1+i2, i1-i2 mod p                               */
/*------------------------------------------------------------*/
static inline void _fft_32_1(mp_hlimb_t * x, const nmod_t mod)
{
    mp_hlimb_t u0, u1, t0, t1;
    
    u0 = x[0];
    u1 = x[1];
    BUTTERFLY32(t0, t1, u0, u1, mod);
    x[0] = t0;
    x[1] = t1;

}

/*------------------------------------------------------------*/
/*  in-place size 2 butterfly                                 */
/*  t1, t2, t3, t4 = i1+i3, i2+i4, i1-i3, w(i2-i4) mod p      */
/*  o1, o2, o3, o4 = t1+t2, t1-t2, t3+t4, t3-t4 mod p         */
/*------------------------------------------------------------*/
static inline void _fft_32_2(mp_hlimb_t * x, mp_hlimb_t w, mp_hlimb_t iw, const nmod_t mod)
{
    mp_hlimb_t u0, u1, u2, u3, v0, v1, v2, v3, z0, z1, z2, z3;

    u0 = x[0];
    u1 = x[1];
    u2 = x[2];
    u3 = x[3];
    
    BUTTERFLY32(v0, v2, u0, u2, mod);
    BUTTERFLY32(v1, v3, u1, u3, mod);
    v3 = mul_mod_precon_32(v3, w, mod.n, iw);
    BUTTERFLY32(z0, z1, v0, v1, mod);
    BUTTERFLY32(z2, z3, v2, v3, mod);
    
    x[0] = z0;
    x[1] = z1;
    x[2] = z2;
    x[3] = z3;
}

/*------------------------------------------------------------*/
/* in-place FFT in size 2^3                                   */
/* powers_w_in are the roots of 1 needed                      */
/* input/output is a vector of 32 bit unsigned ints           */
/*------------------------------------------------------------*/
static void _fft_32_3(mp_hlimb_t * x,
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
static void _fft_32_k(mp_hlimb_t * x,
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



    for (i = 0; i < N; i++)
    {
        if (x[i] >= p)
        {
            printf("input not reduced\n");
            exit(-1);
        }
    }
    
    avx_p = _mm256_set1_epi32(p);
    avx_p2 = _mm256_set1_epi32(p + p);
    avx_p2_minus_1 = _mm256_set1_epi32(p2 - 1);
    
    // the first layer is different: no need to reduce the sums
    x0 = x;
    x1 = x + N/2;
    nb = N/16;
    
    for (i = 0; i < nb; i++)
    {
        __m256i u0, u1, i_w, w, t0, t1, t2;
        u0 = _mm256_load_si256((__m256i const*) x0);
        u1 = _mm256_load_si256((__m256i const*) x1);
        w = _mm256_load_si256((__m256i const*) (powers_w + 8*i));
        i_w = _mm256_load_si256((__m256i const*) (i_powers_w + 8*i));
        
        t0 = _mm256_add_epi32(u0, u1); // [0..2p)
        t1 = mm256_sub_mod_unreduced(u0, u1, avx_p); // [0..2p)
        t2 = mm256_mul_mod_precon_unreduced(t1, w, avx_p, i_w); 
        
        _mm256_store_si256((__m256i*) x0, t0);
        _mm256_store_si256((__m256i*) x1, t2);
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
                __m256i u0, u1, i_w, w, t0, t1, t2;
                u0 = _mm256_load_si256((__m256i const*) x0);
                u1 = _mm256_load_si256((__m256i const*) x1);
                w = _mm256_load_si256((__m256i const*) (powers_w + 8*i));
                i_w = _mm256_load_si256((__m256i const*) (i_powers_w + 8*i));
                
                t0 = mm256_add_mod(u0, u1, avx_p2, avx_p2_minus_1); // [0..2p)
                t1 = mm256_sub_mod_unreduced(u0, u1, avx_p2); // [0..4p)
                t2 = mm256_mul_mod_precon_unreduced(t1, w, avx_p, i_w); // [0..2p)
                
                _mm256_store_si256((__m256i*) x0, t0);
                _mm256_store_si256((__m256i*) x1, t2);
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
void nmod_avx2_32_fft_evaluate(mp_ptr x, const nmod_poly_t poly, const nmod_32_fft_t F, const ulong k)
{
    slong i, N;
    mp_hlimb_t * x_32;

    N = 1L << k;
    x_32 = (mp_hlimb_t *) aligned_alloc(32, 4*N > 32 ? 4*N : 32);

    for (i = 0; i < N; i++)
    {
        x_32[i] = (mp_hlimb_t) nmod_poly_get_coeff_ui(poly, i);
    }

    // special cases
    switch(k)
    {
        case 0:
            break;
        case 1:
            _fft_32_1(x_32, F->mod);
            break;
        case 2:
            _fft_32_2(x_32, F->powers_w[2][1], F->i_powers_w[2][1], F->mod);
            break;
        case 3:
            _fft_32_3(x_32, F->powers_w[k], F->i_powers_w[k], F->mod);
            break;
        default:
            _fft_32_k(x_32, F->powers_w[k], F->i_powers_w[k], F->mod, k);
            break;
    }
    
    for (i = 0; i < N; i++)
    {
        x[i] = (mp_limb_t) x_32[i];
    }
    
    free(x_32);
}


/*------------------------------------------------------------*/
/* tft evaluation                                             */
/* returns a cyclotomic tft evaluation of x at n points       */
/* x must have length >= N                                    */
/* input and output reduced mod p                             */
/*------------------------------------------------------------*/
void nmod_avx2_32_tft_evaluate(mp_ptr x, const nmod_poly_t poly, const nmod_32_fft_t F, const ulong N)
{
    slong i, j, k, a, ell, b, t, b2, lambda, lambda8;
    mp_hlimb_t p, u, v, acc1, acc2, acc3, acc4;
    mp_hlimb_t  * x2, * x2_bak, * powers_rho, * i_powers_rho, * acc;
    __m256i avx_p, avx_p_minus_1;

    
    // N <= 1 -> easy case
    if (N <= 1)
    {
        for (i = 0; i < N; i++)
        {
            x[i] = nmod_poly_get_coeff_ui(poly, i);
        }
        return;
    }

    // else, loop through the reduction process
    p = F->mod.n;
    x2 = (mp_hlimb_t *) aligned_alloc(32, 2*4*(N >= 8 ? N : 8));
    x2_bak = x2;

    avx_p = _mm256_set1_epi32(p);
    avx_p_minus_1 = _mm256_set1_epi32(p - 1);
    acc = (mp_hlimb_t *) aligned_alloc(32, 32);
    
    for (i = 0; i < N; i++)
    {
        x2[i] = nmod_poly_get_coeff_ui(poly, i);
        x2[i + N] = 0;
    }
    
    ell = 0;
    b = 1;
    while (b <= N/2)
    {
        ell++;
        b <<= 1;
    }
    
    do
    {
        a = b;
        k = ell;
        
        b = a >> 1;
        ell = k-1;
        while (b > 0 && (! (N & b)))
        {
            b >>= 1;
            ell--;
        }

        t = 0;
        b2 = b << 1;

        powers_rho = F->powers_w[k+1];
        i_powers_rho = F->i_powers_w[k+1];
        
        // claim: if all entries in x2 are in [0..p) before the switch,
        // they remain in [0..p) after
        switch(b)
        {
            case 0:
                if (a >= 8)
                {
                    for (i = 0; i < a; i+=8)
                    {
                        __m256i x2_1, x2_2, rho, i_rho;
                        x2_1 = _mm256_load_si256((__m256i const*) (x2 + i)); // in [0..p) by assumption
                        x2_2 = _mm256_load_si256((__m256i const*) (x2 + i + a)); // in [0..p) by assumption
                        rho = _mm256_load_si256((__m256i const*) (powers_rho + i));
                        i_rho = _mm256_load_si256((__m256i const*) (i_powers_rho + i));
                        x2_1 = mm256_sub_mod_unreduced(x2_1, x2_2, avx_p);  // in [0..2p)
                        _mm256_store_si256((__m256i*) (x2 + i), mm256_mul_mod_precon(x2_1, rho, avx_p, avx_p_minus_1, i_rho));
                        // in [0..p) 
                    }
                }
                else
                {
                    for (i = 0; i < a; i++)
                    {
                        // intermediate sum in [0..2p)
                        x2[i] = mul_mod_precon_32(x2[i] + p - x2[i + a], powers_rho[i], p, i_powers_rho[i]);
                        // in [0..p)
                    }
                }
                break;
                

            case 1:
                // b = 1: unroll the inner loop
                if (a >= 8)
                {
                    lambda = 1L << (k-1);
                    lambda8 = ((2*lambda)/8)*8;
                    for (i = 0; i < 8; i++)
                    {
                        acc[i] = 0;
                    }
                    for (t = 0; t < lambda8; t+=8)
                    {
                        __m256i u, v, rho, i_rho, avx_acc;
                        u = _mm256_load_si256((__m256i const*) (x2 + t)); // in [0..p) by assumption
                        v = _mm256_load_si256((__m256i const*) (x2 + t + a)); // in [0..p) by assumption
                        rho = _mm256_load_si256((__m256i const*) (powers_rho + t));
                        i_rho = _mm256_load_si256((__m256i const*) (i_powers_rho + t));

                        // intermediate sum in [0..2p)
                        // result in [0..p)
                        _mm256_store_si256((__m256i*) (x2 + t),
                                           mm256_mul_mod_precon(mm256_sub_mod_unreduced(u, v, avx_p), rho, avx_p, avx_p_minus_1, i_rho));
                        
                        u = mm256_add_mod(u, v, avx_p, avx_p_minus_1); // in [0..p)
                        avx_acc = _mm256_load_si256((__m256i const*) acc); // in [0..p) (rec. assumption)
                        u = mm256_add_mod(u, avx_acc, avx_p, avx_p_minus_1); // in [0..p)
                        _mm256_store_si256((__m256i*) acc, u);
                    }
                    acc1 = n_addmod(n_addmod(n_addmod(acc[0], acc[2], p), acc[4], p), acc[6], p); // in [0..p)
                    acc2 = n_addmod(n_addmod(n_addmod(acc[1], acc[3], p), acc[5], p), acc[7], p); // in [0..p)

                    x2[a] = acc1;
                    x2[a+1] = acc2;
                }
                else // (1,2) or (2,4)
                {
                    u = x2[0];
                    v = x2[a];
                    x2[0] = mul_mod_precon_32(u + p - v, powers_rho[0], p, i_powers_rho[0]);
                    acc1 = n_addmod(u, v, p); // in [0..p)

                    u = x2[1];
                    v = x2[a+1];
                    x2[1] = mul_mod_precon_32(u + p - v, powers_rho[1], p, i_powers_rho[1]);
                    acc2 = n_addmod(u, v, p); // in [0..p)
                    
                    for (t = 2; t < a; t+=2)
                    {
                        u = x2[t];
                        v = x2[t + a];
                        x2[t] = mul_mod_precon_32(u + p - v, powers_rho[t], p, i_powers_rho[t]);
                        acc1 = n_addmod(n_addmod(u, v, p), acc1, p); // in [0..p)
                        u = x2[t + 1];
                        v = x2[t + a+1];
                        x2[t + 1] = mul_mod_precon_32(u + p - v, powers_rho[t + 1], p, i_powers_rho[t + 1]);
                        acc2 = n_addmod(n_addmod(u, v, p), acc2, p); // in [0..p)
                    }
                    x2[a] = acc1;
                    x2[a+1] = acc2;
                }
                break;

                // b = 2: unroll the inner loop
            case 2:
                lambda = 1L << (k-2);
                lambda8 = (lambda/2)*8;

                // AVX?
                for (i = 0; i < 8; i++)
                {
                    acc[i] = 0;
                }

                for (t = 0; t < lambda8; t+=8)
                {
                    __m256i u, v, rho, i_rho, avx_acc;
                    u = _mm256_load_si256((__m256i const*) (x2 + t));
                    v = _mm256_load_si256((__m256i const*) (x2 + t + a));
                    rho = _mm256_load_si256((__m256i const*) (powers_rho + t));
                    i_rho = _mm256_load_si256((__m256i const*) (i_powers_rho + t));
                    
                    _mm256_store_si256((__m256i*) (x2 + t),
                                       mm256_mul_mod_precon(mm256_sub_mod_unreduced(u, v, avx_p), rho, avx_p, avx_p_minus_1, i_rho));
                    u = mm256_add_mod(u, v, avx_p, avx_p_minus_1);
                    avx_acc = _mm256_load_si256((__m256i const*) acc);
                    u = mm256_add_mod(u, avx_acc, avx_p, avx_p_minus_1);
                    
                    _mm256_store_si256((__m256i*) acc, u);
                }

                acc1 = n_addmod(acc[0], acc[4], p);
                acc2 = n_addmod(acc[1], acc[5], p);
                acc3 = n_addmod(acc[2], acc[6], p);
                acc4 = n_addmod(acc[3], acc[7], p);

                if (t == 0)
                {
                    u = x2[t];
                    v = x2[a+t];
                    x2[t] = mul_mod_precon_32(u + p - v, powers_rho[t], p, i_powers_rho[t]);
                    acc1 = n_addmod(n_addmod(u, v, p), acc1, p);
                    t++;
                    u = x2[t];
                    v = x2[a+t];
                    x2[t] = mul_mod_precon_32(u + p - v, powers_rho[t], p, i_powers_rho[t]);
                    acc2 = n_addmod(n_addmod(u, v, p), acc2, p);
                    t++;
                    u = x2[t];
                    v = x2[a+t];
                    x2[t] = mul_mod_precon_32(u + p - v, powers_rho[t], p, i_powers_rho[t]);
                    acc3 = n_addmod(n_addmod(u, v, p), acc3, p);
                    t++;
                    u = x2[t];
                    v = x2[a+t];
                    x2[t] = mul_mod_precon_32(u + p - v, powers_rho[t], p, i_powers_rho[t]);
                    acc4 = n_addmod(n_addmod(u, v, p), acc4, p);
                    t++;
                }
                
                x2[a] = acc1;
                x2[a+1] = acc2;
                x2[a+2] = acc3;
                x2[a+3] = acc4;
                
                break;
            default:
                // now we know b >= 4 so b2 >= 8
                // i = 0 to 2^(k-ell-1) = a/(2b), j = 0 to 2b
                // i = 0 : special case
                for (; t < b2; t+=8)
                {
                    __m256i u, v, rho, i_rho;
                    u = _mm256_load_si256((__m256i const*) (x2 + t));
                    v = _mm256_load_si256((__m256i const*) (x2 + t + a));
                    rho = _mm256_load_si256((__m256i const*) (powers_rho + t));
                    i_rho = _mm256_load_si256((__m256i const*) (i_powers_rho + t));

                    // [0..p)
                    _mm256_store_si256((__m256i*) (x2 + t),
                                       mm256_mul_mod_precon(mm256_sub_mod_unreduced(u, v, avx_p), rho, avx_p, avx_p_minus_1, i_rho)); 
                    u = mm256_add_mod(u, v, avx_p, avx_p_minus_1); // [0..p)
                    _mm256_store_si256((__m256i*) (x2 + t + a), u);
                }
                
                lambda = 1L << (k-ell-1);
                for (i = 1; i < lambda; i++)
                {
                    for (j = 0; j < b2; j+=8, t+=8)
                    {
                        __m256i u, v, rho, i_rho;
                        u = _mm256_load_si256((__m256i const*) (x2 + t));
                        v = _mm256_load_si256((__m256i const*) (x2 + t + a));
                        rho = _mm256_load_si256((__m256i const*) (powers_rho + t));
                        i_rho = _mm256_load_si256((__m256i const*) (i_powers_rho + t));

                        // [0..p)
                        _mm256_store_si256((__m256i*) (x2 + t),
                                           mm256_mul_mod_precon(mm256_sub_mod_unreduced(u, v, avx_p), rho, avx_p, avx_p_minus_1, i_rho));
                        u = mm256_add_mod(u, v, avx_p, avx_p_minus_1);  // [0..p)
                        u = mm256_add_mod(u, _mm256_load_si256((__m256i const*) (x2 + j + a)), avx_p, avx_p_minus_1); 
                        _mm256_store_si256((__m256i*) (x2 + j + a), u); // [0..p)
                    }
                }
        }

        switch(k)
        {
            case 0:
                break;
            case 1:
                _fft_32_1(x2, F->mod);
                break;
            case 2:
                _fft_32_2(x2, F->powers_w[2][1], F->i_powers_w[2][1], F->mod);
                break;
            case 3:
                _fft_32_3(x2, F->powers_w[3], F->i_powers_w[3], F->mod);
                break;
            default:
                _fft_32_k(x2, F->powers_w[k], F->i_powers_w[k], F->mod, k);
                break;
        }
        
        x2 += a;
    }
    while (b != 0);
    
    x2 = x2_bak;
    for (i = 0; i < N; i++)
    {
        x[i] = x2[i];
    }
    free(x2);
    free(acc);
}




#endif
