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
    slong r, i;
    mp_hlimb_t * x0, * x1, * x2;

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
    
    powers_inv_w += 2;
    i_powers_inv_w += 2;

    // last layer, with full reduction
    x0 = x;
    x1 = x + 4;
   
    for (i = 0; i < 4; i++)
    {
        mp_hlimb_t u0, u1, q;
        mp_hlimb_signed_t z;
        z = x0[i] - p2;
        u0 = (z < 0) ? z + p2 : z;
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
    ulong N, M;
    slong r, i, nb;
    mp_hlimb_t * x0, * x1, * x2;
    __m256i avx_p, avx_p_minus_1, avx_p2, avx_p2_minus_1;

    N = 4;
    M = 1L << (k-2);    // number of blocks of size 4 = n/4
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
    
    powers_inv_w += 2;
    i_powers_inv_w += 2;
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
            u0 = (z < 0) ? z + p2 : z; // [0,2p)
            u1 = mul_mod_precon_unreduced_32(x1[i], powers_inv_w[i], p, i_powers_inv_w[i]);  // [0,2p)
            x0[i] = u0 + u1; // [0,4p)
            x1[i] = (u0 + p2) - u1; // [0,4p)
        }
    }
    powers_inv_w += 4;
    i_powers_inv_w += 4;
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
        powers_inv_w += N/2;
        i_powers_inv_w += N/2;
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
void fold_minus(mp_hlimb_t * tmp, mp_hlimb_t * x, slong b2, slong a, mp_hlimb_t p)
{
    slong i, j, k, a8;
    mp_hlimb_t * acc;
    mp_hlimb_t acc0, acc1;
    
    switch(b2)
    {
        case 2:
            acc = (mp_hlimb_t *) aligned_alloc(32, 32);
            // TODO: AVX
            acc[0] = acc[1] = acc[2] = acc[3] = acc[4] = acc[5] = acc[6] = acc[7] = 0;
            
            a8 = (a >> 3) << 3;
            for (i = 0; i < a8; i+= 8)
            {
                //TODO: AVX
                for (j = 0; j < 8; j++)
                {
                    acc[j] = n_addmod(acc[j], x[i+j], p);
                }
            }
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
                for (i = 0; i < 4; i++)
                {
                    tmp[i] = x[i];
                }
                break;
            }
            else
            {
                acc = (mp_hlimb_t *) aligned_alloc(32, 32);
                // TODO: AVX
                acc[0] = acc[1] = acc[2] = acc[3] = acc[4] = acc[5] = acc[6] = acc[7] = 0;

                for (i = 0; i < a; i += 8)
                {
                    // TODO: AVX
                    for (j = 0; j < 8; j++)
                    {
                        acc[j] = n_addmod(acc[j], x[i+j], p);
                    }
                }
                
                tmp[0] = n_addmod(acc[0], acc[4], p);
                tmp[1] = n_addmod(acc[1], acc[5], p);
                tmp[2] = n_addmod(acc[2], acc[6], p);
                tmp[3] = n_addmod(acc[3], acc[7], p);

                free(acc);
            }
            break;

        default:
            // now b2 is at least 8

            for (i = 0; i < b2;  i+= 8)
            {
                // TODO: AVX
                for (j = 0; j < 8; j++)
                {
                    tmp[i+j] = x[i+j];
                }
            }
            
            for (i = b2; i < a; i += b2)
            {
                for (j = 0; j < b2; j += 8)
                {
                    // TODO: AVX
                    for (k = 0; k < 8; k++)
                    {
                        tmp[j+k] = n_addmod(tmp[j+k], x[i+j+k], p);
                    }
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
void CRT(mp_hlimb_t * x, mp_hlimb_t * tmp, slong n, mp_hlimb_t p)
{
    slong a, b, b2, b8, n2, n8, i, j, nn;
    mp_hlimb_t half;

    nn = n;
    half = 1 + (p >> 1);
        
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
        for (i = 0; i < b8; i+= 8)
        {
            // TODO: AVX2
            for (j = 0; j < 8; j++)
            {
                x[a+i+j] = n_addmod(n_addmod(tmp[b+i+j], tmp[b+i+j], p), x[a+i+j], p);
            }

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
            //TODO: AVX2
            for (j = 0; j < 8; j++)
            {
                mp_hlimb_t u;
                u = n_addmod(tmp[i+j], tmp[b+i+j], p);
                u = n_submod(x[a+i+j], u, p);
                x[a+i+j] = (u >> 1) + (u & 1)*half;   // u/2 mod p
                x[i+j] = n_addmod(x[i+j], x[a+i+j], p);
            }
        }
        // less than 8 passes
        for (; i < b; i++)
        {
            mp_hlimb_t u;
            u = n_addmod(tmp[i], tmp[b+i], p);
            u = n_submod(x[a+i], u, p);
            x[a+i] = (u >> 1) + (u & 1)*half;   // u/2 mod p
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
                x[a+i] = (x[a+i] >> 1) + (x[a+i] & 1)*half;   // u/2 mod p
                x[i] = n_addmod(x[i], x[a+i], p);
            }
        }
        else
        {
            // less than 8 passes
            for (i = b; i < b8; i++)
            {
                x[a+i] = (x[a+i] >> 1) + (x[a+i] & 1)*half;   // u/2 mod p
                x[i] = n_addmod(x[i], x[a+i], p);
            }
            for (i = b8; i < n8; i+=8)
            {
                //TODO: AVX2
                for (j = 0; j < 8; j++)
                {
                    x[a+i+j] = (x[a+i+j] >> 1) + (x[a+i+j] & 1)*half;   // u/2 mod p
                    x[i+j] = n_addmod(x[i+j], x[a+i+j], p);
                }
            }
            // less than 8 passes
            for (i = n8; i < n; i++)
            {
                x[a+i] = (x[a+i] >> 1) + (x[a+i] & 1)*half;   // u/2 mod p
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
    slong i, N;
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
            _inv_fft_32_2(x_32, F->powers_inv_w[2][3], F->i_powers_inv_w[2][3], F->mod);
            break;
        case 3:
            _inv_fft_32_3(x_32, F->powers_inv_w[3]+2, F->i_powers_inv_w[3]+2, F->mod);
            break;
        default:
            _inv_fft_32_k(x_32, F->powers_inv_w[k]+2, F->i_powers_inv_w[k]+2, F->mod, k);
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
    slong i, nn, k, aa;
    mp_hlimb_t * wk, * wk2, * powers, * i_powers;
    mp_hlimb_t p;
    
    nmod_poly_fit_length(poly, N);
    poly->length = N;

    wk = (mp_hlimb_t *) aligned_alloc(32, 12*N > 32 ? 12*N : 32);
    for (i = 0; i < N; i++)
    {
        wk[i] = x[i];
        wk[i+N] = 0;
        wk[i+2*N] = 0;
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
                _inv_fft_32_2(wk, F->powers_inv_w[2][3], F->i_powers_inv_w[2][3], F->mod);
                break;
            case 3:
                _inv_fft_32_3(wk, F->powers_inv_w[3]+2, F->i_powers_inv_w[3]+2, F->mod);
                break;
            default:
                _inv_fft_32_k(wk, F->powers_inv_w[k]+2, F->i_powers_inv_w[k]+2, F->mod, k);
                break;
        }
        
        powers = F->powers_inv_w_over_2[k];
        i_powers = F->i_powers_inv_w_over_2[k];

        for (i = 0; i < aa; i++)
        {
            wk[i] = mul_mod_precon_32(wk[i], powers[i], p, i_powers[i]);
        }
        
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
    CRT(wk, wk + N, N, p);

    for (i = 0; i < N; i++)
    {
        poly->coeffs[i] = wk[i];
    }
    
    free(wk);
    _nmod_poly_normalise(poly);
}

#endif
