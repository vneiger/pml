#include "nmod_poly_extra.h"
#include "nmod_vec_extra.h"

#ifdef HAS_AVX2
#include <immintrin.h>

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* using AVX2 to do inverse FFT nmod_poly with p < 2^32       */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

static void _inv_fft_3(mp_hlimb_t * x, mp_hlimb_t * powers_inv_w, mp_hlimb_t * i_powers_inv_w, const nmod_t mod)
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
static void _inv_fft_k(mp_hlimb_t * x, mp_hlimb_t * powers_inv_w, mp_hlimb_t * i_powers_inv_w,
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


/*------------------------------------------------------------*/
/* inverse fft                                                */
/* given x[i] = poly(w^bitreverse(i,2^k)), returns poly       */
/*------------------------------------------------------------*/
void nmod_avx2_32_fft_interpolate(nmod_poly_t poly, mp_srcptr x, const nmod_avx2_32_fft_t F, const ulong k)
{
    slong i, N;
    mp_hlimb_t * x_32;

    N = 1L << k;
    nmod_poly_fit_length(poly, N);
    poly->length = N;

    if (k <= 2)
    {
        for (i = 0; i < N; i++)
        {
            poly->coeffs[i] = x[i];
        }
        switch(k)
        {
            case 0:
                break;
            case 1:
                _inv_fft_1(poly->coeffs, F->mod);
                break;
            case 2:
                _inv_fft_2(poly->coeffs, F->powers_inv_w[2][3], F->mod);
                break;
        }
        for (i = 0; i < N; i++)
        {
            poly->coeffs[i] = nmod_mul(poly->coeffs[i], F->powers_inv_2[k], F->mod);
        }
    }
    else
    {
        x_32 = (mp_hlimb_t *) aligned_alloc(32, 4 * N);
        for (i = 0; i < N; i++)
        {
            x_32[i] = (mp_hlimb_t) x[i];
        }

        if (k == 3)
        {
            _inv_fft_3(x_32, F->powers_inv_w[3]+2, F->i_powers_inv_w[3]+2, F->mod);
            
        }
        else
        {
            _inv_fft_k(x_32, F->powers_inv_w[k]+2, F->i_powers_inv_w[k]+2, F->mod, k);
        }
        _nmod_vec_avx2_32_scalar_mul(x_32, x_32, N, F->powers_inv_2[k], F->i_powers_inv_2[k], F->mod);
        
        for (i = 0; i < N; i++)
        {
            poly->coeffs[i] = (mp_limb_t) x_32[i];
        }
        free(x_32);
    }
    _nmod_poly_normalise(poly);
}


#endif
