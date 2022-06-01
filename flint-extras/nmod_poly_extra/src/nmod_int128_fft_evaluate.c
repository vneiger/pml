#include "nmod_poly_extra.h"

#ifdef HAS_INT128

/*------------------------------------------------------------*/
/* in-place FFT in size 2^k, k >= 3                           */
/* powers_w_in are the roots of 1 needed                      */
/*------------------------------------------------------------*/
static void _fft_k(mp_ptr x, const mp_ptr powers_w_in, const mp_ptr i_powers_w_in, const nmod_t mod, const ulong k)
{
    // N=size of block, M=number of blocks
    ulong N, M;
    mp_limb_t w, iw, p, p2;
    slong r, i;
    mp_ptr powers_w, i_powers_w, x0, x1;

    powers_w = powers_w_in;
    i_powers_w = i_powers_w_in;
    
    N = 1L<<k;
    M = 1;

    p = mod.n;
    p2 = p+p;
    
    // the first layer is different: no need to reduce the sums
    x0 = x;
    x1 = x + N/2;
    for (r = 0; r < M; r++, x0 += N, x1 += N)
    {
        for (i = 0; i < N/2; i+=2)
        {
            mp_limb_t u0, u1, t0, t1;

            u0 = x0[i];
            u1 = x1[i];
            t0 = u0 + u1;
            t1 = (u0 + p) - u1;
            
            x0[i] = t0;
            x1[i] = mul_mod_precon_unreduced(t1, powers_w[i], p, i_powers_w[i]);

            u0 = x0[i+1];
            u1 = x1[i+1];

            t0 = u0 + u1;
            t1 = (u0  + p) - u1;
            t1 = mul_mod_precon_unreduced(t1, powers_w[i+1], p, i_powers_w[i+1]);
            x0[i+1] = t0;
            x1[i+1] = t1;
        }
    }
    powers_w += N/2;
    i_powers_w += N/2;
    N /= 2;
    M *= 2;
    
    for (; N > 4; N /= 2, M *= 2)
    {
        x0 = x;
        x1 = x + N/2;

        for (r = 0; r < M; r++, x0 += N, x1 += N)
        {
            slong i;
            for (i = 0; i < N/2; i+=2)
            {
                mp_limb_t u0, u1, t1;
                mp_limb_signed_t t0;
                
                u0 = x0[i];
                u1 = x1[i];
                t0 = u0 + u1 - p2;
                x0[i] = (t0 < 0) ? t0+p2 : t0;
                t1 = (u0 + p2) - u1;
                x1[i] = mul_mod_precon_unreduced(t1, powers_w[i], p, i_powers_w[i]);

                u0 = x0[i+1];
                u1 = x1[i+1];
                t0 = u0 + u1 - p2;
                x0[i+1] = (t0 < 0) ? t0+p2 : t0;
                t1 = (u0 + p2) - u1;
                x1[i+1] = mul_mod_precon_unreduced(t1, powers_w[i+1], p, i_powers_w[i+1]);
            }
        }
        powers_w += N/2;
        i_powers_w += N/2;
    }
    
    // last two layers
    w = powers_w[1];
    iw = i_powers_w[1];

    for (r = 0; r < M; r++, x += 4)
    {
        mp_limb_t u0, u1, u2, u3, v0, v1, v2, v3, z0, z1, z2, z3;
        
        u0 = x[0];
        u1 = x[1];
        u2 = x[2];
        u3 = x[3];

        v0 = u0 + u2;
        v0 -= (v0 >= p2) ? (p2) : 0;
        
        v2 = u0 - u2;
        v2 += ((mp_limb_signed_t) v2 < 0) ? (p2) : 0;
        
        v1 = u1 + u3;
        v1 -= (v1 >= p2) ? (p2) : 0;

        v3 = u1 - u3 + p2;
        v3 = mul_mod_precon_unreduced(v3, w, p, iw);
        
        z0 = v0 + v1;
        z0 -= (z0 >= p2) ? (p2) : 0;
        z0 -= (z0 >= p) ? (p) : 0;

        z1 = v0 - v1;
        z1 += ((mp_limb_signed_t) z1 < 0) ? (p2) : 0;
        z1 -= (z1 >= p) ? (p) : 0;
        
        z2 = v2 + v3;
        z2 -= (z2 >= p2) ? (p2) : 0;
        z2 -= (z2 >= p) ? (p) : 0;

        z3 = v2 - v3;
        z3 += ((mp_limb_signed_t) z3 < 0) ? (p2) : 0;
        z3 -= (z3 >= p) ? (p) : 0;
        
        x[0] = z0;
        x[1] = z1;
        x[2] = z2;
        x[3] = z3;
    }

}


/*------------------------------------------------------------*/
/* fft evaluation                                             */
/* returns x[i] = poly(w^bitreverse_n(i)), n=2^k              */
/* x must have length >= n                                    */
/*------------------------------------------------------------*/
void nmod_int128_fft_evaluate(mp_ptr x, const nmod_poly_t poly, const nmod_int128_fft_t F, const ulong k)
{
    slong i, N;
    
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
        default:
            _fft_k(x, F->powers_w[k], F->i_powers_w[k], F->mod, k);
            break;
    }
}

#endif
