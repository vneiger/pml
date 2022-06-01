#include "nmod_poly_extra.h"

/*------------------------------------------------------------*/
/* in-place FFT in size 2^k, k >= 3                           */
/* powers_w_in are the roots of 1 needed                      */
/*------------------------------------------------------------*/
static void _fft_k(mp_ptr x, const mp_ptr powers_w_in, const nmod_t mod, const ulong k)
{
    // N=size of block, M=number of blocks
    mp_ptr powers_w;
    ulong N, M;
    mp_limb_t w;
    slong r;

    powers_w = powers_w_in;
    for (N = 1L<<k, M = 1; N > 4; N /= 2, M *= 2)
    {
        mp_ptr x0, x1;
        slong r;

        x0 = x;
        x1 = x + N/2;

        for (r = 0; r < M; r++, x0 += N, x1 += N)
        {
            slong i;
            i = N/2 - 2;
            do
            {
                mp_limb_t u0, u1, t0, t1;
                
                u0 = x0[i+1];
                u1 = x1[i+1];
                
                t0 = u0 + u1;
                t1 = u0 - u1;
                CORRECT_0_2P(t0, mod.n);
                CORRECT_MINUSP_P(t1, mod.n);
                t1 = nmod_mul(t1, powers_w[i+1], mod);
                
                x0[i+1] = t0;
                x1[i+1] = t1;
                //
                u0 = x0[i];
                u1 = x1[i];

                t0 = u0 + u1;
                t1 = u0 - u1;
                CORRECT_0_2P(t0, mod.n);
                CORRECT_MINUSP_P(t1, mod.n);
                t1 = nmod_mul(t1, powers_w[i], mod);
                
                x0[i] = t0;
                x1[i] = t1;
                
                i -= 2;
            } while (i >= 0);
        }
        powers_w += N/2;
    }
    
    // last two layers
    w = powers_w[1];

    for (r = 0; r < M; r++, x += 4)
    {
        mp_limb_t u0, u1, u2, u3, v0, v1, v2, v3, z0, z1, z2, z3;
        
        u0 = x[0];
        u1 = x[1];
        u2 = x[2];
        u3 = x[3];
        
        v0 = u0 + u2;
        CORRECT_0_2P(v0, mod.n);
        v2 = u0 - u2;
        CORRECT_MINUSP_P(v2, mod.n);
        
        v1 = u1 + u3;
        CORRECT_0_2P(v1, mod.n);
        v3 = u1 - u3;
        CORRECT_MINUSP_P(v3, mod.n);
        v3 = nmod_mul(v3, w, mod);

        z0 = v0 + v1;
        CORRECT_0_2P(z0, mod.n);
        z1 = v0 - v1;
        CORRECT_MINUSP_P(z1, mod.n);

        z2 = v2 + v3;
        CORRECT_0_2P(z2, mod.n);
        z3 = v2 - v3;
        CORRECT_MINUSP_P(z3, mod.n);
        
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
void nmod_plain_fft_evaluate(mp_ptr x, const nmod_poly_t poly, const nmod_plain_fft_t F, const ulong k)
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
            _fft_k(x, F->powers_w[k], F->mod, k);
            break;
    }
} 
