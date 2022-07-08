#include "nmod_poly_extra.h"

/*------------------------------------------------------------*/
/* in-place FFT in size 2^k, k >= 2                           */
/* powers_w_in are the roots of 1 needed                      */
/*------------------------------------------------------------*/
static void _fft_k(mp_ptr x, const mp_ptr powers_w_in, const nmod_t mod, const ulong k)
{
    // N=size of block, M=number of blocks
    mp_ptr powers_w;
    ulong N, M;

    powers_w = powers_w_in;
    for (N = 1L<<k, M = 1; N > 1; N /= 2, M *= 2)
    {
        mp_ptr x0, x1;
        slong r;

        x0 = x;
        x1 = x + N/2;
        for (r = 0; r < M; r++, x0 += N, x1 += N)
        {
            slong i;
            for (i = 0; i < N/2; i++)
            {
                mp_limb_t u0, u1, t0, t1;
                u0 = x0[i];
                u1 = x1[i];
                BUTTERFLY(t0, t1, u0, u1, mod);
                t1 = nmod_mul(t1, powers_w[i], mod);
                x0[i] = t0;
                x1[i] = t1;
            }
        }
        powers_w += N/2;
    }
}

/*------------------------------------------------------------*/
/* fft evaluation                                             */
/* returns x[i] = poly(w^bitreverse_n(i)), n=2^k              */
/* x must have length >= n                                    */
/*------------------------------------------------------------*/
void nmod_fft_evaluate(mp_ptr x, const nmod_poly_t poly, const nmod_fft_t F, const ulong k)
{
    slong i, N;
    
    N = 1L << k;
    for (i = 0; i < N; i++)
    {
        x[i] = nmod_poly_get_coeff_ui(poly, i);
    }
    
    _fft_k(x, F->powers_w[k], F->mod, k);
}



/*------------------------------------------------------------*/
/* tft evaluation                                             */
/* returns a cyclotomic tft evaluation of x at n points       */
/* x must have length >= N                                    */
/*------------------------------------------------------------*/
void nmod_tft_evaluate(mp_ptr x, const nmod_poly_t poly, const nmod_fft_t F, const ulong N)
{
    slong i, j, k, a, ell, b, t, b2, lambda;
    mp_limb_t p;
    mp_ptr x2, x2_bak, powers_rho;
    mp_limb_t u, v;

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
    x2 = _nmod_vec_init(2*N);
    x2_bak = x2;
    
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
        switch(b)
        {
            case 0:
                for (i = 0; i < a; i++)
                {
                    x2[i] = nmod_mul(n_submod(x2[i], x2[i+a], p), powers_rho[i], F->mod);
                }
                break;

            default:
                // i = 0 to 2^(k-ell-1) = a/(2b), j = 0 to 2b
                // i = 0 : special case
                for (; t < b2; t++)
                {
                    u = x2[t];
                    v = x2[a+t];
                    x2[t] = nmod_mul(n_submod(u, v, p), powers_rho[t], F->mod);
                    x2[a+t] = n_addmod(u, v, p);
                }
                
                lambda = (1L << k) >> (ell+1);
                for (i = 1; i < lambda; i++)
                {
                    for (j = 0; j < b2; j++)
                    {
                        u = x2[t];
                        v = x2[a+t];
                        x2[t] = nmod_mul(n_submod(u,  v, p), powers_rho[t], F->mod);
                        x2[a+j] = n_addmod(n_addmod(u, v, p), x2[a+j], p);
                        t++;
                    }
                }
        }

        _fft_k(x2, F->powers_w[k], F->mod, k);
        
        x2 += a;
    }
    while (b != 0);
    
    x2 = x2_bak;
    for (i = 0; i < N; i++)
    {
        x[i] = x2[i];
    }

    _nmod_vec_clear(x2);
}

