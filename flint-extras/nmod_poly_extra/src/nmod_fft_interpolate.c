#include "nmod_poly_extra.h"

/*------------------------------------------------------------*/
/*  in-place size 2^k inverse butterfly                       */
/*------------------------------------------------------------*/
static void _inv_fft_k(mp_ptr x, const mp_ptr powers_inv_w_in, const nmod_t mod, const ulong k)
{
    // N=size of block, M=number of blocks
    mp_ptr powers_inv_w;
    slong N, M;
    slong r;
 
    powers_inv_w = powers_inv_w_in;
    N = 2;
    M = (1L << k) >> 1;    // number of blocks of size 4 = N/4

    for (; M >= 1; N *= 2, M /= 2)
    {
        mp_ptr x0, x1;
        x0 = x;
        x1 = x + N/2;

        for (r = 0; r < M; r++, x0 += N, x1 += N)
        {
            slong i;
            for (i = 0; i < N/2; i+=1)
            {
                mp_limb_t u0, u1, v0, v1;
                u0 = x0[i];
                u1 = x1[i];
                u1 = nmod_mul(u1, powers_inv_w[i], mod);
                BUTTERFLY(v0, v1, u0, u1, mod);
                x0[i] = v0;
                x1[i] = v1;
            }
        }
        powers_inv_w += N/2;
    }
}


/*-----------------------------------------------------------*/
/* applies the CRT map to x in place                         */
/* x has length n                                            */
/* input is in [0,p), output is in [0,p)                     */
/* tmp is a temporary workspace, of length at least 2n       */
/*-----------------------------------------------------------*/
static inline
void CRT(mp_ptr x, mp_ptr tmp, slong n, mp_limb_t p)
{
    slong a, b, b2, n2, i, j, nn;
    mp_limb_t half;

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
        

        for (i = 0; i < b2; i++)
        {
            tmp[i] = x[i];
        }
        
        while (i < a)
        {
            for (j = 0; j < b2; j++, i++)
            {
                tmp[j] = n_addmod(tmp[j], x[i], p);
            }
        }
        
        for(i = 0; i < b; i++)
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
        
        for (i = 0; i < b; i++)
        {
            mp_limb_t u;
            u = n_addmod(tmp[i], tmp[b+i], p);
            u = n_submod(x[a+i], u, p);
            x[a+i] = (u >> 1) + (u & 1)*half;   // u/2 mod p
            x[i] = n_addmod(x[i], x[a+i], p);
        }
        
        for (i = b; i < n; i++)
        {
            x[a+i] = (x[a+i] >> 1) + (x[a+i] & 1)*half;   // u/2 mod p
            x[i] = n_addmod(x[i], x[a+i], p);
        }
        
        n = n + a;
    }
}


/*------------------------------------------------------------*/
/* inverse fft                                                */
/* given x[i] = poly(w^bitreverse(i,2^k)), returns poly       */
/*------------------------------------------------------------*/
void nmod_fft_interpolate(nmod_poly_t poly, mp_srcptr x, const nmod_fft_t F, const ulong k)
{
    slong i, N;
    
    N = 1L << k;
    nmod_poly_fit_length(poly, N);
    poly->length = N;

    for (i = 0; i < N; i++)
    {
        poly->coeffs[i] = x[i];
    }

    _inv_fft_k(poly->coeffs, F->powers_inv_w[k] + 1, F->mod, k);
            
    for (i = 0; i < N; i++)
    {
        poly->coeffs[i] = nmod_mul(poly->coeffs[i], F->powers_inv_2[k], F->mod);
    }

    _nmod_poly_normalise(poly);
}



/*------------------------------------------------------------*/
/* inverse tft                                                */
/*------------------------------------------------------------*/
void nmod_tft_interpolate(nmod_poly_t poly, mp_srcptr x, const nmod_fft_t F, const ulong N)
{
    slong i, nn, k, aa;
    mp_ptr wk, wk2, powers;
        
    nmod_poly_fit_length(poly, N);
    poly->length = N;

    wk = _nmod_vec_init(3*N);

    for (i = 0; i < N; i++)
    {
        wk[i] = x[i];
        wk[i+N] = 0;
        wk[i+2*N] = 0;
    }

    wk2 = wk;
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
        _inv_fft_k(wk, F->powers_inv_w[k] + 1, F->mod, k);
        powers = F->powers_inv_w_over_2[k];

        for (i = 0; i < aa; i++)
        {
            wk[i] = nmod_mul(wk[i], powers[i], F->mod);
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
    CRT(wk, wk + N, N, F->mod.n);

    for (i = 0; i < N; i++)
    {
        poly->coeffs[i] = wk[i];
    }
    
    _nmod_vec_clear(wk);
    _nmod_poly_normalise(poly);
}
