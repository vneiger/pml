#include "nmod_poly_extra.h"

#ifdef HAS_INT128

/*------------------------------------------------------------*/
/*  in-place size 2^k inverse butterfly, k >= 3               */
/*------------------------------------------------------------*/
static void _inv_fft_k(mp_ptr x, const mp_ptr powers_inv_w_in, const mp_ptr i_powers_inv_w_in,
                       const nmod_t mod, const ulong k)
{
    mp_limb_t w, iw, p, p2;
    ulong N, M;
    slong r, i;
    mp_ptr powers_inv_w, i_powers_inv_w, x0, x1, x2;

    N = 4;
    M = 1L << (k-2);    // number of blocks of size 4 = n/4
    powers_inv_w = powers_inv_w_in;
    i_powers_inv_w = i_powers_inv_w_in;
    w = powers_inv_w[1];
    iw = i_powers_inv_w[1];

    p = mod.n;
    p2 = p+p;
    
    // first two layers: no reduction
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
        v3 = mul_mod_precon_unreduced(v3, w, p, iw); // [0,2p)
        
        x2[0] = v0 + v2; // [0,4p)
        x2[1] = v1 + v3; // [0,4p)
        x2[2] = (v0 + p2) - v2; // [0,4p)
        x2[3] = (v1 + p2) - v3; // [0, 4p) 
    }
    
    powers_inv_w += N/2;
    i_powers_inv_w += N/2;
    N *= 2;
    M /= 2;
    
    // middle layers
    for (; M >= 2; N *= 2, M /= 2)
    {
        x0 = x;
        x1 = x + N/2;
        for (r = 0; r < M; r++, x0 += N, x1 += N)
        {
            for (i = 0; i < N/2; i+=2)
            {
                mp_limb_t u0, u1;
                mp_limb_signed_t z;
                
                z = x0[i] - p2; // (-2p,2p)
                u0 = (z < 0) ? z + p2 : z; // [0,2p)
                u1 = mul_mod_precon_unreduced(x1[i], powers_inv_w[i], p, i_powers_inv_w[i]);  // [0,2p)
                x0[i] = u0 + u1; // [0,4p)
                x1[i] = (u0 + p2) - u1; // [0,4p)

                z = x0[i+1] - p2; // (-2p,2p)
                u0 = (z < 0) ? z + p2 : z; // [0,2p)
                u1 = mul_mod_precon_unreduced(x1[i+1], powers_inv_w[i+1], p, i_powers_inv_w[i+1]);  // [0,2p)
                x0[i+1] = u0 + u1; // [0,4p)
                x1[i+1] = (u0 + p2) - u1; // [0,4p)
            }
        }
        powers_inv_w += N/2;
        i_powers_inv_w += N/2;
    }

    // last layer, with full reduction
    x0 = x;
    x1 = x + N/2;

    for (i = 0; i < N/2; i+=2)
    {
        mp_limb_t u0, u1, q;
        mp_limb_signed_t z;
        
        z = x0[i] - p2;
        u0 = (z < 0) ? z + p2 : z;
        u1 = mul_mod_precon_unreduced(x1[i], powers_inv_w[i], p, i_powers_inv_w[i]);  // [0,2p)
        
        q = u0 + u1; 
        q -= (q >= p2) ? p2 : 0;
        q -= (q >= p) ? p : 0;
        x0[i] = q;
        
        q = (u0 + p2) - u1;
        q -= (q >= p2) ? p2 : 0;
        q -= (q >= p) ? p : 0;
        x1[i] = q;
        
        z = x0[i+1] - p2;
        u0 = (z < 0) ? z + p2 : z; // [0,2p)
        u1 = mul_mod_precon_unreduced(x1[i+1], powers_inv_w[i+1], p, i_powers_inv_w[i+1]);  // [0,2p)
        
        q = u0 + u1; // [0,4p)
        q -= (q >= p2) ? p2 : 0;
        q -= (q >= p) ? p : 0;
        x0[i+1] = q;
        
        q = (u0 + p2) - u1; // [0,4p)
        q -= (q >= p2) ? p2 : 0;
        q -= (q >= p) ? p : 0;
        x1[i+1] = q;
    }
}



/*------------------------------------------------------------*/
/* inverse fft                                                */
/* given x[i] = poly(w^bitreverse(i,2^k)), returns poly       */
/*------------------------------------------------------------*/
void nmod_int128_fft_interpolate(nmod_poly_t poly, mp_srcptr x, const nmod_int128_fft_t F, const ulong k)
{
    slong i, N;
    
    N = 1L << k;
    nmod_poly_fit_length(poly, N);
    poly->length = N;

    for (i = 0; i < N; i++)
    {
        poly->coeffs[i] = x[i];
    }

    // special cases
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
        default:
            _inv_fft_k(poly->coeffs, F->powers_inv_w[k]+2, F->i_powers_inv_w[k]+2, F->mod, k);
    }
 
    for (i = 0; i < N; i++)
    {
        poly->coeffs[i] = mul_mod_precon(poly->coeffs[i], F->powers_inv_2[k], F->mod.n, F->i_powers_inv_2[k]);
    }

    _nmod_poly_normalise(poly);
}

#endif
