#include "nmod_poly_extra.h"

/*------------------------------------------------------------*/
/*  in-place size 2^k inverse butterfly, k >= 3               */
/*------------------------------------------------------------*/
static void _inv_fft_k(mp_ptr x, const mp_ptr powers_inv_w_in, const nmod_t mod, const ulong k)
{
    // N=size of block, M=number of blocks
    mp_ptr powers_inv_w;
    ulong N, M;
    slong r;
    mp_limb_t w;
 
    powers_inv_w = powers_inv_w_in;
    w = powers_inv_w[1];
    N = 4;
    M = 1L << (k-2);    // number of blocks of size 4 = N/4

    // first two layers
    mp_ptr x2 = x;
    for (r = 0; r < M; r++, x2 += 4)
    {
        mp_limb_t u0, u1, u2, u3, v0, v1, v2, v3, z0, z1, z2, z3;
        u0 = x2[0];
        u1 = x2[1];
        u2 = x2[2];
        u3 = x2[3];
      
        BUTTERFLY(v0, v1, u0, u1, mod);
        BUTTERFLY(v2, v3, u2, u3, mod);
        BUTTERFLY(z0, z1, v0, v2, mod);
        v3 = nmod_mul(v3, w, mod);
        BUTTERFLY(z2, z3, v1, v3, mod);

        x2[0] = z0;
        x2[1] = z2;
        x2[2] = z1;
        x2[3] = z3;
    }

    powers_inv_w += N/2;
    N *= 2;
    M /= 2;

    // all remaining layers
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

                /* u0 = x0[i+1]; */
                /* u1 = x1[i+1]; */
                /* u1 = nmod_mul(u1, powers_inv_w[i+1], mod); */
                /* BUTTERFLY(v0, v1, u0, u1, mod); */
                /* x0[i+1] = v0; */
                /* x1[i+1] = v1; */
            }
        }
        powers_inv_w += N/2;
    }
}


/*------------------------------------------------------------*/
/* inverse fft                                                */
/* given x[i] = poly(w^bitreverse(i,2^k)), returns poly       */
/*------------------------------------------------------------*/
void nmod_plain_fft_interpolate(nmod_poly_t poly, mp_srcptr x, const nmod_plain_fft_t F, const ulong k)
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
            _inv_fft_k(poly->coeffs, F->powers_inv_w[k]+2, F->mod, k);
            break;
    }
 
    for (i = 0; i < N; i++)
    {
        poly->coeffs[i] = nmod_mul(poly->coeffs[i], F->powers_inv_2[k], F->mod);
    }

    _nmod_poly_normalise(poly);
}
