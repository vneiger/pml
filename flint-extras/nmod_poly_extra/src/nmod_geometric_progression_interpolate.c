#include "nmod_poly_extra.h"

/*------------------------------------------------------------*/
/* in: coeffs must have size at least d                       */
/* out: interpolating polynomial C s.t. C(q^i) = v[i], i<d    */
/*------------------------------------------------------------*/
void nmod_geometric_progression_interpolate(nmod_poly_t poly, mp_srcptr v, const nmod_geometric_progression_t G)
{
    slong i, N;
    nmod_poly_t f, h;
    nmod_t mod;

    N = G->d;
    nmod_poly_fit_length(poly, N);
    nmod_poly_zero(poly);
    
    if (N == 1)
    {
        nmod_poly_set_coeff_ui(poly, 0, v[0]);
        return;
    }

    mod = G->mod;
    nmod_poly_init2(f, mod.n, N);
    nmod_poly_init2(h, mod.n, N);
    
    for(i = 0; i < N; i++)
    {
        nmod_poly_set_coeff_ui(f, i, nmod_mul(v[i], G->w[i], mod));
    }

    nmod_poly_mullow(h, f, G->g1, N);

    for (i = 0; i < N; i++)
    {
        nmod_poly_set_coeff_ui(f, N - 1 - i, nmod_mul(nmod_poly_get_coeff_ui(h, i), G->y[i], mod));
    }

    nmod_poly_mullow(h, f, G->g2, N);

    for (i = 0; i < N; i++)
    {
        nmod_poly_set_coeff_ui(poly, i, nmod_mul(nmod_poly_get_coeff_ui(h, N - 1 - i), G->z[i], mod));
    }

    nmod_poly_clear(f);
    nmod_poly_clear(h);
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
