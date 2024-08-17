#include <flint/nmod.h>
#include "nmod_poly_extra.h"

/*------------------------------------------------------------*/
/* in: polynomial to evaluate, deg(poly) < d                  */
/* out: v[i] = poly(q^i), i = 0 .. d-1                        */
/*------------------------------------------------------------*/
void nmod_geometric_progression_evaluate(nn_ptr v, const nmod_poly_t poly, const nmod_geometric_progression_t G)
{
    nmod_poly_t a, b;
    slong i, d, dp;

    dp = nmod_poly_degree(poly);
    d = G->d;

    if (dp == -1)
    {
        for (i = 0; i < d; i++)
        {
            v[i] = 0;
        }
        return;
    }

    nmod_poly_init2(a, G->mod.n, d);
    nmod_poly_init(b, G->mod.n);

    for (i = 0; i < d; i++)
    {
        nmod_poly_set_coeff_ui(a, d-1-i, nmod_mul(G->x[i], nmod_poly_get_coeff_ui(poly, i), G->mod));
    }

    nmod_poly_mul(b, a, G->f);

    for (i = 0; i < d; i++)
    {
        v[i] = nmod_mul(G->x[i], nmod_poly_get_coeff_ui(b, i+d-1), G->mod);
    }

    nmod_poly_clear(b);
    nmod_poly_clear(a);
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
