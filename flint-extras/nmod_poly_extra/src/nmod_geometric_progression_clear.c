#include <flint/nmod_vec.h>
#include "nmod_poly_extra.h"

/*------------------------------------------------------------*/
/* frees all memory attached to G                             */
/*------------------------------------------------------------*/
void nmod_geometric_progression_clear(nmod_geometric_progression_t G)
{
    nmod_poly_clear(G->f);
    nmod_poly_clear(G->g2);
    nmod_poly_clear(G->g1);
    _nmod_vec_clear(G->x);
    _nmod_vec_clear(G->z);
    _nmod_vec_clear(G->y);
    _nmod_vec_clear(G->w);
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
