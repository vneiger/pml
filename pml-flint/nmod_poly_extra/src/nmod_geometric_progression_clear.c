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
