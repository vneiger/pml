#include <flint/nmod_mat.h> // for nmod_mat_randtest()
#include "nmod_mat_extra.h" // for nmod_mat_rand()
#include "nmod_mat_poly.h"

void nmod_mat_poly_randtest(nmod_mat_poly_t matp,
                            flint_rand_t state, 
                            slong len)
{
    nmod_mat_poly_fit_length(matp, len);
    _nmod_mat_poly_set_length(matp, len);
    nmod_mat_t cmat;
    for (slong i = 0; i < len; ++i)
    {
        nmod_mat_poly_coeff_attach(cmat, matp, i);
        nmod_mat_randtest(cmat, state);
    }
    _nmod_mat_poly_normalise(matp);
}

void nmod_mat_poly_rand(nmod_mat_poly_t matp,
                        flint_rand_t state, 
                        slong len)
{
    nmod_mat_poly_fit_length(matp, len);
    _nmod_mat_poly_set_length(matp, len);
    nmod_mat_t cmat;
    for (slong i = 0; i < len; ++i)
    {
        nmod_mat_poly_coeff_attach(cmat, matp, i);
        nmod_mat_rand(cmat, state);
    }
    _nmod_mat_poly_normalise(matp);
}
