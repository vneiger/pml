#include "nmod_mat_poly.h"

void nmod_mat_poly_print(const nmod_mat_poly_t matp)
{
    nmod_mat_t cmat;
	for (slong i = 0; i < matp->length; i++)
    {
        nmod_mat_poly_coeff_attach(cmat, matp, i);
		nmod_mat_print(cmat);
    }
}

void nmod_mat_poly_print_pretty(const nmod_mat_poly_t matp)
{
    printf("r: %ld, c: %ld, len: %ld, mod: %ld\n",
           matp->r, matp->c, matp->length, matp->mod.n);
    nmod_mat_t cmat;
	for (slong i = 0; i < matp->length; i++)
    {
        printf("coeff %ld:\n", i);
        nmod_mat_poly_coeff_attach(cmat, matp, i);
        nmod_mat_print_pretty(cmat);
    }
}
