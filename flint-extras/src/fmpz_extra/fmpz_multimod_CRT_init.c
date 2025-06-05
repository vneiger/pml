#include <stdlib.h>
#include <flint/flint.h>
#include <flint/ulong_extras.h>
#include <flint/fmpz.h>

#include "fmpz_extra.h"

/* ------------------------------------------------------------ */
/*                                                              */
/* ------------------------------------------------------------ */
void fmpz_multimod_CRT_init(fmpz_multimod_CRT_t mmod, nn_srcptr primes, ulong num_primes)
{
    double quo;
    ulong size, i, cpt;
    int res;
    nn_ptr cofactors, one;
    fmpz * partial_products;
    fmpz_t linearized_product;
    
    if (num_primes < MULTIMOD_CRT_LEAF_SIZE)
        size = num_primes;
    else
    {
        quo = ceil(num_primes / MULTIMOD_CRT_LEAF_SIZE);
        size = num_primes / pow(2, ceil(log(quo)/log(2)));
    }

    mmod->num_primes = num_primes;
    mmod->size_leaves = size;
    mmod->num_leaves = (num_primes + size - 1) / size;
    mmod->leaves_mod = (fmpz_multimod_naive_t *) malloc(mmod->num_leaves * sizeof(fmpz_multimod_naive_t));
    mmod->products_leaves = (fmpz *) malloc(mmod->num_leaves * sizeof(fmpz));
    mmod->leaves_CRT = (fmpz_CRT_naive_t *) malloc(mmod->num_leaves * sizeof(fmpz_CRT_naive_t));
    mmod->inverse_cofactors = _nmod_vec_init(num_primes);
    fmpz_init(mmod->product_primes);
    
    cofactors = _nmod_vec_init(num_primes);
    one = _nmod_vec_init(num_primes);
    partial_products = (fmpz *) malloc(mmod->num_leaves * sizeof(fmpz));
    fmpz_init(linearized_product);

    for (i = 0; i < num_primes; i++)
        one[i] = 1;

    fmpz_set_ui(mmod->product_primes, 1);
    for (i = 0; i < mmod->num_leaves; i++)
    {
        fmpz_multimod_naive_init(mmod->leaves_mod[i], primes + i * size, FLINT_MIN(num_primes - i * size, size));
        fmpz_init_set(mmod->products_leaves + i, mmod->leaves_mod[i]->prod);
        fmpz_CRT_naive_init(mmod->leaves_CRT[i], primes + i * size, FLINT_MIN(num_primes - i * size, size));
        fmpz_init(partial_products + i);
        // partial_products[i] = result of CRT_combine with leaf values = 1
        // if e.g. moduli[i] are p0,p1,p2, this is p1*p2+p0*p2+p0*p1
        fmpz_CRT_naive_combine(partial_products + i, one + i * size, mmod->leaves_CRT[i]);
        // naive product
        fmpz_mul(mmod->product_primes, mmod->product_primes, mmod->leaves_CRT[i]->prod);
    }

    fmpz_multi_mod_init(mmod->top_mod);
    res = fmpz_multi_mod_precompute(mmod->top_mod, mmod->products_leaves, mmod->num_leaves);
    if (!res)
        flint_throw(FLINT_ERROR, "fmpz_multimod_CRT_init: moduli are not nonzero");

    fmpz_multi_CRT_init(mmod->top_CRT);
    res = fmpz_multi_CRT_precompute(mmod->top_CRT, mmod->products_leaves, mmod->num_leaves);
    if (!res)
        flint_throw(FLINT_ERROR, "fmpz_multimod_CRT_init: moduli are not coprime");

    // combine all partial products 
    fmpz_multi_CRT_combine(linearized_product, mmod->top_CRT, partial_products);

    // reduce it modulo all moduli
    fmpz_multimod_CRT_reduce(cofactors, linearized_product, mmod);
    // modulo-inverts each such remainder
    // todo: Montgomery inversion
    cpt = 0;
    for (i = 0; i < mmod->num_leaves; i++)
    {
        ulong j, nb;
        nb = mmod->leaves_CRT[i]->num_primes;
        for (j = 0; j < nb; j++)
        {
            mmod->inverse_cofactors[cpt] = nmod_inv(cofactors[cpt], mmod->leaves_CRT[i]->mod[j]);
            cpt++;
        }
    }
    
    fmpz_clear(linearized_product);
    for (i = 0; i < mmod->num_leaves; i++)
        fmpz_clear(partial_products + i);
    free(partial_products);
    _nmod_vec_clear(cofactors);
    _nmod_vec_clear(one);
}
