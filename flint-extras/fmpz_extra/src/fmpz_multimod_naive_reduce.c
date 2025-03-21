#include <stdlib.h>
#include <flint/flint.h>
#include <flint/fmpz.h>

#include "fmpz_extra.h"
#include "nmod_vec_extra.h"


/* ------------------------------------------------------------ */
/* computes a mod mod                                           */
/* assumes that powers_of_two are the powers of 2 mod mod       */
/* assumes that number of powers of 2 <= num_limbs(a)           */
/* ------------------------------------------------------------ */
static 
ulong _fmpz_reduce(const fmpz_t a, nn_srcptr powers_of_two, const nmod_t mod)
{
    ulong res;
    if (!COEFF_IS_MPZ(*a))
    {
	fmpz s_a = *a;
	if (s_a < 0)
	{
	    NMOD_RED(res, -s_a, mod);
	    return nmod_neg(res, mod);
	}
	else
	{
	    NMOD_RED(res, s_a, mod);
	    return res;
	}
    }
    else
    {
        __mpz_struct *a_ptr;
        nn_ptr a_coeffs;    
        slong slen;
        a_ptr = COEFF_TO_PTR(*a);
	a_coeffs = a_ptr->_mp_d;
        
        slen = a_ptr->_mp_size;
        if (slen < 0)
            return nmod_neg(nmod_vec_dot_product_unbalanced(powers_of_two, a_coeffs, -slen, FLINT_BIT_COUNT(mod.n), FLINT_BITS, mod), mod);
	else
            return nmod_vec_dot_product_unbalanced(powers_of_two, a_coeffs, slen, FLINT_BIT_COUNT(mod.n), FLINT_BITS, mod);
    }
}



/* ------------------------------------------------------------ */
/* computes A mod mmod[i] for all i                             */
/* assumes all moduli are less than 2^30                        */
/* ------------------------------------------------------------ */
static 
void _fmpz_reduce_small_moduli(nn_ptr out, const fmpz_t a, const fmpz_multimod_naive_t mmod)
{
    if (!COEFF_IS_MPZ(*a))
    {
        ulong res;
        ulong i, num_primes;
	fmpz s_a;

        s_a = *a;
        num_primes = mmod->num_primes;
	if (s_a < 0)
	{
            for (i = 0; i < num_primes; i++)
            {
                NMOD_RED(res, -s_a, mmod->mod[i]);
                out[i] = nmod_neg(res, mmod->mod[i]);
            }
	}
	else
            for (i = 0; i < num_primes; i++)
                NMOD_RED(out[i], s_a, mmod->mod[i]);

        return;
    }
    else
    {
        slong slen;
        ulong i, j, num_limbs, len;
        nn_ptr slice_A, a_coeffs;
        __mpz_struct *a_ptr;
        
        a_ptr = COEFF_TO_PTR(*a);
	a_coeffs = a_ptr->_mp_d;
        slen = a_ptr->_mp_size;
        
        if (slen < 0)
            num_limbs = -slen;
        else
            num_limbs = slen;

        len = 3 * num_limbs;
        len = ((len + 3) >> 2) << 2; // must be a multiple of 4
        slice_A = (nn_ptr) aligned_alloc(32, len * sizeof(ulong));


        j = 0;
        for (i = 0; i < num_limbs; i++)
        {
            ulong limb;

            limb = a_coeffs[i];

            slice_A[j] = limb & 16777215;
            slice_A[j + 1] = (limb >> 24) & 16777215;
            slice_A[j + 2] = limb >> 48;
            j += 3;
        }

        for (; j < len; j++)
            slice_A[j] = 0;

        const ulong pow2 = UWORD(1) << DOT_SPLIT_BITS;
        for (i = 0; i < mmod->num_primes; i++)
        {
            ulong pow2_red;
            NMOD_RED(pow2_red, pow2, mmod->mod[i]);
            ulong dot = _nmod_vec_dot2_split(mmod->powers_of_two[i], slice_A, len, mmod->mod[i], pow2_red);
            if (slen < 0)
                out[i] = nmod_neg(dot, mmod->mod[i]);
            else
                out[i] = dot;
        }
    }
}


    
/* ------------------------------------------------------------ */
/* computes A mod mmod[i] for all i                             */
/* ------------------------------------------------------------ */
void fmpz_multimod_naive_reduce(nn_ptr out, const fmpz_t A, const fmpz_multimod_naive_t mmod)
{
    ulong i;

    if (fmpz_cmpabs(mmod->prod, A) <= 0)
    {
        fmpz_t B;
        fmpz_init(B);
        fmpz_mod(B, A, mmod->prod);
        fmpz_multimod_naive_reduce(out, B, mmod);
        fmpz_clear(B);
        return;
    }

    // small moduli: 
    if (mmod->small_moduli == 1)
        _fmpz_reduce_small_moduli(out, A, mmod);
    else
        for (i = 0; i < mmod->num_primes; i++)
            out[i] = _fmpz_reduce(A, mmod->powers_of_two[i], mmod->mod[i]);
}

