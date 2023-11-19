#include <flint/flint.h>
#include <flint/fmpz.h>

#include "fmpz_extra.h"


/* ------------------------------------------------------------ */
/* number of limbs needed for a dot product of length len       */
/* with one input reduced mod mod, the other one arbitrary      */
/* ------------------------------------------------------------ */
static inline 
int _nmod_vec_dot_bound_limbs_unbalanced(ulong len, nmod_t mod)
{
    mp_limb_t t2, t1, t0, u1, u0;
    
    umul_ppmm(t1, t0, mod.n - 1, -1);
    umul_ppmm(t2, t1, t1, len);
    umul_ppmm(u1, u0, t0, len);
    add_sssaaaaaa(t2, t1, t0,  t2, t1, UWORD(0),  UWORD(0), u1, u0);

    if (t2 != 0) return 3;
    if (t1 != 0) return 2;
    return (t0 != 0);
}

/* ------------------------------------------------------------ */
/* computes a mod mod                                           */
/* assumes that bit_length(mod) = bit_len                       */
/* assumes that powers_of_two are the powers of 2 mod mod       */
/* ------------------------------------------------------------ */
static inline 
mp_limb_t _fmpz_reduce(const fmpz_t a, mp_srcptr powers_of_two, const nmod_t mod)
{
    mp_limb_t res;
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
        mp_ptr a_coeffs;    
        slong slen;
        a_ptr = COEFF_TO_PTR(*a);
	a_coeffs = a_ptr->_mp_d;
        
        slen = a_ptr->_mp_size;
        if (slen < 0)
            return nmod_neg(_nmod_vec_dot(powers_of_two, a_coeffs, -slen, mod,
                                          _nmod_vec_dot_bound_limbs_unbalanced(-slen, mod)), mod);
	else
	    return _nmod_vec_dot(powers_of_two, a_coeffs, slen, mod,
                                 _nmod_vec_dot_bound_limbs_unbalanced(slen, mod));
    }
}

/* ------------------------------------------------------------ */
/* computes A mod mmod[i] for all i                             */
/* ------------------------------------------------------------ */
void fmpz_multimod_naive_reduce(mp_limb_t * out, const fmpz_t A, const fmpz_multimod_naive_t mmod)
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
            
    for (i = 0; i < mmod->num_primes; i++)
        out[i] = _fmpz_reduce(A, mmod->powers_of_two[i], mmod->mod[i]);
}

