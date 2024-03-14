#include <flint/flint.h>
#include <flint/fmpz.h>

#include "nmod_vec_extra.h"
#include "fmpz_extra.h"

/* ------------------------------------------------------------ */
/* compute A s.t. A mod primes[j] = m[j] for all j              */
/* ------------------------------------------------------------ */
void 
fmpz_CRT_naive_CRT(fmpz_t A, mp_srcptr m, const fmpz_CRT_naive_t mCRT)
{
    fmpz_t comb;
    mp_ptr m_premul;
    ulong i;

    fmpz_init(comb);
    m_premul = _nmod_vec_init(mCRT->num_primes);
    for (i = 0; i < mCRT->num_primes; i++)
	m_premul[i] = nmod_mul(m[i], mCRT->inverses[i], mCRT->mod[i]);

    fmpz_CRT_naive_combine(comb, m_premul, mCRT);
    fmpz_mod(A, comb, mCRT->prod);

    fmpz_clear(comb);
    _nmod_vec_clear(m_premul);
}
