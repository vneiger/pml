#include <NTL/vec_lzz_p.h>
#include <NTL/mat_lzz_p.h>

#include "vec_lzz_p_extra.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/*                   Helper functions                         */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* inverts every entry in A                                   */
/*------------------------------------------------------------*/
void inv_naive(Vec<zz_p>& invA, const Vec<zz_p>& A)
{
    invA.SetLength(A.length());
    for (long i = 0; i < A.length(); ++i)
        inv(invA[i], A[i]);
}

/*------------------------------------------------------------*/
/* inverts every entry in A                                   */
/*------------------------------------------------------------*/
void inv(Vec<zz_p>& invA, const Vec<zz_p>& A)
{
    if (&invA == &A)
    {
        invA = inv(A);
        return;
    }

    const long n = A.length();
    invA.SetLength(n);

    if (n == 0)
        return;

    Vec<zz_p> tmp(INIT_SIZE, n);

    // tmp[i] = A[0] * A[1] * ... * A[i]
    tmp[0] = A[0];
    for (long i = 1; i < n; ++i)
        mul(tmp[i], tmp[i-1], A[i]);

    // aux = (A[0] * A[1] * ... * A[i])^{-1}
    zz_p aux = inv(tmp[n-1]);

    for (long i = n-1; i >= 1; --i)
    {
        mul(invA[i], aux, tmp[i-1]);
        mul(aux, aux, A[i]);
    }
    invA[0] = aux;
}

/*------------------------------------------------------------*/
/* builds the vector of mulmod_precon_t                       */
/*------------------------------------------------------------*/
void precomp(Vec<mulmod_precon_t>& out, const Vec<zz_p>& in)
{
    const long p = zz_p::modulus();
    const mulmod_t pinv = zz_p::ModulusInverse();
    out.SetLength(in.length());
    for (long i = 0; i < in.length(); ++i)
        out[i] = PrepMulModPrecon(in[i]._zz_p__rep, p, pinv);
}




// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
