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
/* random vector of length d                                  */
/*------------------------------------------------------------*/
void random_vec_zz_p(Vec<zz_p>& A, long d)
{
    A.SetLength(d);
    for (long i = 0; i < A.length(); i++)
        A[i] = random_zz_p();
}

/*------------------------------------------------------------*/
/* random matrix of size (d,e)                                */
/*------------------------------------------------------------*/
void random_mat_zz_p(mat_zz_p& A, long d, long e)
{
    A.SetDims(d, e);
    for (long i = 0; i < d; i++)
        for (long j = 0; j < e; j++)
            A[i][j] = random_zz_p();
}

/*------------------------------------------------------------*/
/* inverts every entry in A                                   */
/*------------------------------------------------------------*/
void inv_naive(Vec<zz_p>& invA, const Vec<zz_p>& A)
{
    long n = A.length();
    invA.SetLength(n);

    for (long i = 0; i < n; i++)
        invA[i] = 1/A[i];
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
        
    long n = A.length();
    Vec<zz_p> tmp;
    tmp.SetLength(n);
    invA.SetLength(n);

    if (n == 0)
        return;

    tmp[0] = A[0];
    for (long i = 1; i < n; i++)
        tmp[i] = tmp[i-1]*A[i];
    zz_p aux = 1/tmp[n-1];
    for (long i = n-1; i >= 1; i--)
    {
        invA[i] = aux*tmp[i-1];
        aux *= A[i];
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
    long m = in.length();
    out.SetLength(m);
    const zz_p* in_rep = in.elts();
    for (long i = 0; i < m; i++)
    {
        out[i] = PrepMulModPrecon(in_rep[i]._zz_p__rep, p, pinv);
    }
}




// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
