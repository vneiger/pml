#ifndef __VEC_LZZ_P_EXTRA_H
#define __VEC_LZZ_P_EXTRA_H

#include <NTL/mat_lzz_p.h>
#include <NTL/vec_lzz_p.h>
#include <NTL/version.h>

NTL_CLIENT

/*------------------------------------------------------------*/
/* random vector of length d                                  */
/*------------------------------------------------------------*/
void random_vec_zz_p(Vec<zz_p>& A, long d);

#if ( (NTL_MAJOR_VERSION < 10) || ((NTL_MAJOR_VERSION == 10) && (NTL_MINOR_VERSION < 4)) )
inline Vec<zz_p> random_vec_zz_p(long n)
{ 
    Vec<zz_p> x; 
    random_vec_zz_p(x, n); 
    return x;
}
#endif

/*------------------------------------------------------------*/
/* random matrix of size (d,e)                                */
/*------------------------------------------------------------*/
void random_mat_zz_p(mat_zz_p& A, long d, long e);

#if ( (NTL_MAJOR_VERSION < 10) || ((NTL_MAJOR_VERSION == 10) && (NTL_MINOR_VERSION < 4)) )
inline Mat<zz_p> random_mat_zz_p(long d, long e)
{ 
    Mat<zz_p> x; 
    random_mat_zz_p(x, d, e); 
    return x;
}
#endif

/*------------------------------------------------------------*/
/* inverts every entry in A                                   */
/* output can alias input                                     */
/*------------------------------------------------------------*/
void inv_naive(Vec<zz_p>& invA, const Vec<zz_p>& A);

/*------------------------------------------------------------*/
/* inverts every entry in A                                   */
/* output can alias input                                     */
/*------------------------------------------------------------*/
void inv(Vec<zz_p>& invA, const Vec<zz_p>& A);

inline Vec<zz_p> inv(const Vec<zz_p>& A)
{
    Vec<zz_p> x;
    inv(x, A);
    return x;
}

/*------------------------------------------------------------*/
/* builds the vector of mulmod_precon_t                       */
/*------------------------------------------------------------*/
void precomp(Vec<mulmod_precon_t>& out, const Vec<zz_p>& in);

inline Vec<mulmod_precon_t> precomp(const Vec<zz_p>& in)
{
    Vec<mulmod_precon_t> out;
    precomp(out, in);
    return out;
}


#endif

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
