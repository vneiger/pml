#ifndef __VEC_LZZ_P_EXTRA_H
#define __VEC_LZZ_P_EXTRA_H

#include <NTL/mat_lzz_p.h>
#include <NTL/vec_lzz_p.h>

NTL_CLIENT

/*------------------------------------------------------------*/
/* random vector of length d                                  */
/*------------------------------------------------------------*/
void random_vec_zz_p(Vec<zz_p>& A, long d);

/*------------------------------------------------------------*/
/* random matrix of size (d,e)                                */
/*------------------------------------------------------------*/
void random_mat_zz_p(mat_zz_p& A, long d, long e);

/*------------------------------------------------------------*/
/* inverts every entry in A -- TODO: CHECK 0                  */
/*------------------------------------------------------------*/
void inv_naive(Vec<zz_p>& invA, const Vec<zz_p>& A);

/*------------------------------------------------------------*/
/* inverts every entry in A -- TODO: CHECK 0                  */
/*------------------------------------------------------------*/
void inv(Vec<zz_p>& invA, const Vec<zz_p>& A);

/*------------------------------------------------------------*/
/* inverts every entry in A -- TODO: CHECK 0                  */
/*------------------------------------------------------------*/
void inv(Vec<zz_p>& A);

#endif

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
