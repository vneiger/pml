#ifndef __VEC_LZZ_P_EXTRA_H
#define __VEC_LZZ_P_EXTRA_H

#include <NTL/mat_lzz_p.h>
#include <NTL/vec_lzz_p.h>

NTL_CLIENT

/*---------------------------------------------*/
/* random vector of length d                   */
/*---------------------------------------------*/
void random_vec_zz_p(Vec<zz_p>& A, long d);

/*---------------------------------------------*/
/* random matrix of size (d,e)                 */
/*---------------------------------------------*/
void random_mat_zz_p(mat_zz_p& A, long d, long e);

/*---------------------------------------------*/
/* inverts every entry in A -- TODO: CHECK 0   */
/*---------------------------------------------*/
void inv_naive(Vec<zz_p>& invA, const Vec<zz_p>& A);

/*---------------------------------------------*/
/* inverts every entry in A -- TODO: CHECK 0   */
/*---------------------------------------------*/
void inv(Vec<zz_p>& invA, const Vec<zz_p>& A);

/*---------------------------------------------*/
/* inverts every entry in A -- TODO: CHECK 0   */
/*---------------------------------------------*/
void inv(Vec<zz_p>& A);

#ifdef NTL_HAVE_LL_TYPE
/*---------------------------------------------*/
/* Inner product adapted from NTL              */
/*---------------------------------------------*/
long InnerProd_LL(const long *ap, const long *bp, long n);
long InnerProd_L(const long *ap, const long *bp, long n);
#endif

#endif
