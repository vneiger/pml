#include <NTL/vec_lzz_p.h>
#include <NTL/mat_lzz_p.h>

NTL_CLIENT

/* ------------------------------------------------------------ */
/* ------------------------------------------------------------ */
/* ------------------------------------------------------------ */
/*                     Helper functions                         */
/* ------------------------------------------------------------ */
/* ------------------------------------------------------------ */
/* ------------------------------------------------------------ */

/*---------------------------------------------*/
/* random vector of length d                   */
/*---------------------------------------------*/
void random(vec_zz_p& A, long d){
  A.SetLength(d);
  for (long i = 0; i < A.length(); i++)
    A[i] = random_zz_p();
}

/*---------------------------------------------*/
/* random matrix of size (d,e)                 */
/*---------------------------------------------*/
void random(mat_zz_p& A, long d, long e){
  A.SetDims(d, e);
  for (long i = 0; i < d; i++)
    for (long j = 0; j < e; j++)
      A[i][j] = random_zz_p();
}

/*---------------------------------------------*/
/* inverts every entry in A -- TODO: CHECK 0   */
/*---------------------------------------------*/
void inv_naive(vec_zz_p& invA, const vec_zz_p& A){
  long n = A.length();
  invA.SetLength(n);

  for (long i = 0; i < n; i++)
    invA[i] = 1/A[i];
}

/*---------------------------------------------*/
/* inverts every entry in A -- TODO: CHECK 0   */
/*---------------------------------------------*/
void inv(vec_zz_p& invA, const vec_zz_p& A){
  long n = A.length();
  vec_zz_p tmp;
  tmp.SetLength(n);
  invA.SetLength(n);

  if (n == 0)
    return;

  tmp[0] = A[0];
  for (long i = 1; i < n; i++)
    tmp[i] = tmp[i-1]*A[i];
  zz_p aux = 1/tmp[n-1];
  for (long i = n-1; i >= 1; i--){
    invA[i] = aux*tmp[i-1];
    aux *= A[i];
  }
  invA[0] = aux;
}

/*---------------------------------------------*/
/* inverts every entry in A -- TODO: CHECK 0   */
/*---------------------------------------------*/
void inv(vec_zz_p& A){
  vec_zz_p tmp = A;
  inv(A, tmp);
}



/*---------------------------------------------*/
/* Inner product adapted from NTL              */
/*---------------------------------------------*/

#ifdef NTL_HAVE_LL_TYPE

#define NTL_CAST_ULL(x) ((NTL_ULL_TYPE) (x))
#define NTL_MUL_ULL(x,y) (NTL_CAST_ULL(x)*NTL_CAST_ULL(y))

long InnerProd_LL(const long *ap, const long *bp, long n){
   const long BLKSIZE = (1L << min(20, 2*(NTL_BITS_PER_LONG-NTL_SP_NBITS)));

   unsigned long acc0 = 0;
   NTL_ULL_TYPE acc21 = 0;

   long i;
   for (i = 0; i <= n-BLKSIZE; i += BLKSIZE, ap += BLKSIZE, bp += BLKSIZE) {       // sum ap[j]*rep(bp[j]) for j in [0..BLKSIZE)
      NTL_ULL_TYPE sum = 0;
      for (long j = 0; j < BLKSIZE; j += 4) {
         sum += NTL_MUL_ULL(ap[j+0], bp[j+0]);
         sum += NTL_MUL_ULL(ap[j+1], bp[j+1]);
         sum += NTL_MUL_ULL(ap[j+2], bp[j+2]);
         sum += NTL_MUL_ULL(ap[j+3], bp[j+3]);
      }
      sum += acc0; 
      acc0 = sum;
      acc21 += (unsigned long) (sum >> NTL_BITS_PER_LONG);
   }
   if (i < n) {       // sum ap[i]*rep(bp[j]) for j in [0..n-i)
      NTL_ULL_TYPE sum = 0;
      long j = 0;
      for (; j <= n-i-4; j += 4) {
         sum += NTL_MUL_ULL(ap[j+0], bp[j+0]);
         sum += NTL_MUL_ULL(ap[j+1], bp[j+1]);
         sum += NTL_MUL_ULL(ap[j+2], bp[j+2]);
         sum += NTL_MUL_ULL(ap[j+3], bp[j+3]);
      }
      for (; j < n-i; j++)
         sum += NTL_MUL_ULL(ap[j], bp[j]);
      sum += acc0; 
      acc0 = sum;
      acc21 += (unsigned long) (sum >> NTL_BITS_PER_LONG);
   }

   sp_ll_reduce_struct dinv = zz_p::ll_red_struct();
   if (dinv.nbits == NTL_SP_NBITS) 
     return sp_ll_red_31_normalized(acc21 >> NTL_BITS_PER_LONG, acc21, acc0, zz_p::modulus(), dinv);
   else
     return sp_ll_red_31(acc21 >> NTL_BITS_PER_LONG, acc21, acc0, zz_p::modulus(), dinv);
}


long InnerProd_L(const long *ap, const long *bp, long n){
   unsigned long sum = 0;
   long j = 0;

   for (; j <= n-4; j += 4) {
      sum += (ap[j+0]) * (bp[j+0]);
      sum += (ap[j+1]) * (bp[j+1]);
      sum += (ap[j+2]) * (bp[j+2]);
      sum += (ap[j+3]) * (bp[j+3]);
   }

   for (; j < n; j++)
      sum += (ap[j]) * (bp[j]);

   return rem(sum, zz_p::modulus(), zz_p::red_struct());
}

#endif
