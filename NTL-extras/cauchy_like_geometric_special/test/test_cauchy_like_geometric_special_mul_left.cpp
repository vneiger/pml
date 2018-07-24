#include <NTL/vec_lzz_p.h>
#include <assert.h>

#include "vec_lzz_p_extra.h"
#include "cauchy_geometric_special.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* if opt = 1, runs a check                                   */
/* else, runs timings                                         */
/*------------------------------------------------------------*/
void check(int opt){
   zz_p::FFTInit(0);

   for (long i = 1; i < 10000; i += 1){

     zz_p a = random_zz_p();
     long j = i+10;
     long alpha = 4;
     cauchy_geometric_special C(to_zz_p(1), power(a, i), a, i, j);
     mat_zz_p A, B;
     random(A, i, alpha);
     random(B, j, alpha);
     cauchy_like_geometric_special M(A, B, C);
     mat_zz_p Z;

     vec_zz_p in, out;
     random(in, i);

     if (opt == 1){
       vec_zz_p out2;
       mul_left(out, M, in);
       to_dense(Z, M);
       mul(out2, in, Z);
       assert (out2 == out);
     }
     else{
       cout << i << " ";

       double t;

       t = GetTime();
       for (long j = 0; j < 10000; j++)
	 mul_left(out, M, in);
       t = GetTime() - t;
       cout << t << " ";

       cout << endl;
     }
   }
}

/*------------------------------------------------------------*/
/* main just calls check()                                    */
/* if not argument is given, runs timings                     */
/* if the argument 1 is given, runs check                     */
/*------------------------------------------------------------*/
int main(int argc, char** argv){

#ifdef NTL_HAVE_AVX
#define MAX_DBL_INT ((1L << NTL_DOUBLE_PRECISION)-1)
  cout << MAX_DBL_INT << endl;
  return 0;
#endif

  int opt = 0;
  if (argc > 1)
    opt = atoi(argv[1]);
  check(opt);

  return 0;
}
