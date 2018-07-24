#include <NTL/vec_lzz_p.h>
#include <assert.h>

#include "vec_lzz_p_extra.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* if opt = 1, runs a check                                   */
/* else, runs timings                                         */
/*------------------------------------------------------------*/
void check(int opt){
   zz_p::FFTInit(0);

   for (long i = 0; i < 10000; i += 1){
     vec_zz_p A, invA1, invA2;
     random(A, i);
     
     if (opt == 1){
       inv(invA1, A);
       inv_naive(invA2, A);
       for (long j = 0; j < i; j++){
	 assert (invA1[j] == invA2[j]);
	 assert (invA1[j] == 1/A[j]);
       }
     }
     else{
       cout << i << " ";

       double t;

       t = GetTime();
       for (long j = 0; j < 10000; j++)
	 inv(invA1, A);
       t = GetTime() - t;
       cout << t << " ";

       t = GetTime();
       for (long j = 0; j < 10000; j++)
	 inv_naive(invA2, A);
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
  int opt = 0;
  if (argc > 1)
    opt = atoi(argv[1]);
  check(opt);

  return 0;
}
