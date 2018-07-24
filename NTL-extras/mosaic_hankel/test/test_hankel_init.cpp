#include <NTL/vec_lzz_p.h>
#include <assert.h>

#include "mosaic_hankel.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* if opt = 1, runs a check                                   */
/* else, runs timings                                         */
/*------------------------------------------------------------*/
void check(int opt){
   zz_p::FFTInit(0);

   for (long i = 2; i < 10; i += 1){

     zz_p a = random_zz_p();

     if (opt == 1){
       Vec<zz_p> dat;
       dat.SetLength(2*i-1);
       for (long j = 0; j < 2*i-1; j++)
	 dat[j] = j;
       
       hankel h(dat, i, i);
       
       cout << h(0,0) << " " << h(0,1) << endl << h(1,0) << " " << h(1,1) << endl;
       cout << endl;
     }
     else{
       cout << i << " ";

       double t;

       t = GetTime();
       for (long j = 0; j < 10000; j++)
	 ;
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
