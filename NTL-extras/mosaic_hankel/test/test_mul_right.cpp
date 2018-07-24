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

   for (long i1 = 1; i1 < 50; i1 += 1)
     for (long i2 = 1; i2 < 50; i2 += 1){

     if (opt == 1){
       Vec<zz_p> dat;
       dat.SetLength(i1+i2-1);
       for (long j = 0; j < i1+i2-1; j++)
	 dat[j] = j;
       
       hankel h(dat, i1, i2);

       Vec<zz_p> input, output;
       input.SetLength(i2);
       for (long j = 0; j < i2; j++)
	 input[j] = random_zz_p();
       output.SetLength(i1);

       Vec<long> input_long, output_long;
       conv(input_long, input);
       output_long.SetLength(i1);

       mul_right(output_long._vec__rep.rep, h, input_long._vec__rep.rep);

       conv(output, output_long);

       Mat<zz_p> M;
       to_dense(M, h);
       Vec<zz_p> output2 = M*input;

       assert (output == output2);
       cout << i1 << " " << i2 << endl;
     }
     else{
       cout << i1 << " ";

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
