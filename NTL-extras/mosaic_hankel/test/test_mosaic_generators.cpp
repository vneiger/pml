#include <NTL/vec_lzz_p.h>
#include <NTL/mat_lzz_p.h>
#include <assert.h>

#include "mosaic_hankel.h"

NTL_CLIENT

void do_Z(Mat<zz_p>& M, long n, const zz_p& c){
  M.SetDims(n, n);
  for (long i = 0; i < n; i++)
    for (long j = 0; j < n; j++)
      M[i][j] = 0;
  for (long i = 0; i < n-1; i++)
    M[i+1][i] = 1;
  M[0][n-1] = c;
}

void random(Vec<zz_p>& v, long n, long m){
  v.SetLength(n+m-1);
  for (long i = 0; i < n+m-1; i++)
    v[i] = random_zz_p();
}

/*------------------------------------------------------------*/
/* if opt = 1, runs a check                                   */
/* else, runs timings                                         */
/*------------------------------------------------------------*/
void check(int opt){
  zz_p::FFTInit(0);
  for (long i = 2; i < 110; i += 1){

    Vec<zz_p> dat00, dat01, dat02, dat10, dat11, dat12;

    random(dat00, 2, i);
    random(dat01, 2, 2);
    random(dat02, 2, i);
    random(dat10, i-1, i);
    random(dat11, i-1, 2);
    random(dat12, i-1, i);
    
    hankel h00(dat00, 2, i), h01(dat01, 2, 2), h02(dat02, 2, i), h10(dat10, i-1, i), h11(dat11, i-1, 2), h12(dat12, i-1, i);
    
    Vec<hankel> row0, row1;
    
    row0.SetLength(3);
    row0[0] = h00;
    row0[1] = h01;
    row0[2] = h02;
    row1.SetLength(3);
    row1[0] = h10;
    row1[1] = h11;
    row1[2] = h12;
    Vec< Vec<hankel> > H;
    H.SetLength(2);
    H[0] = row0;
    H[1] = row1;
    
    mosaic_hankel MH(H);
    
    if (opt == 1){

      Mat<zz_p> G, H;
      generators(G, H, MH);

      Mat<zz_p> Z1, Z0, M_dense;
      do_Z(Z1, MH.NumRows(), to_zz_p(1));
      do_Z(Z0, MH.NumCols(), to_zz_p(0));
      to_dense(M_dense, MH);

      Mat<zz_p> res = Z1*M_dense-M_dense*transpose(Z0);
      assert (res == G*transpose(H));
      cout << i << " ";
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
