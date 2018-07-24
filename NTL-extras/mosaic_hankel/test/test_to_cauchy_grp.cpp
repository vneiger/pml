#include <NTL/vec_lzz_p.h>
#include <NTL/mat_lzz_p.h>
#include <assert.h>

#include "mosaic_hankel.h"
#include "lzz_pX_CRT.h"
#include "cauchy_geometric_special.h"

NTL_CLIENT

void to_dense(Mat<zz_p>& M, const Vec<zz_p>& e){
  long n = e.length();
  M.SetDims(n, n);
  for (long i = 0; i < n; i++)
    for (long j = 0; j < n; j++)
      M[i][j] = 0;
  for (long i = 0; i < n; i++)
    M[i][i] = e[i];
}

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
  //zz_p::init(9007);

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

      cauchy_like_geometric_special CL;
      Mat<zz_p> X, Y;
      Vec<zz_p> e, f;
      zz_pX_Multipoint_Geometric X_int, Y_int;

      to_cauchy_grp(CL, X_int, Y_int, e, f, MH);

      Mat<zz_p> MM, MC, MX_i, MY_i, Me, Mf;
      to_dense(MM, MH);
      to_dense(MC, CL);
      to_dense(MX_i, X_int);
      to_dense(MY_i, Y_int);
      to_dense(Me, e);
      to_dense(Mf, f);

      Mat<zz_p> T = Me*MX_i*MM*transpose(Mf*MY_i);
      assert (T == MC);

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
