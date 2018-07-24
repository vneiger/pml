#include <NTL/vec_lzz_p.h>
#include <NTL/mat_lzz_p.h>
#include <assert.h>

#include "vec_lzz_p_extra.h"
#include "cauchy_geometric_special.h"

NTL_CLIENT

/*----------------------------------------------------*/
/* ensures that G and H have generic rank profiles    */
/*----------------------------------------------------*/
void generate_generic_rank_profile(Mat<zz_p>& G, Mat<zz_p>& H, unsigned long rk){
  long n = G.NumRows();
  long m = H.NumRows();
  long alpha = G.NumCols();

  for (unsigned long i = 0; i < rk; i++)
    for (long j = 0; j < alpha; j++)
      G[i][j] = random_zz_p();
  for (long i = rk; i < n; i++)
    for (long j = 0; j < alpha; j++)
      G[i][j] = 0;

  for (unsigned long i = 0; i < rk; i++)
    for (long j = 0; j < alpha; j++)
      H[i][j] = random_zz_p();
  for (long i = rk; i < m; i++)
    for (long j = 0; j < alpha; j++)
      H[i][j] = 0;
}


/*----------------------------------------------------*/
/* ensures that G and H have not generic rank profile */
/* puts some rows fulls of zeros, then a non-zero one */
/*----------------------------------------------------*/
void generate_non_generic_rank_profile(Mat<zz_p>& G, Mat<zz_p>& H, unsigned long rk, unsigned long extra){
  long alpha = G.NumCols();

  if (extra <= rk){
    cout << rk <<  " " << extra << endl;
    Error("Matrix will have generic rank profile");
  }
  generate_generic_rank_profile(G, H, rk);

  for (long j = 0; j < alpha; j++){
    G[extra][j] = random_zz_p();
    H[extra][j] = random_zz_p();
  }

}


/*----------------------------------------------------*/
/* test in size (i,j), with rank rk                   */
/* tests a generic rank profile matrix                */
/* tests a non-generic rank profile one               */
/*----------------------------------------------------*/
void test(long i, long j, long alpha, unsigned long rk, int opt){

  if ((long)rk < 0 || (long)rk > i || (long)rk > j){
    cout << "# " << i << " " << j << " nothing to do\n";
    return;
  }

  zz_p a = random_zz_p();
  cauchy_geometric_special C(to_zz_p(1), power(a, i), a, i, j);
  Mat<zz_p> A, B;

  A.SetDims(i, alpha);
  B.SetDims(j, alpha);

  cauchy_like_geometric_special M, invM;
  Mat<zz_p> Z, invZ;
  long r;

  if (opt == 1){

    if (true) {
      generate_generic_rank_profile(A, B, rk);
      M = cauchy_like_geometric_special(A, B, C);
      to_dense(Z, M);
      r = invert(invM, M);
      to_dense(invZ, invM);
      Mat<zz_p> Z2 = Z;
      assert (r == gauss(Z2));
      Mat<zz_p> subZ;
      subZ.SetDims(r, r);
      for (long u = 0; u < r; u++)
      	for (long v = 0; v < r; v++)
      	  subZ[u][v] = Z[u][v];
      subZ = subZ*invZ;
      assert (IsIdent(subZ, r));
      cout << "# " << i << " " << j << " grp\n";
    }
    if ((long)rk < min(i, j)-1) {
      if (((long)rk+5 < i) && ((long)rk+5 < j))
    	generate_non_generic_rank_profile(A, B, rk, rk + 5);
      else
    	generate_non_generic_rank_profile(A, B, rk, min(i-1, j-1));
      M = cauchy_like_geometric_special(A, B, C);
      to_dense(Z, M);
      r = invert(invM, M);
      assert (r == -1);
      cout << "# " << i << " " << j << " not grp\n";
    }
  }
  else{
    cout << i << " ";
    
    double t;
    
    generate_generic_rank_profile(A, B, rk);
    M = cauchy_like_geometric_special(A, B, C);
    t = GetTime();
    for (long j = 0; j < 10; j++)
      invert(invM, M);
    t = GetTime() - t;
    cout << t << " ";
    
    cout << endl;
  }
}


/*------------------------------------------------------------*/
/* if opt = 1, runs a check                                   */
/* else, runs timings                                         */
/*------------------------------------------------------------*/
void check(int opt){
  zz_p::FFTInit(0);

   cout << "# " << zz_p::modulus() << endl;
   for (long i = 3; i < 200; i += 1){
     long alpha = 4;
     if (opt == 1){
       long j;
       j = i-2;
       test(i, j, alpha, 5, 1);
       test(i, j, alpha, min(i,j)/2, 1);
       test(i, j, alpha, min(i,j)-7, 1);
       j = i+2;
       test(i, j, alpha, 5, 1);
       test(i, j, alpha, min(i,j)/2, 1);
       test(i, j, alpha, min(i,j)-7, 1);
     }
     else {
       test(i, i, alpha, i, 0);
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
