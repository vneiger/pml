#include <NTL/lzz_pX.h>
#include <NTL/mat_lzz_p.h>

#include "lzz_pX_middle_product.h"
#include "mosaic_hankel.h"

NTL_CLIENT

/*----------------------------------------------------*/
/*----------------------------------------------------*/
/* Hankel matrices                                    */
/* stored as                                          */
/*       a5 a4 a3 a2                                  */
/*       a4 a3 a2 a1                                  */
/*       a3 a2 a1 a0                                  */
/*----------------------------------------------------*/
/*----------------------------------------------------*/


/*----------------------------------------------------*/
/* turns M into a dense matrix                        */
/*----------------------------------------------------*/
void to_dense(Mat<zz_p>& Mdense, const hankel& M){
  long n = M.NumRows();
  long m = M.NumCols();
  Mdense.SetDims(n, m);

  for (long i = 0; i < n; i++)
    for (long j = 0; j < m; j++)
      Mdense[i][j] = M(i,j);
}

/*----------------------------------------------------*/
/* right multiplication                               */
/*----------------------------------------------------*/
void mul_right(Vec<zz_p>& res, const hankel& M, const Vec<zz_p>& input){

  long nM = M.NumRows();
  long mM = M.NumCols();
  res.SetLength(nM);

  if (min(nM, mM) <= NTL_zz_pX_MUL_CROSSOVER){
    long sp = Kar_stk_size(max(nM, mM));
    zz_p *stk = new zz_p[sp];
    tKarMul_aux(res._vec__rep.rep, nM, input._vec__rep.rep, mM, M.data_rev._vec__rep.rep, nM+mM-1, stk);
    delete[] stk;
  }
  else{
    long K = NextPowerOfTwo(nM+mM-1);
    fftRep fft_input = fftRep(INIT_SIZE, K);
    
    zz_pX input_X, output_X;
    input_X.rep.SetLength(mM);
    zz_p *cf = input_X.rep.elts();
    for (long i = 0; i < mM; i++)
      cf[i] = input[mM-1-i];
    input_X.normalize();


    TofftRep(fft_input, input_X, K);
    mul(fft_input, fft_input, M.fft_data);
    FromfftRep(res._vec__rep.rep, fft_input, mM-1, nM+mM-2);
  }
}

/*----------------------------------------------------*/
/* left multiplication                                */
/*----------------------------------------------------*/
void mul_left(Vec<zz_p>& res, const hankel& M, const Vec<zz_p>& input){

  long nM = M.NumRows();
  long mM = M.NumCols();

  res.SetLength(mM);
  if (min(nM, mM) <= NTL_zz_pX_MUL_CROSSOVER){
    long sp = Kar_stk_size(max(nM, mM));
    zz_p *stk = new zz_p[sp];
    tKarMul_aux(res._vec__rep.rep, mM, input._vec__rep.rep, nM, M.data_rev._vec__rep.rep, nM+mM-1, stk);
    delete[] stk;
  }
  else{
    long K = NextPowerOfTwo(nM+mM-1);
    fftRep fft_input = fftRep(INIT_SIZE, K);
    
    zz_pX input_X, output_X;
    input_X.rep.SetLength(nM);
    zz_p *cf = input_X.rep.elts();
    for (long i = 0; i < nM; i++)
      cf[i] = input[nM-1-i];
    input_X.normalize();

    TofftRep(fft_input, input_X, K);
    mul(fft_input, fft_input, M.fft_data);
    FromfftRep(res._vec__rep.rep, fft_input, nM-1, nM+mM-2);
  }
}

