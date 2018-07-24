#include <NTL/lzz_pX.h>

#include "cauchy_geometric_special.h"

NTL_CLIENT

/*-------------------------------------------------*/
/*-------------------------------------------------*/
/* matrix / matrix product for Cauchy-like special */
/*-------------------------------------------------*/
/*-------------------------------------------------*/

/*-------------------------------------------------*/
/* matrix / matrix product for Cauchy-like special */
/* A += C.B, where C is cauchy-like special        */
/* FFT multiplication                              */
/*-------------------------------------------------*/
void mat_addmul_fft(long *A, // A = the result
		    const long* B, // B = the argument
		    const long* G, const long* H, // generators of C
		    const Vec<mulmod_precon_t> & P_pre, 
		    const Vec<mulmod_precon_t> & Q_pre, 
		    const Vec<long> & P_r, 
		    const Vec<long> & Q_r, 
		    long shift_P, long shift_Q,
		    long N1, long N2, // dimensions (row, column) of C
		    long a, long b){ // displacement rank and number of columns in the input

  const long p = zz_p::modulus();
  const mulmod_t pinv = zz_p::ModulusInverse();

  zz_p *w = new zz_p[N1];

  long K = NextPowerOfTwo(N1+N2-1);
  fftRep fftR = fftRep(INIT_SIZE, K);
  fftRep fftUV = fftRep(INIT_SIZE, K); 

  long d = N1+N2-1;
  zz_pX revA, polU, polV;

  revA.rep.SetLength(d);
  zz_p* revAp = revA.rep.elts();
  for (long i = 0; i < d; i++)
    revAp[i] = Q_r[d-i-1+shift_Q];
  revA.normalize();
  TofftRep(fftR, revA, K);

  polV.rep.SetLength(N2); 
  zz_p *polVp = polV.rep.elts();

  long* stk = new long[N1];  

  for (long k = 0; k < b; k++){

    for (long i = 0; i < N1; i++)
      stk[i] = 0;

    for (long i = 0; i < a; i++){
      for (long j = 0, aj = 0, bj = 0; j < N2; j++, bj += b, aj += a)
	polVp[j]._zz_p__rep = MulMod(B[bj+k], H[aj+i], p, pinv);

      TofftRep(fftUV, polV, K);
      mul(fftUV, fftUV, fftR);
      FromfftRep(polU, fftUV, N2-1, d-1);

      const long Ul = polU.rep.length();
      const zz_p *polUp = polU.rep.elts();

      for (long j = 0, aj = 0; j < Ul; j++, aj += a)
      	stk[j] = AddMod(stk[j], MulMod(rep(polUp[j]), G[aj+i], p, pinv), p);
    }
    for (long j = 0, aj = 0, bj = 0; j < N1; j++, aj += a, bj += b)
      A[bj+k] = AddMod(A[bj+k], MulModPrecon(stk[j], P_r[j+shift_P], p, P_pre[j+shift_P]), p);
  }

  delete[] stk;
  delete[] w;
}

/*-------------------------------------------------*/
/* matrix / matrix product for Cauchy-like special */
/* A -= C.B, where C is cauchy-like special        */
/* FFT multiplication                              */
/*-------------------------------------------------*/
void mat_submul_fft(long *A, // A = the result
		    const long* B, // B = the argument
		    const long* G, const long* H, // generators of C
		    const Vec<mulmod_precon_t> & P_pre, 
		    const Vec<mulmod_precon_t> & Q_pre, 
		    const Vec<long> & P_r, 
		    const Vec<long> & Q_r, 
		    long shift_P, long shift_Q,
		    long N1, long N2, // dimensions (row, column) of C
		    long a, long b){ // displacement rank and number of columns in the input

  const long p = zz_p::modulus();
  const mulmod_t pinv = zz_p::ModulusInverse();

  zz_p *w = new zz_p[N1];

  long K = NextPowerOfTwo(N1+N2-1);
  fftRep fftR = fftRep(INIT_SIZE, K);
  fftRep fftUV = fftRep(INIT_SIZE, K); 

  long d = N1+N2-1;
  zz_pX revA, polU, polV;

  revA.rep.SetLength(d);
  zz_p* revAp = revA.rep.elts();
  for (long i = 0; i < d; i++)
    revAp[i] = Q_r[d-i-1+shift_Q];
  revA.normalize();
  TofftRep(fftR, revA, K);

  polV.rep.SetLength(N2); 
  zz_p *polVp = polV.rep.elts();

  long* stk = new long[N1];  

  for (long k = 0; k < b; k++){

    for (long i = 0; i < N1; i++)
      stk[i] = 0;

    for (long i = 0; i < a; i++){
      for (long j = 0, aj = 0, bj = 0; j < N2; j++, bj += b, aj += a)
	polVp[j]._zz_p__rep = MulMod(B[bj+k], H[aj+i], p, pinv);

      TofftRep(fftUV, polV, K);
      mul(fftUV, fftUV, fftR);
      FromfftRep(polU, fftUV, N2-1, d-1);

      const long Ul = polU.rep.length();
      const zz_p *polUp = polU.rep.elts();

      for (long j = 0, aj = 0; j < Ul; j++, aj += a)
      	stk[j] = AddMod(stk[j], MulMod(rep(polUp[j]), G[aj+i], p, pinv), p);
    }
    for (long j = 0, aj = 0, bj = 0; j < N1; j++, aj += a, bj += b)
      A[bj+k] = SubMod(A[bj+k], MulModPrecon(stk[j], P_r[j+shift_P], p, P_pre[j+shift_P]), p);
  }

  delete[] stk;
  delete[] w;
}

/*-------------------------------------------------*/
/* matrix / matrix product for Cauchy-like special */
/* A += C^t.B, where C is cauchy-like special      */
/* FFT multiplication                              */
/*-------------------------------------------------*/
inline void mat_addmul_t_fft(long *A, // A = the result
			     const long* B, // B = the argument
			     const long* G, const long* H, // generators of C
			     const Vec<mulmod_precon_t> & P_pre, 
			     const Vec<mulmod_precon_t> & Q_pre, 
			     const Vec<long> & P_r, 
			     const Vec<long> & Q_r, 
			     long shift_P, long shift_Q,
			     long N1, long N2, // dimensions (row, column) of C
			     long a, long b){ // displacement rank and number of columns in the input

  const long p = zz_p::modulus();
  const mulmod_t pinv = zz_p::ModulusInverse();

  long K = NextPowerOfTwo(N1+N2-1);
  fftRep fftS = fftRep(INIT_SIZE, K);
  fftRep fftUV = fftRep(INIT_SIZE, K); 

  long d = N1+N2-1;
  zz_pX polA, polV, polU;

  polA.rep.SetLength(d);
  zz_p* Aap = polA.rep.elts();
  for (long i = 0; i < d; i++)
    Aap[i] = Q_r[i+shift_Q];
  polA.normalize();

  TofftRep(fftS, polA, K);

  polV.rep.SetLength(N1); 
  zz_p *polVp = polV.rep.elts();

  long* stk = new long[N2];  
  long* Bk = new long[N1];  

  for (long k = 0; k < b; k++){

    for (long j = 0, bj = 0; j < N1; j++, bj += b)
      Bk[j] = MulModPrecon(B[bj+k], P_r[j+shift_P], p, P_pre[j+shift_P]);
    for (long j = 0; j < N2; j++)
      stk[j] = 0;

    for (long i = 0; i < a; i++){
      for (long j = 0, aj = 0; j < N1; j++, aj += a)
      	polVp[j].LoopHole() = MulMod(Bk[j], H[aj+i], p, pinv);
      
      TofftRep(fftUV, polV, K);
      mul(fftUV, fftUV, fftS);
      FromfftRep(polU, fftUV, N1-1, d-1);

      long Ul = polU.rep.length();
      zz_p *polUp = polU.rep.elts();

      for (long j = 0, aj = 0; j < Ul; j++, aj += a)
      	stk[j] = AddMod(stk[j], MulMod(rep(polUp[j]), G[aj+i], p, pinv), p);
    }
    for (long j = 0, bj = 0; j < N2; j++, bj += b)
      A[bj+k] = AddMod(A[bj+k], stk[j], p);
  }

  delete[] stk;
  delete[] Bk;
}

/*-------------------------------------------------*/
/* matrix / matrix product for Cauchy-like special */
/* A -= C^t.B, where C is cauchy-like special      */
/* FFT multiplication                              */
/*-------------------------------------------------*/
inline void mat_submul_t_fft(long *A, // A = the result
			     const long* B, // B = the argument
			     const long* G, const long* H, // generators of C
			     const Vec<mulmod_precon_t> & P_pre, 
			     const Vec<mulmod_precon_t> & Q_pre, 
			     const Vec<long> & P_r, 
			     const Vec<long> & Q_r, 
			     long shift_P, long shift_Q,
			     long N1, long N2, // dimensions (row, column) of C
			     long a, long b){ // displacement rank and number of columns in the input

  const long p = zz_p::modulus();
  const mulmod_t pinv = zz_p::ModulusInverse();

  long K = NextPowerOfTwo(N1+N2-1);
  fftRep fftS = fftRep(INIT_SIZE, K);
  fftRep fftUV = fftRep(INIT_SIZE, K); 

  long d = N1+N2-1;
  zz_pX polA, polV, polU;

  polA.rep.SetLength(d);
  zz_p* Aap = polA.rep.elts();
  for (long i = 0; i < d; i++)
    Aap[i] = Q_r[i+shift_Q];
  polA.normalize();

  TofftRep(fftS, polA, K);

  polV.rep.SetLength(N1); 
  zz_p *polVp = polV.rep.elts();

  long* stk = new long[N2];  
  long* Bk = new long[N1];  

  for (long k = 0; k < b; k++){

    for (long j = 0, bj = 0; j < N1; j++, bj += b)
      Bk[j] = MulModPrecon(B[bj+k], P_r[j+shift_P], p, P_pre[j+shift_P]);
    for (long j = 0; j < N2; j++)
      stk[j] = 0;

    for (long i = 0; i < a; i++){
      for (long j = 0, aj = 0; j < N1; j++, aj += a)
      	polVp[j].LoopHole() = MulMod(Bk[j], H[aj+i], p, pinv);
      
      TofftRep(fftUV, polV, K);
      mul(fftUV, fftUV, fftS);
      FromfftRep(polU, fftUV, N1-1, d-1);

      long Ul = polU.rep.length();
      zz_p *polUp = polU.rep.elts();

      for (long j = 0, aj = 0; j < Ul; j++, aj += a)
      	stk[j] = SubMod(stk[j], MulMod(rep(polUp[j]), G[aj+i], p, pinv), p);
    }
    for (long j = 0, bj = 0; j < N2; j++, bj += b)
      A[bj+k] = AddMod(A[bj+k], stk[j], p);
  }

  delete[] stk;
  delete[] Bk;
}



/*-------------------------------------------------*/
/*-------------------------------------------------*/
/* wrapper functions                               */
/*-------------------------------------------------*/
/*-------------------------------------------------*/

/*-------------------------------------------------*/
/* matrix / matrix product for Cauchy-like special */
/* A += C.B, where C is cauchy-like special        */
/*-------------------------------------------------*/
void mat_addmul(long *A, // A = the result
		const long* B, // B = the argument
		const long* G, const long* H, // generators of C
		const Vec<mulmod_precon_t> & P_pre, 
		const Vec<mulmod_precon_t> & Q_pre, 
		const Vec<long> & P_r, 
		const Vec<long> & Q_r, 
		long shift_P, long shift_Q,
		long N1, long N2, // dimensions (row, column) of C
		long a){ // displacement rank

  mat_addmul_fft(A, B, G, H, P_pre, Q_pre, P_r, Q_r, shift_P, shift_Q, N1, N2, a, a);
}

/*-------------------------------------------------*/
/* matrix / matrix product for Cauchy-like special */
/* A -= C.B, where C is cauchy-like special        */
/*-------------------------------------------------*/
void mat_submul(long *A, // A = the result
		const long* B, // B = the argument
		const long* G, const long* H, // generators of C
		const Vec<mulmod_precon_t> & P_pre, 
		const Vec<mulmod_precon_t> & Q_pre, 
		const Vec<long> & P_r, 
		const Vec<long> & Q_r, 
		long shift_P, long shift_Q,
		long N1, long N2, // dimensions (row, column) of C
		long a){ // displacement rank

  mat_submul_fft(A, B, G, H, P_pre, Q_pre, P_r, Q_r, shift_P, shift_Q, N1, N2, a, a);
}

/*-------------------------------------------------*/
/* matrix / matrix product for Cauchy-like special */
/* A += C^t.B, where C is cauchy-like special      */
/*-------------------------------------------------*/
void mat_addmul_t(long *A, // A = the result
		  const long* B, // B = the argument
		  const long* G, const long* H, // generators of C
		  const Vec<mulmod_precon_t> & P_pre, 
		  const Vec<mulmod_precon_t> & Q_pre, 
		  const Vec<long> & P_r, 
		  const Vec<long> & Q_r, 
		  long shift_P, long shift_Q,
		  long N1, long N2, // dimensions (row, column) of C
		  long a){ // displacement rank

  mat_addmul_t_fft(A, B, G, H, P_pre, Q_pre, P_r, Q_r, shift_P, shift_Q, N1, N2, a, a);
}

/*-------------------------------------------------*/
/* matrix / matrix product for Cauchy-like special */
/* A -= C^t.B, where C is cauchy-like special      */
/*-------------------------------------------------*/
void mat_submul_t(long *A, // A = the result
		  const long* B, // B = the argument
		  const long* G, const long* H, // generators of C
		  const Vec<mulmod_precon_t> & P_pre, 
		  const Vec<mulmod_precon_t> & Q_pre, 
		  const Vec<long> & P_r, 
		  const Vec<long> & Q_r, 
		  long shift_P, long shift_Q,
		  long N1, long N2, // dimensions (row, column) of C
		  long a){ // displacement rank

  mat_submul_t_fft(A, B, G, H, P_pre, Q_pre, P_r, Q_r, shift_P, shift_Q, N1, N2, a, a);
}

/*---------------------------------------------------*/
/* computes M*input                                  */
/*---------------------------------------------------*/
void mul_right(Vec<zz_p>& res, const cauchy_like_geometric_special& M, const Vec<zz_p>& input){
  long N1 = M.NumRows();
  long N2 = M.NumCols();
  long *A = new long[N1];
  long *B = new long[N2];
  
  for (long i = 0; i < N2; i++)
    B[i] = rep(input[i]);
  for (long i = 0; i < N1; i++)
    A[i] = 0;

  mat_addmul_fft(A, B, M.G._vec__rep, M.H._vec__rep, M.C.P_pre, M.C.Qab_pre, M.C.P_r, M.C.Qab_r, 0, 0, N1, N2, M.Alpha(), 1);

  res.SetLength(N1);
  for (long i = 0; i < N1; i++)
    res[i].LoopHole() = A[i];

  delete[] A;
  delete[] B;
}

/*---------------------------------------------------*/
/* computes M^t*input                                */
/*---------------------------------------------------*/
void mul_left(Vec<zz_p>& res, const cauchy_like_geometric_special& M, const Vec<zz_p>& input){
  long N1 = M.NumRows();
  long N2 = M.NumCols();
  long *A = new long[N2];
  long *B = new long[N1];
  
  for (long i = 0; i < N1; i++)
    B[i] = rep(input[i]);
  for (long i = 0; i < N2; i++)
    A[i] = 0;

  mat_addmul_t_fft(A, B, M.H._vec__rep, M.G._vec__rep, M.C.P_pre, M.C.Qab_pre, M.C.P_r, M.C.Qab_r, 0, 0, N1, N2, M.Alpha(), 1);

  res.SetLength(N2);
  for (long i = 0; i < N2; i++)
    res[i].LoopHole() = A[i];

  delete[] A;
  delete[] B;
}
