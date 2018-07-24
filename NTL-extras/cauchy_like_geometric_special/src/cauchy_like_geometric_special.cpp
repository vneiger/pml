#include <NTL/mat_lzz_p.h>
#include <NTL/lzz_pX.h>

#include "vec_lzz_p_extra.h"
#include "cauchy_geometric_special.h"

NTL_CLIENT


/*-------------------------------------------*/
/* sample runs, with approximate thresholds  */
/* 4  100                                    */
/* 9  250                                    */
/* 14  250                                   */
/* 19  340                                   */
/* 24  400                                   */
/* 29  400                                   */
/* 34  430                                   */
/* 39  550                                   */
/* 44  640                                   */
/* 49  670                                   */
/*-------------------------------------------*/
static 
long threshold(long alpha){
  return 100 + alpha*12;
}

/*------------------------------------------------------------*/
/* A is such that DA - AD' = G.H^t                            */
/*------------------------------------------------------------*/
void to_dense(Mat<zz_p> & A, const cauchy_like_geometric_special& CL){
  long n = CL.NumRows(), m = CL.NumCols(), alpha = CL.Alpha();
  
  A.SetDims(n, m);
  vec_zz_p tmpg, tmph;
  tmpg.SetLength(alpha);
  tmph.SetLength(alpha);
  for (long i = 0; i < n; i++){
    for (long ell = 0; ell < alpha; ell++)
      tmpg[ell] = CL.G[i*alpha+ell];
    for (long j = 0; j < m; j++){
      for (long ell = 0; ell < alpha; ell++)
	tmph[ell] = CL.H[j*alpha+ell];
      InnerProduct(A[i][j], tmpg, tmph);
      A[i][j] /= (CL.C.a*power(CL.C.q, i)-CL.C.b*power(CL.C.q, j));
    }
  }
}

// /*------------------------------------------------------------*/
// /* inverts a special-Cauchy-like matrix                       */
// /* assumes generic rank profile                               */
// /* the generators are Yp, Zp                                  */
// /* alpha = displacement rank                                  */
// /* overwrites Yp, Zp                                          */
// /* i0, j0 are the initial indices in i, j resp.               */
// /* i1, j1 are the end indices in i, j resp.                   */  
// /*------------------------------------------------------------*/
// static long invert_raw(long *Yp, long *Zp, const long alpha,
// 		       const cauchy_geometric_special& C,
// 		       const long i0, const long i1, const long j0, const long j1){


//   // P[i+i0]*Qab[j+j0-i-i0+n-1] = 1/(a*power(q,i+i0)-b*power(q,j+j0))
//   // P[i+i0]*Qa[j+j0-i-i0+n-1] = 1/(a*power(q,i+i0)-a*power(q,j+j0)) = 1/(x[i]-x[j])
//   // P[i+i0]*Qb[j+j0-i-i0+n-1] = 1/(b*power(q,i+i0)-b*power(q,j+j0)) = 1/(y[i]-y[j])

//   const long n = C.n;
//   const long p = zz_p::modulus();
//   mulmod_precon_t *Ypre = new mulmod_precon_t[alpha];
//   mulmod_precon_t *Zpre = new mulmod_precon_t[alpha];
//   const mulmod_t pinv = zz_p::ModulusInverse();

//   long *tmpYk = Yp;
//   long *tmpZk = Zp;

//   cout << "enter with: " << i0 << " " << i1 << " " << j0 << " " << j1 << endl;
//   long N = min(i1-i0, j1-j0);

//   for (long k = 0; k < N; k++){
//     long *tmpYi = Yp;
//     long *tmpZi = Zp;

//     for (long j = 0; j < alpha; j++){
//       Ypre[j] = PrepMulModPrecon(tmpYk[j], p, pinv);
//       Zpre[j] = PrepMulModPrecon(tmpZk[j], p, pinv);
//     }

//     // computes the top-left entry in the current submatrix
//     long rd = 0;
//     for (long j = 0; j < alpha; j++)
//       rd = AddMod(rd, MulMod(tmpZk[j], tmpYk[j], p, pinv), p);


//     // early exit if rank(M) < n
//     if (rd == 0){

//       cout << "with zero pivot with k=" << k<< endl;

//       long *tmpY = tmpYk;
//       for (long a = k; a < i1-i0; a++){
//     	long *tmpZ = tmpZk;
//     	for (long b = k; b < j1-j0; b++){
//     	  long tmp = 0;
//     	  for (long j = 0; j < alpha; j++)
//     	    tmp = AddMod(tmp, MulMod(tmpZ[j], tmpY[j], p, pinv), p);
//     	  if (tmp != 0) {
// 	    cout << "failed " << a << " " << b << endl;
//     	    delete[] Ypre;
//     	    delete[] Zpre;
//     	    return -1;  // not generic rank profile
//     	  }
// 	  cout << "passed " << a << " " << b << endl;
//     	  tmpZ += alpha;
//     	}
//     	tmpY += alpha;
//       }
//       delete[] Ypre;
//       delete[] Zpre;
//       return k; // OK, rank = k
//     }

//     rd = MulModPrecon(rd, C.P_r[k+i0], p, C.P_pre[k+i0]);
//     rd = MulModPrecon(rd, C.Qab_r[j0-i0+n-1], p, C.Qab_pre[j0-i0+n-1]);

//     zz_p d(rd, INIT_LOOP_HOLE);
//     zz_p id = 1/d;
//     zz_p mid = -id;

//     const long rid = rep(id);
//     const long rmid = rep(mid);
//     mulmod_precon_t pd = PrepMulModPrecon(rd, p, pinv);
//     mulmod_precon_t prid = PrepMulModPrecon(rid, p, pinv);

//     for (long i = 0; i < k; i++){

//       long ell = 0, u = 0;
//       for (long j = 0; j < alpha; j++){
//       	ell = AddMod(ell, MulModPrecon(tmpYi[j], tmpZk[j], p, Zpre[j]), p);
//       	u = AddMod(u, MulModPrecon(tmpZi[j], tmpYk[j], p, Ypre[j]), p);
//       }
//       // not so great for small alpha
//       // ell = InnerProd_LL(tmpYi, tmpZk, alpha);   
//       // u = InnerProd_LL(tmpZi, tmpYk, alpha);

//       ell = MulModPrecon(ell, rid, p, prid);
//       ell = MulModPrecon(ell, C.P_r[i+i0], p, C.P_pre[i+i0]);
//       ell = MulModPrecon(ell, C.Qb_r[k+j0-i-i0+n-1], p, C.Qb_pre[k+j0-i-i0+n-1]);

//       u = MulModPrecon(u, rid, p, prid);
//       u = MulModPrecon(u, C.P_r[i+i0], p, C.P_pre[i+i0]);
//       u = MulModPrecon(u, C.Qa_r[k+j0-i-i0+n-1], p, C.Qa_pre[k+j0-i-i0+n-1]);

//       for (long j = 0; j < alpha; j++){
//       	tmpYi[j] = SubMod(tmpYi[j], MulModPrecon(ell, tmpYk[j], p, Ypre[j]), p);
//       	tmpZi[j] = AddMod(tmpZi[j], MulModPrecon(u, tmpZk[j], p, Zpre[j]), p);
//       }

//       tmpYi += alpha;
//       tmpZi += alpha;
//     }

//     for (long j = 0; j < alpha; j++){
//       tmpYk[j] = MulModPrecon(rmid, tmpYk[j], p, Ypre[j]);
//       Ypre[j] = PrepMulModPrecon(tmpYk[j], p, pinv);
//       tmpZk[j] = MulModPrecon(rid, tmpZk[j], p, Zpre[j]);
//       Zpre[j] = PrepMulModPrecon(tmpZk[j], p, pinv);
//     }

//     tmpYi += alpha;
//     tmpZi += alpha;

//     for (long i = k+1; i < N; i++){

//       long ell = 0, u = 0;
//       for (long j = 0; j < alpha; j++){
// 	ell = AddMod(ell, MulModPrecon(tmpYi[j], tmpZk[j], p, Zpre[j]), p);
//         u = AddMod(u, MulModPrecon(tmpZi[j], tmpYk[j], p, Ypre[j]), p);
//       }



//       ell = MulModPrecon(ell, rd, p, pd);  // why ell != u?
//       ell = MulModPrecon(ell, C.P_r[i+i0], p, C.P_pre[i+i0]);
//       ell = MulModPrecon(ell, C.Qab_r[k+j0-i-i0+n-1], p, C.Qab_pre[k+j0-i-i0+n-1]);

//       u = MulModPrecon(u, rd, p, pd); // save one product here ?
//       u = MulModPrecon(u, C.P_r[k+i0], p, C.P_pre[k+i0]);
//       u = MulModPrecon(u, C.Qab_r[i+j0-k-i0+n-1], p, C.Qab_pre[i+j0-k-i0+n-1]);

//       for (long j = 0; j < alpha; j++){
//       	tmpYi[j] = AddMod(tmpYi[j], MulModPrecon(ell, tmpYk[j], p, Ypre[j]), p);
//       	tmpZi[j] = AddMod(tmpZi[j], MulModPrecon(u, tmpZk[j], p, Zpre[j]), p);
//       }

//       tmpYi += alpha;
//       tmpZi += alpha;
//     }

//     tmpYk += alpha;
//     tmpZk += alpha;
//   }

//   delete[] Ypre;
//   delete[] Zpre;
//   return N;
// }

/*------------------------------------------------------------*/
/* inverts a special-Cauchy-like matrix                       */
/* assumes generic rank profile                               */
/* the generators are Yp, Zp                                  */
/* alpha = displacement rank                                  */
/* overwrites Yp, Zp                                          */
/* i0, j0 are the initial indices in i, j resp.               */
/* i1, j1 are the end indices in i, j resp.                   */  
/*------------------------------------------------------------*/
static long invert_raw(long *Yp, long *Zp, const long alpha,
		       const cauchy_geometric_special& C,
		       const long i0, const long i1, const long j0, const long j1){


  // P[i+i0]*Qab[j+j0-i-i0+n-1] = 1/(a*power(q,i+i0)-b*power(q,j+j0))
  // P[i+i0]*Qa[j+j0-i-i0+n-1] = 1/(a*power(q,i+i0)-a*power(q,j+j0)) = 1/(x[i]-x[j])
  // P[i+i0]*Qb[j+j0-i-i0+n-1] = 1/(b*power(q,i+i0)-b*power(q,j+j0)) = 1/(y[i]-y[j])

  const long n = C.n;
  const long p = zz_p::modulus();
  mulmod_precon_t *Ypre = new mulmod_precon_t[alpha];
  mulmod_precon_t *Zpre = new mulmod_precon_t[alpha];
  const mulmod_t pinv = zz_p::ModulusInverse();

  long *tmpYk = Yp;
  long *tmpZk = Zp;

  //  cout << "enter with: " << i0 << " " << i1 << " " << j0 << " " << j1 << endl;
  long N = min(i1-i0, j1-j0);

  for (long k = 0; k < N; k++){
    long *tmpYi = Yp;
    long *tmpZi = Zp;

    for (long j = 0; j < alpha; j++){
      Ypre[j] = PrepMulModPrecon(tmpYk[j], p, pinv);
      Zpre[j] = PrepMulModPrecon(tmpZk[j], p, pinv);
    }

    // computes the top-left entry in the current submatrix
    long rd = 0;
    for (long j = 0; j < alpha; j++)
      rd = AddMod(rd, MulMod(tmpZk[j], tmpYk[j], p, pinv), p);


    // early exit if rank(M) < n
    if (rd == 0){

      //      cout << "with zero pivot with k=" << k<< endl;

      long *tmpY = tmpYk;
      for (long a = k; a < i1-i0; a++){
    	long *tmpZ = tmpZk;
    	for (long b = k; b < j1-j0; b++){
    	  long tmp = 0;
    	  for (long j = 0; j < alpha; j++)
    	    tmp = AddMod(tmp, MulMod(tmpZ[j], tmpY[j], p, pinv), p);
    	  if (tmp != 0) {
	    //	    cout << "failed " << a << " " << b << endl;
    	    delete[] Ypre;
    	    delete[] Zpre;
    	    return -1;  // not generic rank profile
    	  }
	  //	  cout << "passed " << a << " " << b << endl;
    	  tmpZ += alpha;
    	}
    	tmpY += alpha;
      }
      delete[] Ypre;
      delete[] Zpre;
      return k; // OK, rank = k
    }

    rd = MulModPrecon(rd, C.P_r[k+i0], p, C.P_pre[k+i0]);
    rd = MulModPrecon(rd, C.Qab_r[j0-i0+n-1], p, C.Qab_pre[j0-i0+n-1]);

    zz_p d(rd, INIT_LOOP_HOLE);
    zz_p id = 1/d;
    zz_p mid = -id;

    const long rid = rep(id);
    const long rmid = rep(mid);
    mulmod_precon_t pd = PrepMulModPrecon(rd, p, pinv);
    mulmod_precon_t prid = PrepMulModPrecon(rid, p, pinv);

    for (long i = 0; i < k; i++){

      long ell = 0, u = 0;
      for (long j = 0; j < alpha; j++){
      	ell = AddMod(ell, MulModPrecon(tmpYi[j], tmpZk[j], p, Zpre[j]), p);
      	u = AddMod(u, MulModPrecon(tmpZi[j], tmpYk[j], p, Ypre[j]), p);
      }
      // not so great for small alpha
      // ell = InnerProd_LL(tmpYi, tmpZk, alpha);   
      // u = InnerProd_LL(tmpZi, tmpYk, alpha);

      ell = MulModPrecon(ell, rid, p, prid);
      ell = MulModPrecon(ell, C.P_r[i+i0], p, C.P_pre[i+i0]);
      ell = MulModPrecon(ell, C.Qb_r[k+j0-i-i0+n-1], p, C.Qb_pre[k+j0-i-i0+n-1]);

      u = MulModPrecon(u, rid, p, prid);
      u = MulModPrecon(u, C.P_r[i+i0], p, C.P_pre[i+i0]);
      u = MulModPrecon(u, C.Qa_r[k+j0-i-i0+n-1], p, C.Qa_pre[k+j0-i-i0+n-1]);

      for (long j = 0; j < alpha; j++){
      	tmpYi[j] = SubMod(tmpYi[j], MulModPrecon(ell, tmpYk[j], p, Ypre[j]), p);
      	tmpZi[j] = AddMod(tmpZi[j], MulModPrecon(u, tmpZk[j], p, Zpre[j]), p);
      }

      tmpYi += alpha;
      tmpZi += alpha;
    }

    for (long j = 0; j < alpha; j++){
      tmpYk[j] = MulModPrecon(rmid, tmpYk[j], p, Ypre[j]);
      Ypre[j] = PrepMulModPrecon(tmpYk[j], p, pinv);
      tmpZk[j] = MulModPrecon(rid, tmpZk[j], p, Zpre[j]);
      Zpre[j] = PrepMulModPrecon(tmpZk[j], p, pinv);
    }

    tmpYi += alpha;
    tmpZi += alpha;

    for (long i = k+1; i < (i1-i0); i++){

      long ell = 0;
      for (long j = 0; j < alpha; j++)
	ell = AddMod(ell, MulModPrecon(tmpYi[j], tmpZk[j], p, Zpre[j]), p);

      ell = MulModPrecon(ell, rd, p, pd);  // why ell != u?
      ell = MulModPrecon(ell, C.P_r[i+i0], p, C.P_pre[i+i0]);
      ell = MulModPrecon(ell, C.Qab_r[k+j0-i-i0+n-1], p, C.Qab_pre[k+j0-i-i0+n-1]);

      for (long j = 0; j < alpha; j++)
      	tmpYi[j] = AddMod(tmpYi[j], MulModPrecon(ell, tmpYk[j], p, Ypre[j]), p);

      tmpYi += alpha;
    }

    for (long i = k+1; i < (j1-j0); i++){

      long u = 0;
      for (long j = 0; j < alpha; j++)
        u = AddMod(u, MulModPrecon(tmpZi[j], tmpYk[j], p, Ypre[j]), p);

      u = MulModPrecon(u, rd, p, pd); // save one product here ?
      u = MulModPrecon(u, C.P_r[k+i0], p, C.P_pre[k+i0]);
      u = MulModPrecon(u, C.Qab_r[i+j0-k-i0+n-1], p, C.Qab_pre[i+j0-k-i0+n-1]);

      for (long j = 0; j < alpha; j++)
      	tmpZi[j] = AddMod(tmpZi[j], MulModPrecon(u, tmpZk[j], p, Zpre[j]), p);

      tmpZi += alpha;
    }

    tmpYk += alpha;
    tmpZk += alpha;
  }

  delete[] Ypre;
  delete[] Zpre;
  return N;
}

/*---------------------------------------------------*/
/* iCL = CL^(-1), represented for the dual operator  */
/* O(alpha n^2) algorithm                            */
/*---------------------------------------------------*/
long invert(cauchy_like_geometric_special& Cinv,
	    const cauchy_like_geometric_special& CL){

  long r = -1;

  long n = CL.NumRows();
  long m = CL.NumCols();
  long alpha = CL.Alpha();

  long *Yp = new long[n*alpha];
  for (long i = 0; i < n*alpha; i++)
    Yp[i] = CL.G[i];

  long *Zp = new long[m*alpha];
  for (long i = 0; i < m*alpha; i++)
    Zp[i] = CL.H[i];

  r = invert_raw(Yp, Zp, alpha, CL.C, 0, n, 0, m);

  if (r == -1){
    delete[] Yp;
    delete[] Zp;
    return -1;
  }

  build(Cinv.C, CL.C.b, CL.C.a, CL.C.q, r, r);
  Cinv.n = r;
  Cinv.m = r;
  Cinv.alpha = alpha;
  Cinv.G.SetLength(r*alpha);
  Cinv.H.SetLength(r*alpha);

  for (long i = 0; i < r*alpha; i++){
    Cinv.G[i] = Yp[i];
    Cinv.H[i] = Zp[i];
  }
  
  delete[] Yp;
  delete[] Zp;
  return r;
}



/*---------------------------------------------------*/
/* iCL = CL^(-1), represented for the dual operator  */
/* O(alpha^2 M(n)) algorithm                         */
/*---------------------------------------------------*/
long invert_fast_raw(long *Yp, long *Zp, const long alpha,
		     const cauchy_geometric_special& C,
		     const long i0, const long i1, const long j0, const long j1, const long thresh){
  
  const long N = min(i1-i0, j1-j0);
  if (N < thresh)
    return invert_raw(Yp, Zp, alpha, C, i0, i1, j0, j1);

  const long n = C.n;
  const long N1 = N/2;
  const long N2 = N-N1;

  const long alphaN = N*alpha;
  const long alphaN1 = N1*alpha;

  long *G1 = new long[alpha*(i1-i0)];
  long *H1 = new long[alpha*(j1-j0)];
  long *G2 = G1 + alphaN1;
  long *H2 = H1 + alphaN1;
  long *iSG = Yp + alphaN1;
  long *iSH = Zp + alphaN1;

  // copies to BOTH (G1,H1) and (G2,H2)
  for (long i = 0; i < alphaN; i++){
    G1[i] = Yp[i];
    H1[i] = Zp[i];
  }
  
  long r = invert_fast_raw(G1, H1, alpha, C, i0, i0+N1, j0, j0+N1, thresh);

  // case 1: we found that we do not have grp. abort.
  if (r == -1){
    delete[] G1;
    delete[] H1;
    return -1;
  }

  // case 2: we found the rank of the top-left, which was grp. check all the rest.
  if (r < N1){
    G2 = G1 + alpha*r;
    H2 = H1 + alpha*r;
    iSG = Yp + alpha*r;
    iSH = Zp + alpha*r;

    // attempt fix
    mat_addmul(G2, G1, iSG, Zp, C.P_pre, C.Qab_pre, C.P_r, C.Qab_r, i0+r, n+j0-i1, i1-i0-r, r, alpha);
    mat_submul_t(H2, H1, iSH, Yp, C.P_pre, C.Qab_pre, C.P_r, C.Qab_r, i0, n+j0-i0, r, j1-j0-r, alpha);

    // check if G2.H2^t = 0 by computing random_vec.G2.H2^t*random_vec
    long *tmpG = new long[alpha];
    long *tmpH = new long[alpha];

    for (long j = 0; j < alpha; j++)
      tmpG[j] = tmpH[j] = 0;

    for (long i = 0; i < i1-i0-r; i++){
      long cf = rep(random_zz_p());
      for (long j = 0; j < alpha; j++)
	tmpG[j] = AddMod(tmpG[j], MulMod(G2[i*alpha+j], cf, zz_p::modulus()), zz_p::modulus());
    }

    for (long i = 0; i < j1-j0-r; i++){
      long cf = rep(random_zz_p());
      for (long j = 0; j < alpha; j++)
	tmpH[j] = AddMod(tmpH[j], MulMod(H2[i*alpha+j], cf, zz_p::modulus()), zz_p::modulus());
    }

    long sum = 0;
    for (long k = 0; k < alpha; k++)
      sum = AddMod(sum, MulMod(tmpG[k], tmpH[k], zz_p::modulus()), zz_p::modulus());
    if (sum != 0)
      r = -1;

    delete[] tmpG;
    delete[] tmpH;

    if (r != -1){
      for (long i = 0; i < alpha*r; i++){
    	Yp[i] = G1[i];
    	Zp[i] = H1[i];
      }
    }

    delete[] G1;
    delete[] H1;
    return r;
  }

  // case 3: the top block was invertible.
  // look at the whole matrix
  mat_addmul(G2, G1, iSG, Zp, C.P_pre, C.Qab_pre, C.P_r, C.Qab_r, i0+N1, n+j0-i1, i1-i0-N1, N1, alpha);
  mat_submul_t(H2, H1, iSH, Yp, C.P_pre, C.Qab_pre, C.P_r, C.Qab_r, i0, n+j0-i0, N1, j1-j0-N1, alpha);

  // copies to BOTH (Yp,Zp) and (iSG,iSH)
  for (long i = 0; i < alphaN; i++){
    Yp[i] = G1[i];
    Zp[i] = H1[i];
  }

  // looks at the whole matrix
  r = invert_fast_raw(iSG, iSH, alpha, C, i0+N1, i1, j0+N1, j1, thresh);

  if (r == -1){
    delete[] G1;
    delete[] H1;
    return -1;
  }

  mat_addmul(Yp, iSG, G1, H2, C.P_pre, C.Qb_pre, C.P_r, C.Qb_r, i0, n+j0-i0, N1, N2, alpha);
  mat_addmul(Zp, iSH, H1, G2, C.P_pre, C.Qa_pre, C.P_r, C.Qa_r, i0, n+j0-i0, N1, N2, alpha);

  delete[] G1;
  delete[] H1;
  return N1+r;
}

/*---------------------------------------------------*/
/* iCL = CL^(-1), represented for the dual operator  */
/* O(alpha^2 M(n)log(n)) algorithm                   */
/*---------------------------------------------------*/
long invert_fast(cauchy_like_geometric_special& Cinv,
		 const cauchy_like_geometric_special& CL, const long thresh){

  long n = CL.NumRows();
  long m = CL.NumCols();
  long alpha = CL.Alpha();

  long do_thresh;
  if (thresh == -1)
    do_thresh = threshold(alpha);
  else
    do_thresh = thresh;

  long *Yp = new long[n*alpha];
  for (long i = 0; i < n*alpha; i++)
    Yp[i] = CL.G[i];

  long *Zp = new long[m*alpha];
  for (long i = 0; i < m*alpha; i++)
    Zp[i] = CL.H[i];

  long r = invert_fast_raw(Yp, Zp, alpha, CL.C, 0, n, 0, m, do_thresh);

  if (r == -1){
    delete[] Yp;
    delete[] Zp;
    return -1;
  }

  build(Cinv.C, CL.C.b, CL.C.a, CL.C.q, r, r);
  Cinv.n = r;
  Cinv.m = r;
  Cinv.alpha = alpha;
  Cinv.G.SetLength(r*alpha);
  Cinv.H.SetLength(r*alpha);

  for (long i = 0; i < r*alpha; i++){
    Cinv.G[i] = Yp[i];
    Cinv.H[i] = Zp[i];
  }

  delete[] Yp;
  delete[] Zp;
  return r;
}
