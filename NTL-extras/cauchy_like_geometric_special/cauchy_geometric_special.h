#ifndef __CAUCHY_GEOMETRIC_SPECIAL_H
#define __CAUCHY_GEOMETRIC_SPECIAL_H

#include <NTL/lzz_pX.h>
#include <NTL/mat_lzz_p.h>

#define THRESHOLD 100

NTL_CLIENT

/*----------------------------------------------------*/
/*----------------------------------------------------*/
/* special Cauchy matrices on geometric progressions  */
/*----------------------------------------------------*/
/*----------------------------------------------------*/

class cauchy_geometric_special{

public:
  zz_p a, b, q;
  long n, m;   // dimension is n x m
  // the first points are a, aq, aq^2, ..., aq^{n-1}
  // the second points are b, bq, bq^2, ..., bq^{m-1}

  // P[i]*Qab[j-i+n-1] = 1/(a*power(q,i)-b*power(q,j))
  // P[i]*Qa[j-i+n-1] = 1/(a*power(q,i)-a*power(q,j)) = 1/(x[i]-x[j])
  // P[i]*Qb[j-i+n-1] = 1/(b*power(q,i)-b*power(q,j)) = 1/(y[i]-y[j])
  // lengths n, m+n-1, n+m-1, n+m-1
  Vec<mulmod_precon_t> P_pre, Qab_pre, Qa_pre, Qb_pre;
  Vec<long> P_r, Qab_r, Qa_r, Qb_r;

  /*----------------------------------------------------*/
  /* constructor for cauchy_geometric, do nothing       */
  /*----------------------------------------------------*/
  cauchy_geometric_special(){}

  /*----------------------------------------------------*/
  /* constructor for cauchy_geometric: computes all info*/
  /* (calls build below)                                */
  /*----------------------------------------------------*/
  cauchy_geometric_special(const zz_p& a, const zz_p& b, const zz_p& q, const long n, const long m);
};

/*----------------------------------------------------*/
/* prepares the vecs for the whole C                  */
/*----------------------------------------------------*/
void build(cauchy_geometric_special& C, 
	   const zz_p& a, const zz_p& b, const zz_p& q, long n, long m);


/*----------------------------------------------------*/
/*----------------------------------------------------*/
/* special Cauchy-like matrices                       */
/*----------------------------------------------------*/
/*----------------------------------------------------*/

class cauchy_like_geometric_special{

 public:
  long n, m;   // dimension is n x m
  long alpha; // displacement rank

  Vec<long> G, H;

  cauchy_geometric_special C;

  long NumRows() const {
    return n;
  }

  long NumCols() const {
    return m;
  }

  long Alpha() const {
    return alpha;
  }

  /*----------------------------------------------------*/
  /* constructor for cauchy_like_geometric, do nothing  */
  /*----------------------------------------------------*/
  cauchy_like_geometric_special(){}

  /*----------------------------------------------------*/
  /* constructor: copies arguments                      */
  /*----------------------------------------------------*/
  cauchy_like_geometric_special(const mat_zz_p& G_in, const mat_zz_p& H_in, const cauchy_geometric_special& C){
    this->C = C;
    n = G_in.NumRows();
    m = H_in.NumRows();
    alpha = G_in.NumCols();

    G.SetLength(n*alpha);
    H.SetLength(m*alpha);

    long idx;
    idx = 0;
    for (long i = 0; i < n; i++)
      for (long j = 0; j < alpha; j++)
	G[idx++] = rep(G_in[i][j]);
    idx = 0;
    for (long i = 0; i < m; i++)
      for (long j = 0; j < alpha; j++)
	H[idx++] = rep(H_in[i][j]);
  }

  ~cauchy_like_geometric_special(){}
};

/*----------------------------------------------------*/
/* A is such that DA - AD' = G.H^t                    */
/*----------------------------------------------------*/
void to_dense(mat_zz_p& A, const cauchy_like_geometric_special& CL);

/*---------------------------------------------------*/
/* iCL = CL^(-1), represented for the dual operator  */
/* O(alpha n^2) algorithm                            */
/*---------------------------------------------------*/
long invert(cauchy_like_geometric_special& Cinv,
	    const cauchy_like_geometric_special& CL);


/*---------------------------------------------------*/
/* iCL = CL^(-1), represented for the dual operator  */
/* O(alpha^2 M(n)log(n)) algorithm                   */
/*---------------------------------------------------*/
long invert_fast(cauchy_like_geometric_special& Cinv,
		 const cauchy_like_geometric_special& CL, const long thresh = -1);

/*---------------------------------------------------*/
/* computes M*input                                  */
/*---------------------------------------------------*/
void mul_right(Vec<zz_p>& res, const cauchy_like_geometric_special& M, const Vec<zz_p>& input);

/*---------------------------------------------------*/
/* computes M^t*input                                */
/*---------------------------------------------------*/
void mul_left(Vec<zz_p>& res, const cauchy_like_geometric_special& M, const Vec<zz_p>& input);

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
		long a); // displacement rank

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
		long a); // displacement rank

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
		  long a);  // displacement rank

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
		  long a);  // displacement rank

#endif
