#ifndef __MOSAIC_HANKEL_H
#define __MOSAIC_HANKEL_H

#include <NTL/lzz_pX.h>
#include <NTL/mat_lzz_p.h>

#include "lzz_pX_CRT.h"
#include "cauchy_geometric_special.h"

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
class hankel{

  // n rows, m columns
  long n, m;
  Vec<zz_p> data;

 public:

  Vec<zz_p> data_rev;
  fftRep fft_data;

  hankel(){
    data.SetLength(0);
    n = m = 0;
  }

  /*----------------------------------------------------*/
  /* input vector is as showed above                    */
  /*----------------------------------------------------*/
  hankel(Vec<zz_p>& input, long rows, long cols){
    n = rows;
    m = cols;
    data = input;
    data_rev.SetLength(n+m-1);
    zz_pX data_X;
    data_X.rep.SetLength(n+m-1);

    for (long i = 0;i < n+m-1; i++){
      data_rev[i] = input[n+m-2-i];
      data_X.rep[i] = data_rev[i];
    }    
    data_X.normalize();

    long K = NextPowerOfTwo(n+m-1);
    fft_data = fftRep(INIT_SIZE, K);
    TofftRep(fft_data, data_X, K);
  }

  /*----------------------------------------------------*/
  /* getters                                            */
  /*----------------------------------------------------*/
  long NumCols() const{
    return m;
  }
  long NumRows() const{
    return n;
  }
  const zz_p& operator ()(long i, long j) const{
    return data[m+n-2-i-j];
  }

};

/*----------------------------------------------------*/
/* turns M into a dense matrix                        */
/*----------------------------------------------------*/
void to_dense(Mat<zz_p>& Mdense, const hankel& M);

/*----------------------------------------------------*/
/* low-level right multiplication                     */
/* assumes input and res have the right size          */
/*----------------------------------------------------*/
void mul_right(Vec<zz_p>& res, const hankel& M, const Vec<zz_p>& input);

/*----------------------------------------------------*/
/* low-level left multiplication                      */
/* assumes input and res have the right size          */
/*----------------------------------------------------*/
void mul_left(Vec<zz_p>& res, const hankel& M, const Vec<zz_p>& input);




/*----------------------------------------------------*/
/*----------------------------------------------------*/
/* Mosaic Hankel matrices:                            */
/* block matrix where each block is Hankel            */
/*----------------------------------------------------*/
/*----------------------------------------------------*/
class mosaic_hankel{
  long n, m; // numbers of rows / colums
  long nb, mb; // number of blocks in rows / columns

 public:

  // data[i][j] is the block at i-th row, j-th column
  Vec< Vec<hankel> > data;

  /*----------------------------------------------------*/
  /* dummy constructor                                  */
  /*----------------------------------------------------*/
  mosaic_hankel(){
    data.SetLength(0);
    n = m = nb = mb = 0;
  }

  /*----------------------------------------------------*/
  /* copies all data                                    */
  /*----------------------------------------------------*/
  mosaic_hankel(Vec< Vec<hankel> > init){
    data = init;
    nb = init.length();
    mb = init[0].length();

    n = 0;
    m = 0;

    for(long i = 0; i < nb; i++)
      n += init[i][0].NumRows();
    for (long j = 0; j < mb; j++)
      m += init[0][j].NumCols();
  }

  /*----------------------------------------------------*/
  /* getters                                            */
  /*----------------------------------------------------*/
  long NumRows() const{
    return n;
  }
  long NumCols() const{
    return m;
  }
  long NumBlockRows() const{
    return nb;
  }
  long NumBlockCols() const{
    return mb;
  }
  long NumRows_of_block(long i) const{
    return data[i][0].NumRows();
  }
  long NumCols_of_block(long i) const{
    return data[0][i].NumCols();
  }

  const zz_p& operator ()(long i, long j) const{

    if (i >= n || j >= m)
      Error("matrix indices out of bounds\n");
    long idx, jdx;

    idx = 0;
    while (i >= data[idx][0].NumRows()){
      i -= data[idx][0].NumRows();
      idx++;
    }

    jdx = 0;
    while (j >= data[0][jdx].NumCols()){
      j -= data[0][jdx].NumCols();
      jdx++;
    }
    return data[idx][jdx](i, j);
  }

};

/*----------------------------------------------------*/
/* turns M into a dense matrix                        */
/*----------------------------------------------------*/
void to_dense(Mat<zz_p>& Mdense, const mosaic_hankel& M);

/*----------------------------------------------------*/
/* right multiplication                               */
/*----------------------------------------------------*/
void mul_right(Vec<zz_p>& res, const mosaic_hankel& M, const Vec<zz_p>& input);

/*----------------------------------------------------*/
/* left multiplication                                */
/*----------------------------------------------------*/
void mul_left(Vec<zz_p>& res, const mosaic_hankel& M, const Vec<zz_p>& input);

/*----------------------------------------------------*/
/* access to particular rows and columns              */
/*----------------------------------------------------*/
void first_column_of_block(Vec<zz_p>& res, long i, const mosaic_hankel& M);
void last_column_of_block(Vec<zz_p>& res, long i, const mosaic_hankel& M);
void first_row_of_block(Vec<zz_p>& res, long i, const mosaic_hankel& M);
void last_row_of_block(Vec<zz_p>& res, long i, const mosaic_hankel& M);

/*----------------------------------------------------*/
/* G, H such that Z1 M - Z0^t M = G H^t               */
/*----------------------------------------------------*/
void generators(Mat<zz_p>& G, Mat<zz_p>& H, const mosaic_hankel& M); 

/*------------------------------------------------------------------*/
/* preconditions M                                                  */
/* builds the matrix CL = (D_e X_int) M (D_f Y_int)^t, where:       */
/* - X_int, Y_int are geometric interpolation                       */
/* - D_e, D_f are diagonal matrix built on vectors e and f          */
/* - CL is cauchy-like special                                      */
/* - CL is expected to have generic rank profile                    */
/*------------------------------------------------------------------*/
void to_cauchy_grp(cauchy_like_geometric_special& CL, 
		   zz_pX_Multipoint_Geometric& X_int, zz_pX_Multipoint_Geometric& Y_int,
		   Vec<zz_p> &e, Vec<zz_p> &f,
		   const mosaic_hankel& M);

#endif
