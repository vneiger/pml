#include <NTL/BasicThreadPool.h>
#include <NTL/lzz_pX.h>
#include <NTL/matrix.h>
#include <NTL/vector.h>
#include <iomanip>

#include "util.h"
#include "mat_lzz_pX_extra.h"

NTL_CLIENT

#ifndef __EXAMPLE_H__

// zero matrix
Mat<zz_pX> mat_zero_square; // square zero matrix
Mat<zz_pX> mat_zero_64x63;
Mat<zz_pX> mat_zero_64x48;
Mat<zz_pX> mat_zero_64x32;
Mat<zz_pX> mat_zero_64x16;
Mat<zz_pX> mat_zero_64x1;

// random matrices
Mat<zz_pX> mat_square; // square matrix 64x64
Mat<zz_pX> mat_rect_64x63;
Mat<zz_pX> mat_rect_64x48;
Mat<zz_pX> mat_rect_64x32;
Mat<zz_pX> mat_rect_64x16;
Mat<zz_pX> mat_rect_64x1;

// rank deficient (square 64x64) matrices
Mat<zz_pX> mat_square_rank1;
Mat<zz_pX> mat_square_rank16;
Mat<zz_pX> mat_square_rank32;
Mat<zz_pX> mat_square_rank48;
Mat<zz_pX> mat_square_rank63;

static void rand_rows (Mat<zz_pX> &mat, const long r, const long c, const long rows, const long deg)
{
  mat = Mat<zz_pX>(r,c);
  for (long i = 0; i < rows; i++)
    for (long j = 0; j < c; j++)
      random(mat[i][j],deg+1);
}

// call these to setup matrices
void mat_setup_zero()
{
   mat_zero_square = Mat<zz_pX>(64,64);
   mat_zero_64x63 = Mat<zz_pX>(64,63);
   mat_zero_64x48 = Mat<zz_pX>(64,48);
   mat_zero_64x32 = Mat<zz_pX>(64,32);
   mat_zero_64x16 = Mat<zz_pX>(64,16);
   mat_zero_64x1 = Mat<zz_pX>(64,1);
}

void mat_setup_rand(long deg = 20)
{
    mat_square = Mat<zz_pX>(64,64); // square matrix 64x64
    random(mat_square,64,64,deg+1);
    mat_rect_64x63 = Mat<zz_pX>(64,63);
    random(mat_rect_64x63,64,63,deg+1);
    mat_rect_64x48 = Mat<zz_pX>(64,48);
    random(mat_rect_64x48,64,48,deg+1);
    mat_rect_64x32 = Mat<zz_pX>(64,32);
    random(mat_rect_64x32,64,32,deg+1);
    mat_rect_64x16 = Mat<zz_pX>(64,16);
    random(mat_rect_64x16,64,16,deg+1);
    mat_rect_64x1 = Mat<zz_pX>(64,1);
    random(mat_rect_64x1,64,1,deg+1);
}

void mat_setup_rank_deficient(long deg=20)
{
    rand_rows(mat_square_rank1,64,64,1,20);
    rand_rows(mat_square_rank16,64,64,16,20);
    rand_rows(mat_square_rank32,64,64,32,20);
    rand_rows(mat_square_rank48,64,64,48,20);
    rand_rows(mat_square_rank63,64,64,63,20);
}

#define __EXAMPLE_H__




















// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
