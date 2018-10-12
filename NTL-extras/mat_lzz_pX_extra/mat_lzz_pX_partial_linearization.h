#ifndef MAT_LZZ_PX_PARTIAL_LINEARIZATION__H
#define MAT_LZZ_PX_PARTIAL_LINEARIZATION__H

#include <NTL/matrix.h>
#include <NTL/lzz_pX.h>
#include <iostream>
#include <vector>
#include <numeric> // for 'accumulate'

#include "mat_lzz_pX_extra.h"

NTL_CLIENT



/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* COMPLETE LINEARIZATION                                     */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

// Column linearization
// transforms matrix into constant matrix of coefficients
// TODO improve description
// TODO not implemented yet
std::vector<long> column_linearization(
                                       Mat<zz_p> & lin, 
                                       const Mat<zz_pX> & pmat
                                      );



/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* PARTIAL LINEARIZATION                                      */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* Basic linearizations:                                      */
/* expand columns of the matrix according to a                */
/* specified degree profile                                   */
/*------------------------------------------------------------*/

// Column partial linearization
// returns vector of linearization parameters
// parlin_degree and target_degree must have length pmat.NumCols()
// TODO improve description
std::vector<long> column_partial_linearization(
                                               Mat<zz_pX> & parlin, 
                                               const Mat<zz_pX> & pmat, 
                                               const DegVec & parlin_degree,
                                               const DegVec & target_degree
                                              );

// same, using cdeg
// TODO improve description
inline std::vector<long> column_partial_linearization_cdeg(
                                                           Mat<zz_pX> &parlin, 
                                                           const Mat<zz_pX> & pmat, 
                                                           const DegVec & parlin_degree
                                                          )
{
    DegVec cdeg(pmat.NumCols());
    column_degree(cdeg, pmat);
    return column_partial_linearization(parlin, pmat, parlin_degree, cdeg);
}

// Column partial linearization, uniform parlin degree
// returns vector of linearization parameters
// target_degree default: set as the column_degree of pmat
// TODO improve description
inline std::vector<long> column_partial_linearization(
                                                      Mat<zz_pX> &parlin, 
                                                      const Mat<zz_pX> & pmat, 
                                                      const long parlin_degree,
                                                      const DegVec & target_degree
                                                     )
{
    DegVec parlin_degrees(pmat.NumCols(), parlin_degree);
    return column_partial_linearization(parlin, pmat, parlin_degrees, target_degree);
}

// same using cdeg
// TODO improve description
inline std::vector<long> column_partial_linearization_cdeg(
                                                           Mat<zz_pX> &parlin, 
                                                           const Mat<zz_pX> & pmat, 
                                                           const long parlin_degree
                                                          )
{
    DegVec cdeg(pmat.NumCols());
    column_degree(cdeg, pmat);
    DegVec parlin_degrees(pmat.NumCols(), parlin_degree);
    return column_partial_linearization(parlin, pmat, parlin_degrees, cdeg);
}


// Column partial linearization, uniform target and parlin degrees
// returns vector of linearization parameters
// TODO improve description
inline std::vector<long> column_partial_linearization(
                                                      Mat<zz_pX> &parlin, 
                                                      const Mat<zz_pX> & pmat, 
                                                      const long parlin_degree,
                                                      const long target_degree
                                                     )
{
    DegVec parlin_degrees(pmat.NumCols(), parlin_degree);
    DegVec target_degrees(pmat.NumCols(), target_degree);
    return column_partial_linearization(parlin, pmat, parlin_degrees, target_degrees);
}

// Column partial linearization, default degrees
// returns vector of linearization parameters
// TODO improve description
inline std::vector<long> column_partial_linearization(
                                                      Mat<zz_pX> &parlin, 
                                                      const Mat<zz_pX> & pmat, 
                                                      const DegVec & target_degree
                                                     )
{
    long parlin_degree = std::accumulate(target_degree.begin(), target_degree.end(), 0);
    parlin_degree = 1 + parlin_degree / target_degree.size();
    DegVec parlin_degrees(pmat.NumCols(), parlin_degree);
    return column_partial_linearization(parlin, pmat, parlin_degrees, target_degree);
}

// same with cdeg
// TODO improve description
inline std::vector<long> column_partial_linearization_cdeg(
                                                           Mat<zz_pX> &parlin, 
                                                           const Mat<zz_pX> & pmat
                                                          )
{
    DegVec cdeg(pmat.NumCols());
    column_degree(cdeg, pmat);
    long parlin_degree = std::accumulate(cdeg.begin(), cdeg.end(), 0);
    parlin_degree = 1 + parlin_degree / cdeg.size();
    DegVec parlin_degrees(pmat.NumCols(), parlin_degree);
    return column_partial_linearization(parlin, pmat, parlin_degrees, cdeg);
}

std::vector<long> column_partial_linearization(
                                               Mat<zz_pX> & parlin, 
                                               const Mat<zz_pX> & pmat, 
                                               const DegVec & parlin_degree,
                                               const DegVec & target_degree,
                                               const long dinf
                                              );

/*------------------------------------------------------------*/
/* Basic row partial linearizations:                          */
/* expand rows of the matrix according to a                   */
/* specified degree profile                                   */
/*------------------------------------------------------------*/
// TODO


// TODO: doc for the functions below (parlin_multiply / middleprod)
// here, column_degree must be an (non-strict) upper bound on the column degree of b
void right_parlin_multiply(
                           Mat<zz_pX> &c,
                           const Mat<zz_pX> &a,
                           const Mat<zz_pX> &b,
                           const long parlin_degree,
                           const DegVec & column_degree
                          );

// column degree of b not provided, compute it
inline void right_parlin_multiply(
                                  Mat<zz_pX> &c,
                                  const Mat<zz_pX> &a,
                                  const Mat<zz_pX> &b,
                                  const long parlin_degree
                                 )
{
    DegVec cdeg(b.NumCols());
    column_degree(cdeg, b);
    right_parlin_multiply(c, a, b, parlin_degree, cdeg);
}

// global degree provided
inline void right_parlin_multiply(
                                  Mat<zz_pX> &c,
                                  const Mat<zz_pX> &a,
                                  const Mat<zz_pX> &b,
                                  const long parlin_degree,
                                  const long degree
                                 )
{
    DegVec target_degrees(b.NumCols(), degree);
    right_parlin_multiply(c, a, b, parlin_degree, target_degrees);
}

// uses deg(a) as linearization degree parameter
inline void right_parlin_multiply(
                                  Mat<zz_pX> &c,
                                  const Mat<zz_pX> &a,
                                  const Mat<zz_pX> &b
                                 )
{
    DegVec cdeg(b.NumCols());
    column_degree(cdeg, b);
    right_parlin_multiply(c, a, b, deg(a), cdeg);
}

/* c = trunc( trunc(a, dA+1)*b div x^dA, dB+1 )           */
void right_parlin_middle_product(
                                 Mat<zz_pX> &c,
                                 const Mat<zz_pX> &a,
                                 const Mat<zz_pX> &b,
                                 const long parlin_degree,
                                 const DegVec & column_degree,
                                 long dA,
                                 long dB
                                );

inline void right_parlin_middle_product(
                                        Mat<zz_pX> &c,
                                        const Mat<zz_pX> &a,
                                        const Mat<zz_pX> &b,
                                        const long parlin_degree,
                                        const long degree,
                                        long dA,
                                        long dB
                                       )
{
    DegVec column_degree(b.NumCols(), degree);
    right_parlin_middle_product(c, a, b, parlin_degree, column_degree, dA, dB);
}

#endif // MAT_LZZ_PX_PARTIAL_LINEARIZATION__H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
