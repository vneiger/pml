#ifndef MAT_LZZ_PX_LINEARIZATION__H
#define MAT_LZZ_PX_LINEARIZATION__H

/** \brief Experimental code for partial linearization.
 *
 * \file mat_lzz_pX_linearization.h
 * \author Seung Gyu Hyun, Vincent Neiger, Eric Schost
 * \version 0.1
 * \date 2019-02-01
 *
 * \todo Not documented yet. Will be when the functions become more stable.
 */

#include <numeric> // for 'accumulate'
#include "mat_lzz_pX_forms.h" // for VecLong and col_degree()

NTL_CLIENT


/*------------------------------------------------------------*/
/* horizontal join                                            */
/* requires a.NumRows() == b.NumRows()                        */
/*------------------------------------------------------------*/
void horizontal_join(Mat<zz_pX>& c, const Mat<zz_pX>& a, const Mat<zz_pX>& b);

inline Mat<zz_pX> horizontal_join(const Mat<zz_pX>& a, const Mat<zz_pX>& b)
{ Mat<zz_pX> c; horizontal_join(c, a, b); return c; }

// TODO vertical join
// TODO vertical/horizontal splits (then update kernel basis)

/*------------------------------------------------------------*/
/* collapses s consecutive columns of a into one column of c  */
/* let t=a.NumCols(). For i=0..t/s-1, the i-th column of c is */
/* a[i*s] + x^d a[i*s+1] + ... + x^{(s-1)*d} a[i*s+s-1)]      */
/* requires that s divides t exactly                          */
/*------------------------------------------------------------*/
void collapse_consecutive_columns(Mat<zz_pX>& c, const Mat<zz_pX>& a, long d, long s);

inline Mat<zz_pX> collapse_consecutive_columns(const Mat<zz_pX>& a, long d, long s)
{ Mat<zz_pX> c; collapse_consecutive_columns(c, a, d, s); return c; }

/*------------------------------------------------------------*/
/* collapses columns with stepsize s of a into a column of c  */
/* let t=a.NumCols(). For i=0..s-1, the i-th column of c is   */
/* a[i] + x^d a[i+s] + ... + x^{(t/s-1)*d} a[i+(t/s-1)*s)]    */
/* requires that s divides t exactly                          */
/*------------------------------------------------------------*/
void collapse_nonconsecutive_columns(Mat<zz_pX>& c, const Mat<zz_pX>& a, long d, long s);

inline Mat<zz_pX> collapse_nonconsecutive_columns(const Mat<zz_pX>& a, long d, long s)
{ Mat<zz_pX> c; collapse_nonconsecutive_columns(c, a, d, s); return c; }




// TODO all the functions below are experimental and have not well
// been tested; they are not used in other functions for the moment



/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* COMPLETE LINEARIZATION                                     */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

// Column linearization
// transforms matrix into constant matrix of coefficients
// TODO improve description
// TODO not implemented yet
VecLong column_linearization(
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
VecLong column_partial_linearization(
                                     Mat<zz_pX> & parlin, 
                                     const Mat<zz_pX> & pmat, 
                                     const VecLong & parlin_degree,
                                     const VecLong & target_degree
                                    );

// same, using cdeg
// TODO improve description
inline VecLong column_partial_linearization_cdeg(
                                                 Mat<zz_pX> &parlin, 
                                                 const Mat<zz_pX> & pmat, 
                                                 const VecLong & parlin_degree
                                                )
{
    VecLong cdeg(pmat.NumCols());
    col_degree(cdeg, pmat);
    return column_partial_linearization(parlin, pmat, parlin_degree, cdeg);
}

// Column partial linearization, uniform parlin degree
// returns vector of linearization parameters
// target_degree default: set as the col_degree of pmat
// TODO improve description
inline VecLong column_partial_linearization(
                                            Mat<zz_pX> &parlin, 
                                            const Mat<zz_pX> & pmat, 
                                            const long parlin_degree,
                                            const VecLong & target_degree
                                           )
{
    VecLong parlin_degrees(pmat.NumCols(), parlin_degree);
    return column_partial_linearization(parlin, pmat, parlin_degrees, target_degree);
}

// same using cdeg
// TODO improve description
inline VecLong column_partial_linearization_cdeg(
                                                 Mat<zz_pX> &parlin, 
                                                 const Mat<zz_pX> & pmat, 
                                                 const long parlin_degree
                                                )
{
    VecLong cdeg(pmat.NumCols());
    col_degree(cdeg, pmat);
    VecLong parlin_degrees(pmat.NumCols(), parlin_degree);
    return column_partial_linearization(parlin, pmat, parlin_degrees, cdeg);
}


// Column partial linearization, uniform target and parlin degrees
// returns vector of linearization parameters
// TODO improve description
inline VecLong column_partial_linearization(
                                            Mat<zz_pX> &parlin, 
                                            const Mat<zz_pX> & pmat, 
                                            const long parlin_degree,
                                            const long target_degree
                                           )
{
    VecLong parlin_degrees(pmat.NumCols(), parlin_degree);
    VecLong target_degrees(pmat.NumCols(), target_degree);
    return column_partial_linearization(parlin, pmat, parlin_degrees, target_degrees);
}

// Column partial linearization, default degrees
// returns vector of linearization parameters
// TODO improve description
inline VecLong column_partial_linearization(
                                            Mat<zz_pX> &parlin, 
                                            const Mat<zz_pX> & pmat, 
                                            const VecLong & target_degree
                                           )
{
    long parlin_degree = std::accumulate(target_degree.begin(), target_degree.end(), 0);
    parlin_degree = 1 + parlin_degree / target_degree.size();
    VecLong parlin_degrees(pmat.NumCols(), parlin_degree);
    return column_partial_linearization(parlin, pmat, parlin_degrees, target_degree);
}

// same with cdeg
// TODO improve description
inline VecLong column_partial_linearization_cdeg(
                                                 Mat<zz_pX> &parlin, 
                                                 const Mat<zz_pX> & pmat
                                                )
{
    VecLong cdeg(pmat.NumCols());
    col_degree(cdeg, pmat);
    long parlin_degree = std::accumulate(cdeg.begin(), cdeg.end(), 0);
    parlin_degree = 1 + parlin_degree / cdeg.size();
    VecLong parlin_degrees(pmat.NumCols(), parlin_degree);
    return column_partial_linearization(parlin, pmat, parlin_degrees, cdeg);
}

VecLong column_partial_linearization(
                                     Mat<zz_pX> & parlin, 
                                     const Mat<zz_pX> & pmat, 
                                     const VecLong & parlin_degree,
                                     const VecLong & target_degree,
                                     const long dinf
                                    );

/*------------------------------------------------------------*/
/* Basic row partial linearizations:                          */
/* expand rows of the matrix according to a                   */
/* specified degree profile                                   */
/*------------------------------------------------------------*/
// TODO


// TODO: doc for the functions below (parlin_multiply / middleprod)
// here, col_degree must be an (non-strict) upper bound on the column degree of b
void right_parlin_multiply(
                           Mat<zz_pX> &c,
                           const Mat<zz_pX> &a,
                           const Mat<zz_pX> &b,
                           const long parlin_degree,
                           const VecLong & col_degree
                          );

// column degree of b not provided, compute it
inline void right_parlin_multiply(
                                  Mat<zz_pX> &c,
                                  const Mat<zz_pX> &a,
                                  const Mat<zz_pX> &b,
                                  const long parlin_degree
                                 )
{
    VecLong cdeg(b.NumCols());
    col_degree(cdeg, b);
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
    VecLong target_degrees(b.NumCols(), degree);
    right_parlin_multiply(c, a, b, parlin_degree, target_degrees);
}

// uses deg(a) as linearization degree parameter
inline void right_parlin_multiply(
                                  Mat<zz_pX> &c,
                                  const Mat<zz_pX> &a,
                                  const Mat<zz_pX> &b
                                 )
{
    VecLong cdeg(b.NumCols());
    col_degree(cdeg, b);
    right_parlin_multiply(c, a, b, deg(a), cdeg);
}

/* c = trunc( trunc(a, dA+1)*b div x^dA, dB+1 )           */
void right_parlin_middle_product(
                                 Mat<zz_pX> &c,
                                 const Mat<zz_pX> &a,
                                 const Mat<zz_pX> &b,
                                 const long parlin_degree,
                                 const VecLong & col_degree,
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
    VecLong col_degree(b.NumCols(), degree);
    right_parlin_middle_product(c, a, b, parlin_degree, col_degree, dA, dB);
}

#endif // MAT_LZZ_PX_LINEARIZATION__H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
