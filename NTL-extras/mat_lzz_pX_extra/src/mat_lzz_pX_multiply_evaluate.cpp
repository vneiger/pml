#include <NTL/matrix.h>
#include <NTL/mat_lzz_p.h>
#include <NTL/lzz_pX.h>

#include "lzz_p_extra.h"
#include "mat_lzz_pX_extra.h"
#include "lzz_pX_CRT.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* c = a*b                                                    */
/* geometric points                                           */
/*------------------------------------------------------------*/
void multiply_evaluate_geometric(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b)
{
    long dA = deg(a);
    long dB = deg(b);
    long dC = dA+dB;
    long sz = dC+1;

    zz_pX_Multipoint_Geometric ev_geom = get_geometric_points(sz);

    ev_geom.prepare_degree(dA);
    ev_geom.prepare_degree(dB);

    Vec<Mat<zz_p>> valA, valB, valC;
    ev_geom.evaluate_matrix(valA, a);
    ev_geom.evaluate_matrix(valB, b);

    long len = ev_geom.length();
    valC.SetLength(len);
    for (long i = 0; i < len; i++)
    {
        mul(valC[i], valA[i], valB[i]);
    }

    ev_geom.interpolate_matrix(c, valC);
}

/*------------------------------------------------------------*/
/* c = a*b                                                    */
/* FFT points                                                 */
/*------------------------------------------------------------*/
void multiply_evaluate_FFT(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b)
{
  
    long dA = deg(a);
    long dB = deg(b);
    long dC = dA+dB;
    long sz = dC+1;

    zz_pX_Multipoint_FFT ev_FFT = get_FFT_points(sz);

    Vec<Mat<zz_p>> valA, valB, valC;
    ev_FFT.evaluate_matrix(valA, a);
    ev_FFT.evaluate_matrix(valB, b);

    long len = ev_FFT.length();
    valC.SetLength(len);
    for (long i = 0; i < len; i++)
    {
        mul(valC[i], valA[i], valB[i]);
    }

    ev_FFT.interpolate_matrix(c, valC);
}


/*------------------------------------------------------------*/
/* c = a*b                                                    */
/* chooses the kind of points                                 */
/*------------------------------------------------------------*/
void multiply_evaluate(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b)
{
    if (is_FFT_ready())
    {
        multiply_evaluate_FFT(c, a, b);
    }
    else
    {
        multiply_evaluate_geometric(c, a, b);
    }
}


/*------------------------------------------------------------*/
/* transpose of b mapsto c = a*b. output is                   */
/*    trunc( rev(a, dA)*c div x^dA, dB+1 )                    */
/* geometric points                                           */
/*------------------------------------------------------------*/
void t_multiply_evaluate_geometric(Mat<zz_pX> & b, const Mat<zz_pX> & a, const Mat<zz_pX> & c, long dA, long dB)
{
    long dC = dA + dB;
    if (deg(a) > dA || deg(c) > dC)
        LogicError("bad degree for t_multiply");
    long sz = dC + 1;

    zz_pX_Multipoint_Geometric ev_geom = get_geometric_points(sz);

    ev_geom.prepare_degree(dA);
    ev_geom.prepare_degree(dB);

    Vec<Mat<zz_p>> valA, valB, valC;
    ev_geom.evaluate_matrix(valA, a);
    ev_geom.t_interpolate_matrix(valC, c);

    long len = ev_geom.length();
    valB.SetLength(len);
    for (long i = 0; i < len; i++)
    {
        mul(valB[i], valA[i], valC[i]);
    }

    ev_geom.t_evaluate_matrix(b, valB);
    trunc(b, b, dB + 1);
}
