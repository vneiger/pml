#include <NTL/matrix.h>
#include <NTL/mat_lzz_p.h>
#include <NTL/lzz_pX.h>

#include "util.h"
#include "lzz_p_extra.h"
#include "mat_lzz_pX_extra.h"
#include "lzz_pX_CRT.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* c = a*b                                                    */
/*------------------------------------------------------------*/
void multiply_evaluate_do_it(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b, zz_pX_Multipoint& ev)
{
    long s = a.NumRows();
    long t = a.NumCols();
    long u = b.NumCols();
    c.SetDims(s, u);

    long n = ev.length();

    Vec<zz_p> mat_valA, mat_valB;
    Vec<Vec<zz_p>> mat_valC;

    mat_valA.SetLength(n * s * t);
    mat_valB.SetLength(n * t * u);

    Vec<zz_p> tmp;
    long st = s*t;
    for (long i = 0; i < s; i++)
    {
        for (long k = 0; k < t; k++)
        {
            ev.evaluate(tmp, a[i][k]);
            for (long r = 0, rst = 0; r < n; r++, rst += st)
                mat_valA[rst + i*t + k] = tmp[r];
        }
    }

    long tu = t*u;
    for (long i = 0; i < t; i++)
    {
        for (long k = 0; k < u; k++)
        {
            ev.evaluate(tmp, b[i][k]);
            for (long r = 0, rtu = 0; r < n; r++, rtu += tu)
                mat_valB[rtu + i*u + k] = tmp[r];
        }
    }

    Mat<zz_p> va, vb, vc;
    va.SetDims(s, t);
    vb.SetDims(t, u);

    mat_valC.SetLength(s * u);
    for (long i = 0; i < s * u; i++)
        mat_valC[i].SetLength(n);

    for (long j = 0, jst = 0, jtu = 0; j < n; j++, jst += st, jtu += tu)
    {
        for (long i = 0; i < s; i++)
        {
            for (long k = 0; k < t; k++)
            {
                va[i][k] = mat_valA[jst + i*t + k];
            }
        }
        for (long i = 0; i < t; i++)
        {
            for (long k = 0; k < u; k++)
            {
                vb[i][k] = mat_valB[jtu + i*u + k];
            }
        }

        vc = va * vb;

        for (long i = 0; i < s; i++)
        {
            for (long k = 0; k < u; k++)
            {
                mat_valC[i*u + k][j] = vc[i][k];
            }
        }
    }

    for (long i = 0; i < s; i++)
    {
        for (long k = 0; k < u; k++)
        {
            ev.interpolate(c[i][k], mat_valC[i*u + k]);
        }
    }
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

    zz_pX_Multipoint *ev;
    zz_pX_Multipoint_FFT ev_FFT = get_FFT_points(sz);

    ev = &ev_FFT;
    multiply_evaluate_do_it(c, a, b, *ev);
}

/*------------------------------------------------------------*/
/* c = a*b                                                    */
/* geometric points, uses FFTs directly                       */
/*------------------------------------------------------------*/
void multiply_evaluate_geometric_using_FFT(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b)
{
    long dA = deg(a);
    long dB = deg(b);
    long dC = dA+dB;
    long sz = dC+1;

    zz_pX_Multipoint *ev;
    zz_pX_Multipoint_Geometric ev_geometric = get_geometric_points(sz);

    ev_geometric.prepare_degree(dA);
    ev_geometric.prepare_degree(dB);

    ev_geometric.set_FFT_evaluate();
    ev_geometric.set_FFT_interpolate();

    ev = &ev_geometric;
    multiply_evaluate_do_it(c, a, b, *ev);
}

/*------------------------------------------------------------*/
/* c = a*b                                                    */
/* geometric points, uses middle product                      */
/*------------------------------------------------------------*/
void multiply_evaluate_geometric_no_FFT(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b)
{
    long dA = deg(a);
    long dB = deg(b);
    long dC = dA+dB;
    long sz = dC+1;

    zz_pX_Multipoint *ev;
    zz_pX_Multipoint_Geometric ev_geometric = get_geometric_points(sz);

    ev_geometric.prepare_degree(dA);
    ev_geometric.prepare_degree(dB);

    ev_geometric.unset_FFT_evaluate();
    ev_geometric.unset_FFT_interpolate();

    ev = &ev_geometric;
    multiply_evaluate_do_it(c, a, b, *ev);
}

/*------------------------------------------------------------*/
/* c = a*b                                                    */
/* geometric points, thresholded                              */
/*------------------------------------------------------------*/
void multiply_evaluate_geometric(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b)
{
    long dA = deg(a);
    long dB = deg(b);
    long dC = dA+dB;
    long sz = dC+1;

    zz_pX_Multipoint *ev;
    zz_pX_Multipoint_Geometric ev_geometric = get_geometric_points(sz);

    ev_geometric.prepare_degree(dA);
    ev_geometric.prepare_degree(dB);

    ev = &ev_geometric;
    multiply_evaluate_do_it(c, a, b, *ev);
}

/*------------------------------------------------------------*/
/* c = a*b                                                    */
/* chooses the kind of points                                 */
/* assumes that the field is large enough                     */
/*------------------------------------------------------------*/
void multiply_evaluate(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b)
{
    if (is_FFT_ready(NextPowerOfTwo(deg(a) + deg(b) + 1)))
    {
        multiply_evaluate_FFT(c, a, b);
    }
    else
    {
        multiply_evaluate_geometric(c, a, b);
    }
}


// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
