#include "util.h"
#include "mat_lzz_pX_multiply.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* transpose of b mapsto c = a*b. output is                   */
/*    trunc( rev(a, dA)*c div x^dA, dB+1 )                    */
/* output may alias input; b does not have to be zero matrix  */
/*------------------------------------------------------------*/
static void t_multiply_evaluate_do_it(Mat<zz_pX> & b, const Mat<zz_pX> & a, const Mat<zz_pX> & c, long dA, long dB, zz_pX_Multipoint& ev)
{
    long s = a.NumRows();
    long t = a.NumCols();
    long u = c.NumCols();

    long n = ev.length();

    Vec<zz_p> mat_valA, mat_valC;
    Vec<Vec<zz_p>> mat_valB;

    mat_valA.SetLength(n * s * t);
    mat_valC.SetLength(n * t * u);

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
            ev.t_interpolate(tmp, c[i][k]);
            for (long r = 0, rtu = 0; r < n; r++, rtu += tu)
                mat_valC[rtu + i*u + k] = tmp[r];
        }
    }

    Mat<zz_p> va, vb, vc;
    va.SetDims(s, t);
    vc.SetDims(t, u);

    mat_valB.SetLength(s * u);
    for (long i = 0; i < s * u; i++)
        mat_valB[i].SetLength(n);

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
                vc[i][k] = mat_valC[jtu + i*u + k];
            }
        }

        vb = va * vc;

        for (long i = 0; i < s; i++)
        {
            for (long k = 0; k < u; k++)
            {
                mat_valB[i*u + k][j] = vb[i][k];
            }
        }
    }

    b.SetDims(s, u);
    for (long i = 0; i < s; i++)
    {
        for (long k = 0; k < u; k++)
        {
            ev.t_evaluate(b[i][k], mat_valB[i*u + k], dB + 1);
        }
    }
}

/*------------------------------------------------------------*/
/* transpose of b mapsto c = a*b. output is                   */
/*    trunc( rev(a, dA)*c div x^dA, dB+1 )                    */
/* geometric points                                           */
/* output may alias input; b does not have to be zero matrix  */
/*------------------------------------------------------------*/
void t_multiply_evaluate_geometric(Mat<zz_pX> & b, const Mat<zz_pX> & a, const Mat<zz_pX> & c, long dA, long dB)
{
    long dC = dA + dB;
    if (deg(a) > dA || deg(c) > dC)
        LogicError("bad degree for t_multiply");
    long sz = dC + 1;

    zz_pX_Multipoint *ev;
    zz_pX_Multipoint_Geometric ev_geometric = get_geometric_points(sz);

    ev_geometric.prepare_degree(dA);
    ev_geometric.prepare_degree(dB);

    ev = &ev_geometric;
    t_multiply_evaluate_do_it(b, a, c, dA, dB, *ev);
}


// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
