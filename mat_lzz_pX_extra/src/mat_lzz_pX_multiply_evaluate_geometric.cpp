#include "mat_lzz_pX_multiply.h"

NTL_CLIENT


/*------------------------------------------------------------*/
/* c = a*b                                                    */
/* output may alias input; c does not have to be zero matrix  */
/* uses Mat<zz_p> matrix multiplication                       */
/*------------------------------------------------------------*/
void multiply_evaluate_geometric(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b)
{
    // dimensions
    const long s = a.NumRows();
    const long t = a.NumCols();
    const long u = b.NumCols();
    // degrees
    const long dA = deg(a);
    const long dB = deg(b);
    const long sz = dA+dB+1;

    // if a or b is zero, return
    if (dA < 0 || dB < 0)
    {
        c.SetDims(s, t);
        clear(c);
        return;
    }

    // prepare multipoint evaluation at geometric points
    zz_pX_Multipoint_Geometric ev = get_geometric_points(sz);
    ev.prepare_degree(dA);
    ev.prepare_degree(dB);
    const long n = ev.length();

    // vector to store evaluations
    Vec<zz_p> evs(INIT_SIZE, n);

    // stores evaluations in a single vector for a, same for b
    const long st = s*t;
    Vec<zz_p> mat_valA(INIT_SIZE, n*st);
    const long tu = t*u;
    Vec<zz_p> mat_valB(INIT_SIZE, n*tu);

    // mat_valA[r*s*t + i*t + k] is a[i][k] evaluated at the r-th point
    for (long i = 0; i < s; ++i)
        for (long k = 0; k < t; ++k)
        {
            ev.evaluate(evs, a[i][k]);
            for (long r = 0, rst = 0; r < n; ++r, rst += st)
                mat_valA[rst + i*t + k] = evs[r];
        }

    // mat_valB[r*s*t + i*t + k] is b[i][k] evaluated at the r-th point
    for (long i = 0; i < t; ++i)
        for (long k = 0; k < u; k++)
        {
            ev.evaluate(evs, b[i][k]);
            for (long r = 0, rtu = 0; r < n; ++r, rtu += tu)
                mat_valB[rtu + i*u + k] = evs[r];
        }
    
    // will store a evaluated at the r-th point, same for b
    Mat<zz_p> va(INIT_SIZE, s, t);
    Mat<zz_p> vb(INIT_SIZE, t, u);
    // will store the product evaluated at the r-th point, i.e. va*vb
    Mat<zz_p> vc(INIT_SIZE, s, u);

    Vec<Vec<zz_p>> mat_valC(INIT_SIZE, s*u);
    for (long i = 0; i < s * u; ++i)
        mat_valC[i].SetLength(n);
    
    // compute pairwise products
    for (long j = 0, jst = 0, jtu = 0; j < n; ++j, jst += st, jtu += tu)
    {
        for (long i = 0; i < s; ++i)
            for (long k = 0; k < t; ++k)
                va[i][k].LoopHole() = mat_valA[jst + i*t + k]._zz_p__rep;
        for (long i = 0; i < t; ++i)
            for (long k = 0; k < u; ++k)
                vb[i][k].LoopHole() = mat_valB[jtu + i*u + k]._zz_p__rep;
        
        mul(vc, va, vb);

        for (long i = 0; i < s; ++i)
            for (long k = 0; k < u; ++k)
                mat_valC[i*u + k][j].LoopHole() = vc[i][k]._zz_p__rep;
    }

    // interpolate the evaluations stored in mat_valC back into c
    c.SetDims(s, u);

    for (long i = 0; i < s; ++i)
        for (long k = 0; k < u; ++k)
            ev.interpolate(c[i][k], mat_valC[i*u + k]);
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
