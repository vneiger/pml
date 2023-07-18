#include "mat_lzz_pX_extra.h"

NTL_CLIENT

long reduced_form_gjv(
                      Mat<zz_pX> & reduced,
                      const Mat<zz_pX> & pmat
                     )
{
    if (determinant(coeff(pmat,0)) == 0)
        return 0;

    // dimensions
    const long dim = pmat.NumCols();
    const long d = deg(pmat);
    std::cout << dim << "\t" << d << "\t";

    // inv_trunc version
    // compute the 2d high-degree coefficients of the truncated inverse of
    // pmat, precisely, those of degree (dim-1)*d+1, ..., (dim-1)*d + 2*d
    double t=GetWallTime();
    Mat<zz_pX> series;
    inv_trunc(series, pmat, (dim+1)*d + 1);
    RightShift(series, series, (dim-1)*d + 1);
    std::cout << GetWallTime()-t << "\t";

    // HOL version
    t=GetWallTime();
    Mat<zz_pX> buf;
    inv_trunc(buf, pmat, 2*d+2);
    RightShift(series,buf,1);
    trunc(buf, buf, d+1);

    // prepare left multipliers mod X^d
    std::unique_ptr<mat_lzz_pX_lmultiplier> minv = get_lmultiplier(buf,d);
    std::unique_ptr<mat_lzz_pX_lmultiplier> mpmat = get_lmultiplier(pmat,d);

    // get 2d+1 terms of the expansion starting from degree at least (dim-1)*d+1
    // for the moment we have the terms i-(2d+1)...i-1 for i = 2d+2
    // we want to lift until reaching i-(2d+1) >= (dim-1)*d+1, that is, i >= (dim+1)*d + 2
    for (long i=2*d+2; i < (dim+1)*d+2; i = 2*i-d-1)
        high_order_lift_inverse_odd(series, series, mpmat, minv, d+1);
    std::cout << GetWallTime()-t << "\t";

    // Reconstruction via approximants: add identity below series
    series.SetDims(2*dim, dim);
    for (long i = 0; i < dim; ++i)
        set(series[dim+i][i]);

    // Reconstruction
    Mat<zz_pX> appbas;
    VecLong shift(2*dim);
    t=GetWallTime();
    pmbasis(appbas, series, 2*d, shift);
    std::cout << GetWallTime()-t << "\t";

    // Extract the reduced form
    reduced.SetDims(dim, dim);
    for (long i = 0; i < dim; ++i)
        for (long j = 0; j < dim; ++j)
            reduced[i][j] = appbas[i][j];

    return 1;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
