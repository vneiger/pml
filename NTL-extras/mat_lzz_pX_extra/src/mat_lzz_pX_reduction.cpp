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

    // compute the 2d high-degree coefficients of the truncated inverse of
    // pmat, precisely, those of degree (dim-1)*d+1, ..., (dim-1)*d + 2*d
    // TODO try with high-order lifting, compare speed
    Mat<zz_pX> series;
    inv_trunc(series, pmat, (dim+1)*d + 1);
    RightShift(series, series, (dim-1)*d + 1);

    // Reconstruction via approximants: add identity below series
    series.SetDims(2*dim, dim);
    for (long i = 0; i < dim; ++i)
        set(series[dim+i][i]);

    // Reconstruction
    Mat<zz_pX> appbas;
    VecLong shift(2*dim);
    std::cout << "appbas: " << 2*dim << "," << dim << "," << 2*d << std::endl;
    pmbasis(appbas, series, 2*d, shift);

    // Extract the reduced form
    reduced.SetDims(dim, dim);
    for (long i = 0; i < dim; ++i)
        for (long j = 0; j < dim; ++j)
            reduced[i][j] = appbas[i][j];

    return 1;
}


/** Helper for high order lifting. As input, we have a polynomial matrix `pmat`
 * of degree `d`, and its truncated inverse `inv` at order `d`, that is, `inv =
 * pmat^{-1} mod x^d` (see #TruncatedInverse).
 *
 * For a given nonnegative integer `i`, we define `Si` to be the slice of the
 * expansion of `pmat^{-1}` of its terms of degree between `i-(2d-1)` and `i-1`
 * (both included).  Given `src = Si`, this function computes `next =
 * S_{2i-d}`. Note that deg(Si) < 2d-1.
 *
 * The OUT parameter `next` can alias the IN parameter `src`.
 */
//void high_order_lift_inverse_odd(
//                                 Mat<zz_pX> & next,
//                                 const Mat<zz_pX> & src, 
//                                 const std::unique_ptr<mat_lzz_pX_lmultiplier> & pmat, 
//                                 const std::unique_ptr<mat_lzz_pX_lmultiplier> & inv,
//                                 long d
//                                );



// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
