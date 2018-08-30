#include <NTL/matrix.h>
#include <NTL/mat_lzz_p.h>
#include <NTL/lzz_pX.h>

#include "mat_lzz_pX_extra.h"

void determinant_via_diagonal_of_hermite(zz_pX & det, const Mat<zz_pX> & pmat)
{
    Vec<zz_pX> diag;
    diagonal_of_hermite(diag, pmat);
    set(det); // det = 1
    for (long i = 0; i < diag.length(); ++i)
        det *= diag[i];
}

bool determinant_generic_knowing_degree(zz_pX & det, const Mat<zz_pX> & pmat, long degree)
{
    long dim = pmat.NumRows();
    if (dim==1)
    {
        det = pmat[0][0];
        return (degree==deg(det));
    }
    else
    {
        long dimm = (dim>>1);
        // TODO write and use "generic" kernel
        Mat<zz_pX> ker;
        Mat<zz_pX> pmat_l;
        Mat<zz_pX> pmat_r;
        ker.SetDims(dim,dim);
        pmat_l.SetDims(dim,dim-dimm);
        pmat_r.SetDims(dim,dimm);

        for (long i = 0; i < dim; ++i)
        {
            for (long j = 0; j < dim-dimm; ++j)
                pmat_l[i][j] = pmat[i][j];
            for (long j = dimm; j < dim; ++j)
                pmat_r[i][j-dimm] = pmat[i][j];
        }

        // compute the kernel via approximant basis at high order
        // TODO is computing the degree at each recursion level necessary?
        // (goes over the whole matrix each time... info could be transmitted
        // through the successive calls)
        pmbasis(ker, pmat_r, 2*deg(pmat_r), Shift(dim,0));

        // only keep the first rows (minimal left kernel basis of pmat_r),
        // and compute the product
        ker.SetDims(dim-dimm,dim);
        Mat<zz_pX> pmatt;
        multiply(pmatt, ker, pmat_l);

        return determinant_generic_knowing_degree(det,pmatt,degree);
    }
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
