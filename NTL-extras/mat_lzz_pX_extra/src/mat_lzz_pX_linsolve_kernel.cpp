#include <NTL/matrix.h>
#include <NTL/mat_lzz_p.h>
#include <NTL/lzz_pX.h>
#include <cmath>
#include <algorithm> // for manipulating std::vector (min, max, ..)
#include <numeric> // for std::iota

#include "lzz_p_extra.h"
#include "mat_lzz_pX_extra.h"
#include "lzz_pX_CRT.h"

// TODO make decision about input type (mat vs vec)

// solve aM = b via kernel basis
// return a and denominator d
// assumes M is invertible
void linsolve_via_kernel(
                         Vec<zz_pX> &a,
                         zz_pX &d,
                         const Mat<zz_pX> & pmat,
                         const Vec<zz_pX> & b
                        )
{
    const long m = pmat.NumRows();
    const long n = pmat.NumCols();

    // TODO this kind of check should probably be in higher level function
    if (m != n)
        throw std::logic_error("==linsolve_via_kernel== pmat must be square");
    if (b.length() != m)
        throw std::logic_error("==linsolve_via_kernel== length of b != pmat.NumRows()");

    // compute augmented matrix (block with input matrix and system 'b')
    Mat<zz_pX> augmented_pmat;
    augmented_pmat.SetDims(m+1, m);
    for (long i = 0; i < m; ++i)
        augmented_pmat[i] = pmat[i];
    augmented_pmat[m] = b;

    // compute shift to make sure kernel corresponds to solution
    // --> row degree of augmented matrix, with large value added to last entry
    VecLong shift(m+1);
    row_degree(shift, augmented_pmat);
    shift[m] += deg(b);

    // compute kernel
    Mat<zz_pX> kerbas;
    VecLong pivind, pivdeg;
    kernel_basis_zls_via_approximation(kerbas, augmented_pmat, shift, pivind, pivdeg);

    // deduce solution
    a.SetLength(m);
    for (long i = 0; i < m; ++i)
        a[i] = -kerbas[0][i];
    d = kerbas[0][m];
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
