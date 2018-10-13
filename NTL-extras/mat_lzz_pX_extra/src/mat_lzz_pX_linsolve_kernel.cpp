#include <NTL/matrix.h>
#include <NTL/mat_lzz_p.h>
#include <NTL/lzz_pX.h>
#include <cmath>
#include <algorithm> // for manipulating std::vector (min, max, ..)
#include <numeric> // for std::iota

#include "lzz_p_extra.h"
#include "mat_lzz_pX_extra.h"
#include "mat_lzz_pX_partial_linearization.h"
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
    long m = pmat.NumRows();
    long n = pmat.NumCols();

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
    DegVec shift(m+1);
    // --> take strictly larger than row degree of augmented pmat
    // (requirement of kernel algo, TODO although it should not require strict)
    row_degree(shift, augmented_pmat);
    std::transform(shift.begin(), shift.end(), shift.begin(), [&](long d){return d+1;});
    // --> and add sufficiently large value to last entry to ensure that we
    // retrieve a solution of small degree
    //shift[m] += deg(b);

    // compute kernel
    Mat<zz_pX> kerbas;
    kernel_basis_zls_via_approximation(kerbas, augmented_pmat, shift);

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
