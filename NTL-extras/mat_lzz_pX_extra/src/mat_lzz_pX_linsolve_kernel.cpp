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
                         Mat<zz_pX> pmat,
                         const Vec<zz_pX> &b
                        )
{
    if (b.length() != pmat.NumRows())
        throw std::logic_error("length of b != pmat.NumRows()");

    long m = pmat.NumRows();
    long n = pmat.NumCols();
    if (m != n)
        throw std::logic_error("pmat must be square");

    long degb = deg(b[0]);
    for (long i = 1; i < n; i++)
        if (degb < deg(b[i])) degb = deg(b[i]);

    Shift shift;
    shift.resize(m+1);

    DegVec rdeg;
    rdeg.resize(m);
    row_degree(rdeg,pmat);
    long max_deg = rdeg[0];
    for (unsigned long i = 0; i < rdeg.size(); i++)
    {
        shift[i] = rdeg[i]+1;
        if (max_deg < rdeg[i]) max_deg = rdeg[i];
    }
    shift[m] = (m+1)*(max_deg+1);

    pmat.SetDims(m+1, n);
    pmat[m] = b;

    // compute kernel
    Mat<zz_pX> kerbas;
    auto deg = kernel_basis(kerbas, pmat, shift);
    //auto deg = kernel_basis_zls_via_approximation(kerbas, pmat, shift);

    Mat<zz_pX> res;
    multiply(res,kerbas,pmat);
    cout << "kernel prod: " << degree_matrix(res) << endl;

    a.SetLength(n);
    for (long i = 0; i < n; i++)
        a[i] = -kerbas[0][i];
    d = kerbas[0][n];
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
