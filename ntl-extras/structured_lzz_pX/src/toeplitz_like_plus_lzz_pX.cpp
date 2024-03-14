#include <NTL/lzz_pX.h>
#include <NTL/mat_lzz_p.h>

#include "lzz_p_extra.h"
#include "mat_lzz_p_extra.h"
#include "mat_lzz_pX_extra.h"
#include "structured_lzz_p.h"
#include "structured_lzz_pX.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* returns Z1 A - A Z0                                        */
/*------------------------------------------------------------*/
void toeplitz_lzz_pX_phi_plus(Mat<zz_pX> & res, const Mat<zz_pX>& A)
{
    if (&res == &A)
    {
        res = toeplitz_lzz_pX_phi_plus(A);
        return;
    }

    long m = A.NumRows();
    long n = A.NumCols();
    res = Z_lzz_p(m, to_zz_p(1)) * A - A * Z_lzz_p(n, to_zz_p(0));
}
