#include <NTL/matrix.h>
#include <NTL/mat_lzz_p.h>
#include <NTL/lzz_pX.h>

#include "mat_lzz_pX_extra.h"

void diagonal_of_hermite(Vec<zz_pX> & diag, const Mat<zz_pX> & pmat)
{
    diag.SetLength(pmat.NumRows());
    // TODO
    throw std::invalid_argument("==diagonal_of_hermite== Not implemented yet");
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
