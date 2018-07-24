#include <NTL/lzz_pX.h>
#include <NTL/vector.h>

#include "lzz_pX_CRT.h"

NTL_CLIENT


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* multipoint evaluation over zz_p                            */
/* FFT points                                                 */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* constructor for geometric progressions                     */
/* n has to be a power of 2, p must be FFT prime              */
/*------------------------------------------------------------*/
zz_pX_Multipoint_FFT::zz_pX_Multipoint_FFT(long n)
{
    k = NextPowerOfTwo(n);
    if (n != (1L << k))
    {
	cerr << "Wrong length for zz_pX_Multipoint_FFT\n";
	exit(-1);
    }
    
    this->n = n;
    if (zz_pInfo->p_info == NULL)
    {
	cerr << "Attempt to init a zz_pX_Multipoint_FFT without zz_p::FFTInit\n";
	exit(-1);
    }
    
    zz_pX X;
    X = 0;
    wk = fftRep(INIT_SIZE, k);
    TofftRep(wk, X, k); // makes sure that wk has length n (needed for NTL >= 11)
}

/*------------------------------------------------------------*/
/* does a forward FFT                                         */
/*------------------------------------------------------------*/
void zz_pX_Multipoint_FFT::evaluate(Vec<zz_p>& val, const zz_pX& f) const 
{
    fftRep frep(INIT_SIZE, k);
    TofftRep(frep, f, k);
    long *frept = &frep.tbl[0][0];

    val.SetLength(n);
    for (long i = 0; i < n; i++)
    {
	val[i] = frept[i];
    }
}

/*------------------------------------------------------------*/
/* does an inverse FFT                                        */
/*------------------------------------------------------------*/
void zz_pX_Multipoint_FFT::interpolate(zz_pX& f, const Vec<zz_p>& val) {
    long *frept = &wk.tbl[0][0];

    for (long i = 0; i < n; i++)
    {
	frept[i] = rep(val[i]);
    }

    FromfftRep(f, wk, 0, n-1);
}


/*------------------------------------------------------------*/
/* returns a zz_pX_Multipoint_FFT with at least n points      */
/*------------------------------------------------------------*/
zz_pX_Multipoint_FFT get_FFT_points(long n)
{
    long k = NextPowerOfTwo(n);
    return zz_pX_Multipoint_FFT(1L << k);
}
