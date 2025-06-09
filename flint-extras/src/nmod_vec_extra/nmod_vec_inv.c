#include <flint/nmod.h>
#include <flint/nmod_vec.h>
#include <flint/ulong_extras.h>

#include "nmod_vec_extra.h"

void _nmod_vec_inv(nn_ptr res, nn_srcptr vec, ulong len, nmod_t mod)
{
    if (len > 0)
    {
        /* fill tmp with tmp[k] = vec[0] * ... * vec[k] */
        nn_ptr tmp = _nmod_vec_init(len);
        tmp[0] = vec[0];
        for (ulong k = 1; k < len; k++)
            tmp[k] = nmod_mul(vec[k], tmp[k-1], mod);

        /* current inverse is (vec[0] * ... * vec[len-1])^{-1} */
        ulong inv = nmod_inv(tmp[len-1], mod);

        /* at the end of iteration k, 
         * res[k] = (vec[0] * ... * vec[k-1]) * (vec[0] * ... * vec[k])^{-1}
         *        = vec[k]^{-1}
         * inv = (vec[0] * ... * vec[k-1])^{-1}
         **/
        for (long k = len-1; k > 0; k--)
        {
            res[k] = nmod_mul(inv, tmp[k-1], mod);
            inv = nmod_mul(inv, vec[k], mod);
        }
        res[0] = inv;
        _nmod_vec_clear(tmp);
    }
}

void _nmod_vec_inv2(nn_ptr res, nn_srcptr vec, ulong len, nmod_t mod)
{
    if (len > 0)
    {
        /* fill tmp with 
         *    tmp[4*k+0] = precomp quo for mul by vec[k]
         *    tmp[4*k+1] = precomp rem for mul by vec[k]
         *    tmp[4*k+2] = vec[0] * ... * vec[k]
         *    tmp[4*k+3] = precomp quo for mul by tmp[4*k+2]
         */
        nn_ptr tmp = _nmod_vec_init(4*len);
        n_mulmod_precomp_shoup_quo_rem(tmp+0, tmp+1, vec[0], mod.n);
        tmp[2] = vec[0];
        tmp[3] = tmp[0];
        for (ulong k = 1; k < len; k++)
        {
            n_mulmod_precomp_shoup_quo_rem(tmp+4*k, tmp+4*k+1, vec[k], mod.n);
            n_mulmod_and_precomp_shoup(tmp+4*k+2, tmp+4*k+3, vec[k], tmp[4*k-2], tmp[4*k], tmp[4*k+1], tmp[4*k-1], mod.n);
        }

        /* current inverse is (vec[0] * ... * vec[len-1])^{-1} */
        ulong inv = nmod_inv(tmp[4*len-2], mod);

        /* at the end of iteration k, 
         * res[k] = (vec[0] * ... * vec[k-1]) * (vec[0] * ... * vec[k])^{-1}
         *        = vec[k]^{-1}
         * inv = (vec[0] * ... * vec[k-1])^{-1}
         **/
        for (long k = len-1; k > 0; k--)
        {
            res[k] = n_mulmod_shoup(tmp[4*k-2], inv, tmp[4*k-1], mod.n);
            inv = n_mulmod_shoup(vec[k], inv, tmp[4*k], mod.n);
        }
        res[0] = inv;
        _nmod_vec_clear(tmp);
    }
}
