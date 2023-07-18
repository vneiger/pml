#include <NTL/lzz_pX.h>
#include <NTL/vector.h>

#include "lzz_p_extra.h"

NTL_CLIENT

bool test_order(const long order, const zz_p & a)
{
    if (a==1 && order==0)
        return true;

    long ainv;
    long status = InvModStatus(ainv, a._zz_p__rep, zz_p::modulus());
    if (status == 1)
        return (order == -1) ? true : false;

    zz_p pow;
    power(pow, a, order);
    if (pow != 1)
        return false;

    pow = a;
    for (long k = 1; k < order; ++k)
    {
        if (pow==1)
            return false;
        pow = pow*a;
    }

    return true;
}

/*------------------------------------------------------------*/
/* computes a few orders                                      */
/*------------------------------------------------------------*/
void check()
{
    long p = 65537;
    zz_p::init(p);

    {
        std::cout << "p = " << zz_p::modulus() << std::endl;
        zz_p a;

        a = to_zz_p(-1);
        std::cout << "order of " << a << ", " << test_order(order(a), a) << std::endl;

        a = random_zz_p();
        std::cout << "order of " << a << ", " << test_order(order(a), a) << std::endl;

        a = 50;
        std::cout << "order of " << a << ", " << test_order(order(a), a) << std::endl;
    }

    p = 65537*2;
    zz_p::init(p);

    {
        std::cout << "p = " << zz_p::modulus() << std::endl;
        zz_p a;

        a = to_zz_p(-1);
        std::cout << "order of " << a << ", " << test_order(order(a), a) << std::endl;

        a = random_zz_p();
        std::cout << "order of " << a << ", " << test_order(order(a), a) << std::endl;

        a = 50;
        std::cout << "order of " << a << ", " << test_order(order(a), a) << std::endl;
    }
}  

/*------------------------------------------------------------*/
/* main just calls check()                                    */
/*------------------------------------------------------------*/
int main(int argc, char** argv)
{
    check();
    return 0;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
