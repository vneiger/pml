#include <NTL/ZZX.h>
#include <NTL/lzz_pX.h>
#include <NTL/vector.h>
#include <NTL/matrix.h>

#include "sage_output.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* initializes GF(p)                                          */
/*------------------------------------------------------------*/
void sage_init()
{
    cout << "k = GF(" << zz_p::modulus() << ")\n";
}

/*------------------------------------------------------------*/
/* initializes U.<x>=GF(p)[x]                                 */
/*------------------------------------------------------------*/
void sage_init_X()
{
    cout << "U.<x> = PolynomialRing(k)\n";
}

/*------------------------------------------------------------*/
/* initializes U<XX>=Z[XX]                                    */
/*------------------------------------------------------------*/
void sage_init_ZZX()
{
    cout << "U.<XX> = PolynomialRing(Integers())\n";
}


/*------------------------------------------------------------*/
/* prints a vector                                            */
/*------------------------------------------------------------*/
void sage_output(const Vec<ZZ> & v)
{
    if (v.length() == 0)
    {
        cout << "[]";
        return;
    }
    cout << "[";
    for (long i = 0; i < v.length()-1; i++)
    {
        cout << v[i] << ", ";
    }
    cout << v[v.length()-1] << "]";
}

void sage_output(const Vec<long> & v)
{
    if (v.length() == 0)
    {
        cout << "[]";
        return;
    }
    cout << "[";
    for (long i = 0; i < v.length()-1; i++)
    {
        cout << v[i] << ", ";
    }
    cout << v[v.length()-1] << "]";
}

void sage_output(const Vec<zz_p> & v)
{
    if (v.length() == 0)
    {
        cout << "[]";
        return;
    }
    cout << "[";
    for (long i = 0; i < v.length()-1; i++)
    {
        cout << v[i] << ", ";
    }
    cout << v[v.length()-1] << "]";
}

void sage_output(const Vec<unsigned long> & v)
{
    if (v.length() == 0)
    {
        cout << "[]";
        return;
    }
    cout << "[";
    for (long i = 0; i < v.length()-1; i++)
    {
        cout << v[i] << ", ";
    }
    cout << v[v.length()-1] << "]";
}


/*------------------------------------------------------------*/
/* assigns a vector to variable "name"                        */
/*------------------------------------------------------------*/
void sage_assign(const Vec<ZZ> & v, const string & name)
{
    cout << name << " = ";
    sage_output(v);
    cout  << endl;
}

void sage_assign(const Vec<long> & v, const string & name)
{
    cout << name << " = ";
    sage_output(v);
    cout  << endl;
}

void sage_assign(const Vec<unsigned long> & v, const string & name)
{
    cout << name << " = ";
    sage_output(v);
    cout  << endl;
}

void sage_assign(const Vec<zz_p> & v, const string & name)
{
    cout << name << " = ";
    sage_output(v);
    cout  << endl;
}


/*------------------------------------------------------------*/
/* prints a matrix                                            */
/*------------------------------------------------------------*/
void sage_output(const Mat<zz_p> & v)
{
    if (v.NumRows() == 0)
    {
        cout << "Matrix(GF(" << zz_p::modulus() << "), [[]])";
        return;
    }
    cout << "Matrix(GF(" << zz_p::modulus() << "), [";
    for (long i = 0; i < v.NumRows()-1; i++)
    {
        sage_output(v[i]);
        cout << ", ";
    }
    sage_output(v[v.NumRows()-1]);
    cout << "])";
}

/*------------------------------------------------------------*/
/* assigns a matrix to variable "name"                        */
/*------------------------------------------------------------*/
void sage_assign(const Mat<zz_p> & v, const string & name)
{
    cout << name << " = ";
    sage_output(v);
    cout << endl;
}



/*------------------------------------------------------------*/
/* prints a poly / vector / .. with indeterminate var         */
/*------------------------------------------------------------*/
void sage_output(const ZZX & v, const string & var)
{
    cout << "(" << var << ".parent()(0)";
    for (long i = 0; i <= deg(v); i++)
    {
        cout << "+(" << coeff(v, i) << ")*" << var << "^" << i;
    }
    cout << ")";
}

void sage_output(const zz_pX & v, const string & var)
{
    cout << "(" << var << ".parent()(0)";
    for (long i = 0; i <= deg(v); i++)
    {
        cout << "+(" << coeff(v, i) << ")*" << var << "^" << i;
    }
    cout << ")";
}

void sage_output(const Vec<zz_pX> & v, const string & var)
{
    cout << "[";
    for (long i = 0; i < v.length(); i++)
    {
        sage_output(v[i]);
        if (i < v.length()-1)
            cout << ", ";
    }
    cout << "]";
}

void sage_output(const Mat<zz_pX>& a, const string & var)
{
    cout << "Matrix([";
    for (long i = 0; i < a.NumRows(); i++)
    {
        cout << "[";
        for (long j = 0; j < a.NumCols(); j++)
        {
            sage_output(a[i][j], var);
            if (j < a.NumCols()-1)
            {
                cout << ", ";
            }
        }
        cout << "]";
        if (i < a.NumRows()-1)
        {
            cout << ", ";
        }
    }
    cout << "])";
}

/*------------------------------------------------------------*/
/* prints a poly / .. in the right variable                   */
/*------------------------------------------------------------*/
void sage_output(const ZZX & v)
{
    cout << "UZ(";
    sage_output(v, "XX");
    cout << ")";
}

void sage_output(const zz_pX & v)
{
    cout << "U(";
    sage_output(v.rep);
    cout << ")";
}

void sage_output(const Vec<zz_pX> & v)
{
    sage_output(v, "x");
}

void sage_output(const Mat<zz_pX> & v)
{
    sage_output(v, "x");
}

/*------------------------------------------------------------*/
/* assign a poly with indeterminate var to variable "name"    */
/*------------------------------------------------------------*/
void sage_assign(const ZZX & v, const string & var, const string & name)
{
    cout << name << " = ";
    sage_output(v, var);
    cout << endl;
} 

void sage_assign(const zz_pX & v, const string & var, const string & name)
{
    cout << name << " = ";
    sage_output(v, var);
    cout << endl;
} 

void sage_assign(const Vec<zz_pX> & v, const string & var, const string & name)
{
    cout << name << "=";
    sage_output(v, var);
    cout << endl;
}

void sage_assign(const Mat<zz_pX>& a, const string & var, const string & name)
{
    cout << name << "=";
    sage_output(a, var);
    cout << endl;
}

/*------------------------------------------------------------*/
/* assign a poly with indeterminate XX to variable "name"     */
/*------------------------------------------------------------*/
void sage_assign(const ZZX & v, const string & name)
{
    sage_assign(v, "XX", name);
}

void sage_assign(const zz_pX & v, const string & name)
{
    sage_assign(v, "x", name);
}

void sage_assign(const Vec<zz_pX> & v, const string & name)
{
    sage_assign(v, "x", name);
}

void sage_assign(const Mat<zz_pX>& v, const string & name)
{
    sage_assign(v, "x", name);
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
