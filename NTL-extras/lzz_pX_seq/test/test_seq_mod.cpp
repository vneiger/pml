#include "lzz_pX_seq.h"
#include <NTL/lzz_pX.h>
#include <iostream>
#include <sstream>
#include <NTL/vector.h>

NTL_CLIENT
using namespace std;

int main(){
	long p = 9001;
	zz_p::init(p);

	bool verbose = true;

	long d;
	if (verbose)
		cout << "Enter d: ";
	cin >> d;

	long ell;
	if (verbose)
		cout << "Enter ell: ";
	cin >> ell;

	long tau;
	if (verbose)
		cout << "Enter tau: ";
	cin >> tau;


	cout << "d: " << d << endl;
	cout << "ell: " << ell << endl;
	cout << "tau: " << tau << endl;

	Vec<Vec<zz_pX>> S;
	for (long i = 0; i < ell; i++){
		if (verbose) cout << "At entry: " << i << endl;
		Vec<zz_pX> cur;
		for (long j = 0; j < tau; j++){
			zz_pX curpol;
			cout << "subentry: " << j << endl;
			cin >> curpol;
			cur.append(curpol);
		}
		S.append(cur);
	}
	Vec<zz_pXY> gens;
	kurakin(d,S,gens);
	cout << "gens:" << endl;
	for (long i = 0; i < gens.length(); i++)
		cout << gens[i] << endl;

}
