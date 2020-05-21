#include "lzz_pX_seq.h"
#include <NTL/lzz_pX.h>
#include <iostream>
#include <sstream>
#include <NTL/vector.h>

NTL_CLIENT
using namespace std;

int main(){
	long p = 13;
	zz_p::init(p);

	long d;
	cout << "Enter d: ";
	cin >> d;

	long ell;
	cout << "Enter ell: ";
	cin >> ell;

	Vec<zz_pX> seq;
	seq.SetLength(ell);
	for (long i = 0; i < ell; i++){
		cout << "enter seq[" << i << "]:"<< endl;
		cin >> seq[i];
		trunc(seq[i],seq[i],d);
	}
	cout << "SEQ: " << endl;
	for (long i = 0; i < ell; i++)
		cout << seq[i] << endl;

	cout << "STARTING KURAKIN" << endl;
	Vec<zz_pXY> gens;
	kurakin(d,seq,gens);

	cout << "GENS:" << endl;
	for (long i = 0; i < gens.length(); i++)
		cout << gens[i] << endl;

	Vec<zz_pXY> gens2;
	modified_kurakin(d,seq,gens2);

	cout << "GENS2:"<< endl;
	for (long i = 0; i < gens2.length(); i++)
		cout << gens2[i] << endl;
	
	cout << "GENS2 filled in:" << endl;
	fill_in(gens2);
	for (long i  = 0; i < gens2.length(); i++)
		cout << gens2[i] << endl;

	cout << "checking that gens cancels the seq" << endl;
	if (check_cancel(seq,gens,d)) cout << "OKAY!" << endl;
	else cout << "BAD!" << endl;

	cout << "checking that gens2 cancels the seq" << endl;
	if (check_cancel(seq,gens2,d)) cout << "OKAY!" << endl;
	else cout << "BAD!" << endl;

}
