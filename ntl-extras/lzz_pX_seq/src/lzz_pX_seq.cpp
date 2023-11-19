#include "lzz_pX_seq.h"
#include "mat_lzz_pX_approximant.h"
#include <map>
#include <vector>
#include <iostream>
#include <cmath>
#include <queue>
#include <utility>
#include "mat_lzz_pX_kernel.h"
using namespace std;

const bool verbose = false;
//const bool verbose = true;


Vec<zz_pX> mul(const Vec<zz_pX> &S, const zz_pX &a){
  Vec<zz_pX> res;
  res.SetLength(S.length());

  for (long i = 0; i < S.length(); i++)
    res[i] = S[i]*a;

  return res;
}

Vec<zz_pX> mul(const Vec<zz_pX> &S, const zz_pXY &f){
  Vec<zz_pX> res;
  long dy = f.degY();
  res.SetLength(min(0, S.length() - dy));
  for (long i = 0; i < res.length(); i++){
    zz_pX val{0};
    for (long t = 0; t <= dy; t++)
      val += f.rep[t] * S[i+t];
    res[i] = val;
  }

  return res;
}

Vec<Vec<zz_pX>> mul(const Vec<Vec<zz_pX>> &S, const zz_pX &a){
  Vec<Vec<zz_pX>> res;
  res.SetLength(S.length());

  for (long i = 0; i < S.length(); i++)
    res[i] = mul(S[i], a);

  return res;
}

Vec<Vec<zz_pX>> mul(const Vec<Vec<zz_pX>> &S, const zz_pXY &f){
  Vec<Vec<zz_pX>> res;
  long dy = f.degY();
  res.SetLength(min(0, S.length() - dy));
  for (long i = 0; i < res.length(); i++){
    Vec<zz_pX> val;
    val.SetLength(S[i].length());
    for (long t = 0; t <= dy; t++)
      val = add(val, mul(S[i+t], f.rep[t]));
    res[i] = val;
  }

  return res;
}

Vec<zz_pX> mulTrunc(const Vec<zz_pX> &S, const zz_pX &a, const long d){
  Vec<zz_pX> res;
  res.SetLength(S.length());

  for (long i = 0; i < S.length(); i++)
    res[i] = MulTrunc(S[i],a,d);

  return res;

}

Vec<zz_pX> mulTrunc(const Vec<zz_pX> &S, const zz_pXY &f, const long d){
  Vec<zz_pX> res;
  long dy = f.degY();
  res.SetLength(min(0, S.length() - dy));
  for (long i = 0; i < res.length(); i++){
    zz_pX val{0};
    for (long t = 0; t <= dy; t++)
      MulTrunc(val, f.rep[t], S[i+t], d);
    res[i] = val;
  }

  return res;
}

Vec<Vec<zz_pX>> mulTrunc(const Vec<Vec<zz_pX>> &S, const zz_pX &a,
    const long d){
  Vec<Vec<zz_pX>> res;
  res.SetLength(S.length());

  for (long i = 0; i < S.length(); i++)
    res[i] = mulTrunc(S[i], a, d);

  return res;
}

Vec<Vec<zz_pX>> mulTrunc(const Vec<Vec<zz_pX>> &S, const zz_pXY &f,
    const long d){
  Vec<Vec<zz_pX>> res;
  long dy = f.degY();
  if (dy < 0) return Vec<Vec<zz_pX>>();

  res.SetLength(max(0, S.length() - dy));
  for (long i = 0; i < res.length(); i++){
    Vec<zz_pX> val;
    val.SetLength(S[i].length());
    for (long t = 0; t <= dy; t++)
      val = add(val, mulTrunc(S[i+t], f.rep[t],d));
    res[i] = val;
  }

  return res;
}

Vec<zz_pX> add(const Vec<zz_pX> &S, const Vec<zz_pX> &T){
  Vec<zz_pX> res;
  long n = min(S.length(), T.length());
  res.SetLength(n);
  for (long i = 0; i < n; i++)
    res[i] = S[i] + T[i];
  return res;
}

Vec<Vec<zz_pX>> add(const Vec<Vec<zz_pX>> &S, const Vec<Vec<zz_pX>> &T){
  Vec<Vec<zz_pX>> res;
  long n = min(S.length(), T.length());
  res.SetLength(n);
  for (long i = 0; i < n; i++)
    res[i] = add(S[i], T[i]);
  return res;
}


long index_non_zero (const Vec<zz_pX> &S){
  if (S.length() == 0) return -1;
  for (long i = 0; i < S.length(); i++)
    if (S[i] != zz_pX(0)) return i;
  return S.length();
}

long index_non_zero (const Vec<Vec<zz_pX>> &S){
  if (S.length() == 0) return -1;
  for (long i = 0; i < S.length(); i++)
    if (index_non_zero(S[i]) != S[i].length())
      return i;
  return S.length();
}

bool is_zero_seq(const Vec<zz_pX> &S){
  if (S.length() == 0) return true;
  for (long i = 0; i < S.length(); i++)
    if (S[i] != zz_pX(0)) return false;
  return true;
}

bool is_zero_seq(const Vec<Vec<zz_pX>> &S){
  if (S.length() == 0) return true;
  return (index_non_zero(S) == S.length());
}

long valuation(const zz_pX &a){
  if (deg(a) == -1) return -1;
  for (long i = 0; i <= deg(a); i++){
    if (a[i] != zz_p(0)) return i;
  }
  return 0; // should never happen
}

void shift(Vec<zz_pX> &S, const long s){
  if (s >= S.length()){
    S.SetLength(0);
    return;
  }
  for (long i = 0; i + s < S.length(); i++)
    S[i] = S[i+s];
  S.SetLength(S.length()-s);

}

void shift(Vec<Vec<zz_pX>> &S, const long s){
  if (s > S.length()){
    S.SetLength(0);
    return;
  }
  for (long i = 0; i + s < S.length(); i++)
    S[i] = S[i+s];
  S.SetLength(S.length()-s);
}

// solves for a in l = a*r mod x^d
void solve(zz_pX &a, const zz_pX &l, const zz_pX &r, const long d){
  long s = valuation(r);
  auto l_shifted = RightShift(l,s);
  auto r_shifted = RightShift(r,s);

  zz_pX mod;
  SetCoeff(mod,d-s,1);
  auto r_inv = InvMod(r_shifted, mod);
  MulMod(a, l_shifted, r_inv, mod);


  // MulMod(a, l_shifted, InvMod(r_shifted, mod), mod);
}

bool check_cancel(const Vec<zz_pX> &S, const Vec<zz_pXY> &gens,
    const long d){
  for (auto &f : gens){
    Vec<zz_pX> res;
    if (d == -1)
      res = mul(S,f);
    else
      res = mulTrunc(S,f,d);
    if (!is_zero_seq(res)) return false;
  }
  return true;
}

bool check_cancel(const Vec<Vec<zz_pX>> &S, const Vec<zz_pXY> &gens,
    const long d){
  for (auto &f : gens){
    Vec<Vec<zz_pX>> res;
    if (d == -1)
      res = mul(S,f);
    else
      res = mulTrunc(S,f,d);
    if (!is_zero_seq(res)) return false;
  }
  return true;
}

void kurakin(const long d, const Vec<zz_pX> &S, Vec<zz_pXY> &gens){
  Vec<Vec<zz_pX>> us;
  us.SetLength(d);
  gens.SetLength(d);

  map<long, Vec<zz_pX>> I_seq;
  map<long, zz_pXY> I_fs;

  // phase 1
  zz_pX running;
  SetCoeff(running, 0, 1);
  for (long t = 0; t < d; t++){
    gens[t] = zz_pXY(running);
    us[t] = mulTrunc(S, running, d);
    LeftShift(running, running, 1);

    long k = index_non_zero(us[t]);
    auto it = I_seq.find(k);
    if (it == I_seq.end() && k >= 0){
      I_seq[k] = us[t];
      I_fs[k] = gens[t];
    } 

    if (verbose){
      cout << "at: " << t << endl;
      cout << "f_t: " << gens[t] << endl;
      cout << "u_s[t]" << endl;
      for (long i = 0; i < us[t].length(); i++)
	cout << us[t][i] << endl;
      cout << endl;
    }

    // we don't need to worry about the case where there
    // is already an entry at k, since the previous entry
    // will always have lower valuation
  }

  // phase 2
  for (long s = 1; s < S.length(); s++){
    for (long t = 0; t < d; t++){
      auto &f = gens[t];
      auto &u = us[t];
      if(verbose){
	cout << "\n\nNEW ITER" << endl;
	cout << "at (s,t): " << s << ", " << t << endl;
	cout << "f: " << f << endl;
	cout << "u:" << endl;

	for (long i = 0; i < us[t].length(); i++)
	  cout << us[t][i] << endl;
      }


      if (u.length() == 0) continue;
      if (is_zero_seq(u)) continue;

      shift(u,1);
      f = shift_y(f,1);

      while (true){
	if (u.length() == 0) break;
	if (is_zero_seq(u)) break;

	if(verbose){

	  cout << "new subiter" << endl;
	  cout << "at (s,t): " << s << ", " << t << endl;
	  cout << "f: " << f << endl;
	  cout << "u:" << endl;

	  for (long i = 0; i < us[t].length(); i++)
	    cout << us[t][i] << endl;
	}

	long k = index_non_zero(u);
	if (k == u.length()) break;

	auto it = I_seq.find(k);
	if (it != I_seq.end() && k >= 0){
	  // means there exists a sequece that
	  // appeared before with the same
	  // number of leading zeros; try
	  // using this to cancel out
	  //the leading non-zero term
	  auto &l = u[k];
	  auto &r = I_seq[k][k];
	  if (valuation(l) < valuation(r))
	    break;
	  else{
	    zz_pX a{1};
	    solve(a,l,r,d);
	    if(verbose){
	      cout << "l: " << l << endl;
	      cout << "r: " << r << endl;
	      cout << "a: " << a << endl;
	    }

	    // update u
	    auto temp = mulTrunc(I_seq[k], -a, d);
	    u = add(u, temp);
	    if(verbose){
	      cout << "u after:" << endl;
	      for (auto & i : u)
		cout << i << endl;
	    }

	    // update f
	    f = trunc_x(f - a * I_fs[k], d);
	    if(verbose){
	      cout << "f after: " << f << endl << endl;
	    }
	  }
	}else{
	  // means there exists no prev seq
	  // that has the same number of leading
	  // zeros
	  break;

	}
      }
    }
    // update the ideals
    for (long t = 0; t < d; t++){
      auto &u = us[t];
      auto &f = gens[t];
      long k = index_non_zero(u);
      if (k < u.length() && k >= 0){
	auto it = I_seq.find(k);
	if (it != I_seq.end()){
	  if (valuation(u[k]) < valuation(I_seq[k][k])){
	    I_seq[k] = u;
	    I_fs[k] = f;
	  }
	}else{
	  I_seq[k] = u;
	  I_fs[k] = f;
	}
      }
    }

  }
}

void modified_kurakin(const long d, const Vec<zz_pX> &S, Vec<zz_pXY> &gens){
  gens.SetLength(0);
  map<long, Vec<zz_pX>> us; // keeps track of the current sequences
  map<long, zz_pXY> fs; // keeps track of the current poly

  // ideal generated by non-zero elements at index k from previous calcs
  map<long, Vec<zz_pX>> I_seq;
  // polynomial that generated the sequence with non-zero at index k
  map<long, zz_pXY> I_fs;

  vector<bool> bitmap; // bitmap[i] = true if x^t is potentially useful
  bitmap.resize(d,false);

  queue<long> Q;
  Q.push(0);
  bitmap[0] = true;
  while(!Q.empty()){
    long t = Q.front();
    Q.pop();
    zz_pX f_x;
    SetCoeff(f_x,t,1);
    zz_pXY f_xy{f_x};
    auto u = mulTrunc(S, f_x, d);

    if(verbose){
      cout << "t: " << t << endl;
      cout << "u: " << u << endl;
      cout << "f_xy: " << f_xy << endl;
    }

    // add to what we are tracking
    us[t] = u;
    fs[t] = f_xy;

    // update the ideals (checking isn't really needed)
    auto k = index_non_zero(u);
    auto it = I_seq.find(k);
    if (it == I_seq.end() && k >= 0){
      I_seq[k] = u;
      I_fs[k] = f_xy;
    }

    // check if we can find t' st (t'+t) < d and x^t' u[k] = 0
    long tprime = d - valuation(u[k]);
    if (tprime+t < d){ 
      // this means that x^(t'+t) is non-zero and x^(t'+t) *S 
      // will create more zero than prev
      bitmap[tprime+t] = true;
      Q.push(tprime+t);
    }
  }

  for (long s = 1; s < S.length(); s++){
    for (long b = 0; b < d; b++) 
      if (bitmap[b]) Q.push(b);
    while (!Q.empty()){
      auto t = Q.front();
      Q.pop();
      auto &f = fs[t];
      auto &u = us[t];

      auto fc = fs[t];
      auto uc = us[t];

      if(verbose){
	cout << "start at (s,t): (" << s << ", " << t << ")" << endl;
	cout << "u: " << u << endl;
	cout << "f: " << f << endl;
      }

      if (is_zero_seq(u)) continue;
      // shift by y
      f = shift_y(f,1);
      shift(u,1);

      // start the subiterations
      while (true){
	if (u.length() == 0) break;
	long k = index_non_zero(u);
	if (k == u.length()) break;
	long tprime = d - valuation(u[k]);
	if (tprime + t < d && bitmap[tprime+t] == false){
	  bitmap[tprime+t] = true;
	  Q.push(tprime+t);
	  fs[tprime+t] = shift_x(fc, tprime);
	  zz_pX x_tprime;
	  SetCoeff(x_tprime,tprime,1);
	  us[t+tprime] = mulTrunc(uc, x_tprime, d);
	}
	auto it = I_seq.find(k);
	if (it != I_seq.end()){
	  // means there is a prev sequence that starts
	  // with same number of non-zeros
	  auto &l = u[k];
	  auto &r = I_seq[k][k];
	  if (valuation(l) < valuation(r)){
	    // see if there is t' such that x^t' l has
	    // valuation >= r
	    tprime = valuation(r) - valuation(l);
	    if (tprime + t < d && bitmap[tprime+t] == false){
	      bitmap[tprime+t] = true;
	      Q.push(tprime+t);
	      fs[tprime+t] = shift_x(fc, tprime);
	      zz_pX x_tprime;
	      SetCoeff(x_tprime,tprime,1);
	      us[t+tprime] = mulTrunc(uc, x_tprime, d);
	    }
	    break;
	  }else{
	    zz_pX a{1}; 
	    solve(a,l,r,d);
	    if (verbose){
	      cout << "a: " << a << endl;
	    }
	    auto temp = mulTrunc(I_seq[k],-a,d);

	    if (verbose){
	      cout << "temp: " << temp << endl;
	      cout << "u before: " << u << endl;
	    }

	    u = add(u,temp);

	    f = trunc_x(f-a*I_fs[k],d);

	    if (verbose){
	      cout << "u after: " << u << endl;
	      cout << "f after: " << f << endl;
	    }
	  }
	}else{
	  // we already checked if there is t' such that
	  // x^t' u_t[k] = 0 earlier, so we can just break
	  break;
	}
      }
    }
    for (long t = 0; t < d; t++){
      if (bitmap[t] == false) continue;
      auto &u = us[t];
      auto &f = fs[t];
      auto k = index_non_zero(u);
      if (k < u.length() && k >= 0){
	auto it = I_seq.find(k);
	if (it != I_seq.end()){
	  if (valuation(u[k]) < valuation(I_seq[k][k])){
	    I_seq[k] = u;
	    I_fs[k] = f;
	  }
	}else{
	  I_seq[k] = u;
	  I_fs[k] = f;
	}
      }
    }

  }
  gens.SetLength(d);
  for (long t = 0; t < d; t++){
    if (bitmap[t] == true) gens[t] = fs[t];
  }
}

void fill_in(const long d, Vec<zz_pXY> &gens){
  zz_pXY zero{zz_pX{0}};
  for (long t = 1; t < gens.length(); t++){
    if (gens[t] == zero){
      gens[t] = shift_x(gens[t-1],1);
      trunc_x(gens[t],gens[t],d); 
    }
  }
}

typedef Vec<zz_pX> Module;
typedef pair<Vec<Module>, zz_pXY> mypair;

long find_dependency(const Vec<Module> &basis, const long d){
  Mat<zz_pX> F;
  long n = basis.length();
  long m = basis[0].length();
  F.SetDims(n,m);
  for (long i = 0; i < n; i++)
    F[i] = basis[i];
  Mat<zz_pX> appbas;
  vector<long> shifts;
  for (long i = 0; i < n; i++) shifts.emplace_back(i);
  pmbasis(appbas,F,d,shifts);
  for (long i = 0; i < n; i++){
    for (long j = 0; j < n; j++){
      if (valuation(appbas[i][j]) == 0){
	return j;
      }
    }
  }
  return -1;
}

Vec<zz_pX> solve(const Vec<zz_pX> &l, const Vec<Vec<zz_pX>> &rs, const long d){
  // construct [rs_0 ... rs_j | -l]^tr
  Mat<zz_pX> F;
  long n = rs.length();
  long m = rs[0].length();
  F.SetDims(rs.length()+1, m);
  for (long i = 0; i < n; i++)
    F[i] = rs[i];
  F[n] = -l;
  if (verbose)
    cout << "F: " << F << endl;

  // compute the left approximant basis at order d
  vector<long> shifts;
  for (long i = 0; i < n+1; i++) shifts.emplace_back(0);
  Mat<zz_pX> appbas;
  pmbasis(appbas,F,d,shifts);
  if (verbose){
    cout << "appbas: " << appbas << endl;
  }

  // read off the solution
  Vec<zz_pX> soln;
  long ind = -1;
  for (long r = 0; r < n+1; r++){
    if (verbose){
      cout << "(r,n+1): " << r << " " << n << endl;
      cout << "looking at: " << appbas[r][n] << endl;
    }
    if (valuation(appbas[r][n]) == 0){
      ind = r;
      break;
    }
  }
  if (ind == -1) return soln;
  zz_pX invC = InvTrunc(appbas[ind][n], d);
  soln.SetLength(n);
  for (long i = 0; i < n; i++)
    soln[i] = MulTrunc(appbas[ind][i],invC,d);
  return soln;
}

Vec<zz_pX> solve(const Vec<zz_pX> &l, const Vec<Vec<zz_pX>> &rs, const long d,
    long &tprime){
  // construct [rs_0 ... rs_j | -l]^tr
  Mat<zz_pX> F;
  long n = rs.length();
  long m = rs[0].length();
  F.SetDims(rs.length()+1, m);
  for (long i = 0; i < n; i++)
    F[i] = rs[i];
  F[n] = -l;
  if(verbose)
    cout << "F: " << F << endl;

  // compute the left approximant basis at order d
  vector<long> shifts;
  for (long i = 0; i < n+1; i++) shifts.emplace_back(0);
  Mat<zz_pX> appbas;
  pmbasis(appbas,F,d,shifts);
  if(verbose){
    cout << "appbas: " << appbas << endl;
  }

  // read off the solution
  Vec<zz_pX> soln;
  long ind = -1;
  for (long r = 0; r < n+1; r++){
    if (verbose){
      cout << "(r,n+1): " << r << " " << n << endl;
      cout << "looking at: " << appbas[r][n] << endl;
    }
    if (valuation(appbas[r][n]) == 0){
      ind = r;
      break;
    }
  }
  if (ind == -1){
    // check for the minimal valuation on the last column
    long min_val = valuation(appbas[0][n]);
    for (long r = 1; r < n+1; r++)
      if (min_val < 0  || valuation(appbas[r][n]) < min_val)
	min_val = valuation(appbas[r][n]);
    if (verbose) cout << "minval: " << min_val << endl;
    tprime = min_val;
    if (verbose) cout << "tprime: " << tprime << endl;
    return soln;
  }
  zz_pX invC = InvTrunc(appbas[ind][n], d);
  soln.SetLength(n);
  for (long i = 0; i < n; i++)
    soln[i] = MulTrunc(appbas[ind][i],invC,d);
  return soln;
}
void kurakin(const long d, const Vec<Module> &S, Vec<zz_pXY> &gens){ 
  map<long, vector<mypair>> I;
  Vec<Vec<Module>> us;
  us.SetLength(d);
  gens.SetLength(d);

  zz_pX running;
  SetCoeff(running, 0, 1);
  for (long t = 0; t < d; t++){
    gens[t] = zz_pXY(running);
    us[t] = mulTrunc(S, running, d);
    LeftShift(running, running, 1);

    long k = index_non_zero(us[t]);
    auto it = I.find(k);
    if (it == I.end() && k >= 0){
      // create a pair and add to I
      mypair p{us[t],gens[t]};
      I[k].emplace_back(p);
    }
  }

  for (long s = 1; s < S.length(); s++){
    // stage 1
    for (long t = 0; t < d; t++){
      auto &f = gens[t];
      auto &u = us[t];


      if (u.length() == 0) continue;
      if (is_zero_seq(u)) continue;

      shift(u,1);
      f = shift_y(f,1);

      if (verbose){
	cout << endl << endl <<  "AT: " << s << " " << t << endl;
	cout << "f: " << f << endl;
	cout << "u: " << u << endl;
      }

      long temp_counter = 0;
      while (temp_counter++ < 2){
	if (u.length() == 0) break;
	if (is_zero_seq(u)) break;

	long k = index_non_zero(u);
	if (k == u.length()) break;

	auto it = I.find(k);
	if (it != I.end() && k >= 0){
	  auto &l = u[k];
	  if(verbose)
	    cout << "l: " << l << endl;
	  Vec<Module> rs;
	  Vec<zz_pXY> rfs;
	  Vec<Vec<Module>> rus;
	  for (unsigned long i = 0; i < I[k].size(); i++){
	    rs.append(I[k][i].first[k]);
	    rus.append(I[k][i].first);
	    rfs.append(I[k][i].second);
	  }
	  if (verbose)
	    cout << "rs: " << rs << endl;
	  // try to solve for l = \sum a_i rs_i
	  auto soln = solve(l,rs,d);
	  if (verbose)
	    cout << "soln: " << soln << endl;
	  if (soln.length() == 0) // no solutions
	    break;
	  // update u and f for the next subiteration
	  for (long i = 0; i < rfs.length(); i++){
	    u = add(u, mulTrunc(rus[i], -soln[i], d));
	    f = f - rfs[i]*zz_pXY{soln[i]};
	    f = trunc_x(f,d);
	  }
	}else{
	  break;
	}
      }
    }
    // stage 2
    for (long t = 0; t < d; t++){
      auto &u = us[t];
      auto &f = gens[t];
      long k = index_non_zero(u);
      if (k < u.length() && k >= 0){
	auto it = I.find(k);
	long nu = I[k].size();
	if (it != I.end()){
	  Vec<Module> basis;
	  for (unsigned long i = 0; i < I[k].size(); i++){
	    basis.append(I[k][i].first[k]);
	  }
	  basis.append(u[k]);
	  long ind = find_dependency(basis, d);
	  if (ind == -1){ // no dependency, just add to basis
	    mypair p{u,f};
	    I[k].emplace_back(p);
	  }else if(ind < nu){
	    // this means that the dependency is one of the vectors
	    // inside I[k]
	    if (verbose){
	      cout << "vector size: " << I[k].size() << endl;
	      cout << "deleting: " << ind << endl;
	    }
	    I[k].erase(I[k].begin()+ind);
	    mypair p{u,f};
	    I[k].emplace_back(p);
	  }
	}else{
	  mypair p{u,f};
	  I[k].emplace_back(p);
	}
      }
    }
  }
}

// returns the minimal (non-zero) valuation of a, returns -1
// if a is zero
long valuation(const Module &a){
  if (is_zero_seq(a)) return -1;
  long min_val = valuation(a[0]);
  for(long i = 1; i < a.length(); i++){
    if (min_val > valuation(a[i])) min_val = valuation(a[i]);
  }
  return min_val;
}

void modified_kurakin(const long d, const Vec<Module> &S, Vec<zz_pXY> &gens){
  gens.SetLength(0);
  map<long, vector<mypair>> I;
  map<long, Vec<Module>> us; // current sequences;
  map<long, zz_pXY> fs; // current polys

  vector<bool> bitmap; // bitmap[i] = true if x^t is potentially useful
  bitmap.resize(d,false);

  if (verbose) cout << "\nStarting mo. kurakin" << endl;

  queue<long> Q;
  Q.push(0); // 0 is always useful
  bitmap[0] = true;
  while(!Q.empty()){
    long t = Q.front();
    Q.pop();
    zz_pX f_x;
    SetCoeff(f_x,t,1);
    zz_pXY f_xy{f_x};
    auto u = mulTrunc(S, f_x, d);

    // add to what we are tracking;
    us[t] = u;
    fs[t] = f_xy;

    // update the ideals
    auto k = index_non_zero(u);
    auto it = I.find(k);
    if (it == I.end() && k >= 0){
      mypair p{u,f_xy};
      I[k].emplace_back(p);
    }

    // check if we can find t' st (t'+t) < d and x^t' u[k] = 0
    long tprime = d - valuation(u[k]);
    if (tprime+t < d && !bitmap[tprime+t]){
      bitmap[tprime+t] = true;
      Q.push(tprime+t);
    }
  }

  for (long s = 1; s < S.length(); s++){
    for (long b = 0; b < d; b++)
      if (bitmap[b]) Q.push(b);
    // stage 1
    while(!Q.empty()){
      auto t = Q.front();
      Q.pop();
      auto &f = fs[t];
      auto &u = us[t];

      // these are the original ones
      auto fc = fs[t]; 
      auto uc = us[t];

      if (is_zero_seq(u)) continue;
      // shift by y
      f = shift_y(f,1);
      shift(u,1);

      if (verbose){ 
	cout << "\nat " << s << " " << t << endl;
	cout << "f: " << f << endl;
	cout << "u: " << u << endl;
      }

      // start of the subiterations
      while(true){
	if (u.length() == 0) break;
	long k = index_non_zero(u);
	if (k == u.length()) break;
	long tprime = d - valuation(u[k]);
	// check if we can cancel the first element just by
	// multiply it with x^tprime
	if (tprime + t < d && bitmap[tprime+t] == false){
	  bitmap[tprime+t] = true;
	  Q.push(tprime+t);
	  fs[tprime+t] = shift_x(fc, tprime);
	  zz_pX x_tprime;
	  SetCoeff(x_tprime,tprime,1);
	  us[t+tprime] = mulTrunc(uc, x_tprime, d);
	}
	// check if we can cancel the element by previous sequences
	auto it = I.find(k);
	if (it != I.end()){
	  auto &l = u[k];
	  Vec<Module> rs;
	  Vec<zz_pXY> rfs;
	  Vec<Vec<Module>> rus;
	  for (unsigned long i = 0; i < I[k].size(); i++){
	    rs.append(I[k][i].first[k]);
	    rus.append(I[k][i].first);
	    rfs.append(I[k][i].second);
	  }
	  long tprime;
	  auto soln = solve(l,rs,d,tprime);
	  if (soln.length() == 0){
	    if (tprime+t >= d || bitmap[tprime+t]) break;
	    bitmap[tprime+t] = true;
	    Q.push(tprime+t);
	    fs[tprime+t] = shift_x(fc,tprime);
	    zz_pX x_tprime;
	    SetCoeff(x_tprime, tprime, 1);
	    us[t+tprime] = mulTrunc(uc, x_tprime, d);
	    break;
	  }
	  for (long i = 0; i < rfs.length(); i++){
	    u = add(u, mulTrunc(rus[i], -soln[i], d));
	    f = f - rfs[i]*zz_pXY{soln[i]};
	    f = trunc_x(f,d);
	  }
	}
	if (it != I.end()){
	}else{
	  break;
	}
      }
    }
    // stage 2
    for (long t = 0; t < d; t++){
      if (bitmap[t] == false) continue; // nothing important here
      auto &u = us[t];
      auto &f = fs[t];
      auto k = index_non_zero(u);
      if (k < u.length() && k >= 0){
	auto it = I.find(k);
	long nu = I[k].size();
	if (it != I.end()){
	  // there is something in the ideal
	  Vec<Module> basis;
	  for (unsigned long i = 0; i < I[k].size(); i++)
	    basis.append(I[k][i].first[k]);
	  basis.append(u[k]);
	  long ind = find_dependency(basis,d);
	  if (ind == -1){
	    mypair p{u,f};
	    I[k].emplace_back(p);
	  }else if (ind < nu){
	    I[k].erase(I[k].begin()+ind);
	    mypair p{u,f};
	    I[k].emplace_back(p);
	  }
	}else{
	  // there's nothing in the ideal
	  mypair p{u,f};
	  I[k].emplace_back(p);
	}
      }
    }
  }
  gens.SetLength(d);
  for (long t = 0;  t< d; t++)
    if (bitmap[t]) gens[t] = fs[t];
}

void minpoly_DAC(const long d, const zz_pXY &SXY, const zz_pX &S0X, 
    const zz_pXY &rhs, zz_pXY &P){
  if (d == 1){
    vector<long> shift;
    shift.emplace_back(0);
    shift.emplace_back(0);
    Mat<zz_pX> appbas;
    Mat<zz_pX> pmat;
    long order = deg(S0X)+1;
    pmat.SetDims(3,1);
    pmat[0][0] = S0X;
    pmat[1][0] = zz_pX{-1};

    // construct rhs mod x
    zz_pX rhs_x;
    for (long i = 0; i <= rhs.degY(); i++){
      if (rhs.rep[0] != zz_p(0))
	SetCoeff(rhs_x, i, rhs.rep[i][0]);
    }
    LeftShift(rhs_x,rhs_x, rhs.degY());
    pmat[2][0] = -rhs_x;
    pmbasis(appbas, pmat, order, shift);
    long r_cor = 0;
    for (long i = 0; i < 3; i++){
      if (deg(appbas[i][2]) == 0){
	r_cor = i;
	break;
      }
    }
    zz_pX px = appbas[r_cor][0]*InvTrunc(appbas[r_cor][2],1);

    P = zz_pX{0};
    for (long i = 0; i <= deg(px); i++)
      P = P + shift_y(zz_pXY{zz_pX{px[i]}}, i);
    if (verbose){
      cout << "rhs_x:" << rhs_x << endl;
      cout << "appbas: " << appbas << endl;
      cout << "px: " << px << endl;
      cout << "P: " << P << endl;
    }
    return;
  }
  zz_pXY P0;
  minpoly_DAC(d/2, trunc_x(SXY,d/2), S0X, trunc_x(rhs,d/2),P0);

  // compute the residue
  zz_pXY r;
  mul(r, SXY, P0);
  r = trunc_x(r,d);
  r = shift_y(r,-rhs.degY());
  r = trunc_y(r,rhs.degY()+1);
  r = rhs - r;

  r = shift_x(r,-d/2);

  zz_pXY P1;
  minpoly_DAC(d/2, trunc_x(SXY,d/2), S0X, r, P1);
  shift_x(P1,P1,d/2);
  P = P0 + P1;
}

void minpoly_nondegenerate(const long d, const Vec<zz_pX> &S, zz_pXY &P){
  // make bivariate polynomial of S for multiplication
  zz_pXY SXY;
  for (long i = 0; i < S.length()-1; i++)
    SXY  = SXY + shift_y(zz_pXY{S[i]}, i);
  zz_pX S0X;
  for (long i = 0; i < S.length()-1; i++)
    SetCoeff(S0X, i, S[i][0]);
  // construct rhs
  zz_pXY rhs;
  for (long i = S.length()/2; i < S.length(); i++)
    rhs = rhs + shift_y(zz_pXY{S[i]},(i-S.length()/2));
  zz_pXY P_res;
  minpoly_DAC(d,SXY,S0X,zz_pXY{zz_pX{0}}-rhs,P_res);
  P = zz_pXY{zz_pX{0}};
  for (long i = 0; i < P_res.degY(); i++){
    P.rep.append(P_res.rep[i]);
  }
  P = P + shift_y(zz_pXY{zz_pX{1}}, S.length()/2);
}

void berlekamp_massey_pmbasis(const long d, const Vec<zz_pX> &S,
    Vec<zz_pXY> &gens){
  long n = S.length()/2;
  Mat<zz_pX> T;
  T.SetDims(n+1,n);
  for (long c = 0; c < n; c++){
    for (long r = 0; r < n+1; r++)
      T[r][c] = S[c+r];
  }
  Mat<zz_pX> appbas;
  VecLong rdeg;
  for (long i = 0; i < n+1; i++) rdeg.emplace_back(0);
  pmbasis(appbas,T,d,rdeg);
  for (long r = 0; r < n+1; r++){
    zz_pXY res;
    for (long c = 0; c < n+1; c++)
      res = res + shift_y(zz_pXY{appbas[r][c]}, c);
    gens.append(res);
  }
}

void berlekamp_massey_pmbasis(const long d, const Vec<Module> &S,
    Vec<zz_pXY> & gens){
  long e = S.length()/2;
  long n = S[0].length();

  Mat<zz_pX> H;
  H.SetDims(e+1, e*n);
  for (long s = 0; s < n; s++){
    for (long c = 0; c < e; c++)
      for (long r = 0; r < e+1; r++)
	H[r][c+s*e] = S[c+r][s];
  }
  Mat<zz_pX> appbas;
  VecLong rdeg;
  for (long i = 0; i < e+1; i++) rdeg.emplace_back(0);
  pmbasis(appbas, H, d, rdeg);
  for (long r= 0; r < e+1; r++){
    zz_pXY res;
    for (long c = 0; c < e+1; c++)
      res = res + shift_y(zz_pXY(appbas[r][c]), c);
    gens.append(res);
  }
  bool rank_check = false;
  if (rank_check){
    Mat<zz_pX> kern;
    kernel_basis(kern, H, rdeg);
    long K = kern.NumRows();
    cout << "rank of F: " << e+1-K << endl;

    Mat<zz_p> Hc;
    Hc.SetDims(e+1,e*n);
    for (long r = 0; r < e+1; r++)
      for (long c = 0; c < e*n; c++)
	Hc[r][c] = coeff(H[r][c], 0);
    long rank = gauss(Hc);
    cout << "rank of F(0): " << rank << endl;
  }
}

void berlekamp_massey_pmbasis_hnf(const long d, const Vec<Module> &S,
    Vec<zz_pXY> & gens){
  long e = S.length()/2;
  long n = S[0].length();

  Mat<zz_pX> H;
  H.SetDims(e+1, e*n);
  for (long s = 0; s < n; s++){
    for (long c = 0; c < e; c++)
      for (long r = 0; r < e+1; r++)
	H[r][c+s*e] = S[c+r][s];
  }
  Mat<zz_pX> appbas;
  VecLong rdeg;
  for (long i = 0; i < e+1; i++){
    rdeg.emplace_back(i*d);
  }
  pmbasis(appbas, H, d, rdeg);
  for (long r= 0; r < e+1; r++){
    zz_pXY res;
    for (long c = 0; c < e+1; c++)
      res = res + shift_y(zz_pXY(appbas[r][c]), c);
    gens.append(res);
  }
}

void generate_right_seq(Vec<zz_pX> &S, const structured_lzz_pX &A,
    const Mat<zz_pX> b, const bool prec){
}


void right_mul(Vec<zz_p> &out, const Vec<zz_p> seq, const Vec<zz_p> v){
  zz_pX A;
  zz_pX B;
  long e = seq.length()/2;
  for (long i = 0; i < 2*e; i++)
    SetCoeff(A, i, seq[i]);
  for (long i = 0; i < e; i++)
    SetCoeff(B, i, v[e-i-1]);
  double t = GetWallTime();
  auto C = A*B;
  cout << "time single mul: " << GetWallTime()-t << endl;
  out.SetLength(e+1);
  for (long i = 0; i < e+1;i++)
    out[i] = coeff(C,i+(e-1));
}

void right_mul(Mat<zz_p> &out, const Vec<zz_p> seq, const Mat<zz_p> V){
  long e = seq.length()/2;
  out.SetDims(e+1,e+1);
  for (long i = 0; i < e+1; i++){
    Vec<zz_p> cur;
    Vec<zz_p> v;
    for (long j = 0; j < e; j++){
      v.append(V[j][i]);
    }
    right_mul(cur, seq, v);
    for (long j = 0; j < e+1; j++){
      out[j][i] = cur[j];
    }
  }
}


// computes res = sum Hs[i] * mats[i]
void right_mul(Mat<zz_p> & res, const Vec<zz_pX> &polys_H, const Mat<zz_pX> &polmat){
  
  long e = res.NumCols()-1; 
  /*
  Vec<zz_pX> polys_H; // the polynomials of Hs

  for (long i = 0; i < n; i++){
    zz_pX poly;
    for (long j = 0; j < 2*e; j++)
      SetCoeff(poly, j, Hs[i][j]);
    polys_H.append(poly);
  }
  Mat<zz_pX> polmat; // polmat[i][j] stores the i-th column of j-th matrix in mats
  polmat.SetDims(e+1,n);
  for (long j = 0; j < n; j++){
    for (long i = 0; i < e+1; i++){
      zz_pX poly;
      for (long c = 0; c < e; c++)
	SetCoeff(poly, c, mats[j][e-c-1][i]);
      polmat[i][j] = poly;
    }
  }
  */
  auto prod = polmat * polys_H;

  res.SetDims(e+1,e+1);
  for (long r = 0; r < e+1; r++){
    for (long c = 0; c < e+1; c++){
      res[r][c] = coeff(prod[c], r+(e-1));
    }
  }

}

void compress(Mat<zz_pX> &mat, const long d, const Vec<Module> &S){
  long e = S.length()/2;
  long n = S[0].length();

  for (long i = 0; i < e+1; i++)
    for (long j = 0; j < e+1; j++)
      mat[i][j] = zz_pX(zz_p(0));

  Vec<Mat<zz_p>> rands;
  for (long i = 0; i < n; i++)
    rands.append(random_mat_zz_p(e,e+1));

  Mat<zz_pX> rand_polmat;
  rand_polmat.SetDims(e+1,n);
  for (long j = 0; j < n; j++){
    for (long i = 0; i < e+1; i++){
      zz_pX poly;
      for (long c = 0; c < e; c++)
        SetCoeff(poly, c, rands[j][e-c-1][i]);
      rand_polmat[i][j] = poly;
    }
  }

  Mat<zz_pX> polHs;
  polHs.SetDims(n, d);
  for (long i = 0; i < d; i++){
    for (long j = 0; j < n; j++){
      for (long s = 0; s < 2*e; s++){
	SetCoeff(polHs[j][i], s, coeff(S[s][j], i));
      }
    }
  }

  auto prod = rand_polmat * polHs;
  
  mat.SetDims(e+1,e+1);
  for (long i = 0; i < d; i++)
    for (long r = 0; r < e+1; r++)
      for (long c = 0; c < e+1; c++)
	  SetCoeff(mat[r][c], i, coeff( prod[c][i], r+(e-1)));

/*
  for (long i = 0; i < d; i++){
    Vec<zz_pX> seqs; // for sequences of coefficient at x^i
    for (long j = 0; j < n; j++){
      zz_pX seq;
      for (long s = 0; s < 2*e; s++){
	SetCoeff(seq, s, coeff(S[s][j], i));
      }
      seqs.append(seq);
    }
    // create the Hankel matrices for each seq in seqs
    Mat<zz_p> mat_i;
    mat_i.SetDims(e+1,e+1);
    right_mul(mat_i, seqs, rand_polmat);
    for (long r = 0; r < e+1; r++)
      for (long c = 0; c < e+1; c++)
	SetCoeff(mat[r][c], i, mat_i[r][c]);
  }
*/
 /*
  double time_ver = GetWallTime();
  Mat<zz_pX> H;
  H.SetDims(e+1, e*n);
  for (long s = 0; s < n; s++){
    for (long c = 0; c < e; c++)
      for (long r = 0; r < e+1; r++)
	H[r][c+s*e] = S[c+r][s];
  }

  Mat<zz_p> proj;
  proj.SetDims(n*e, e+1);
  for (long i = 0; i < n; i++){
    for (long r = 0; r < e; r++)
      for (long c = 0; c < e+1; c++)
	proj[i*n+r][c] = rands[i][r][c];
  }
  prod = H*proj;
  cout << "time to verify: " << GetWallTime()-time_ver << endl;
*/
}

void berlekamp_massey_pmbasis_compressed(const long d, const Vec<Module> &S,
    Vec<zz_pXY> & gens){

  long e = S.length()/2;

  Mat<zz_pX> compressed;
  compressed.SetDims(e+1,e+1);
  compress(compressed, d, S);

  Mat<zz_pX> appbas;
  VecLong rdeg;
  for (long i = 0; i < e+1; i++) rdeg.emplace_back(0);
  pmbasis(appbas, compressed, d, rdeg);
  for (long r= 0; r < e+1; r++){
    zz_pXY res;
    for (long c = 0; c < e+1; c++)
      res = res + shift_y(zz_pXY(appbas[r][c]), c);
    gens.append(res);
  }

}





























