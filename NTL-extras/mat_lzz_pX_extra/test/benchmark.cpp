#include <iostream>
#include <vector>
#include <limits>
#include <fstream>
#include <sstream>

using namespace std;

struct data{
  static bool return_size;  
  
  long p;
  long size;
  long deg;
  vector<double> timings;
  
  data(long p, long size, long deg, 
       double t1, double t2, double t3, double t4):
  p{p},size{size},deg{deg}{
    timings.emplace_back(t1);
    timings.emplace_back(t2);
    timings.emplace_back(t3);
    timings.emplace_back(t4);
  }
  
  long val() const{
    if (return_size)
      return size;
    else
      return deg;
  }
};

bool data::return_size = true;

struct tree{
  const data d;
  tree *left = nullptr;
  tree *right = nullptr;
  tree *subtree = nullptr;
  
  tree(const data& d, 
       tree* left, 
       tree* right,
       tree* subtree)
       :d(d),left(left),right(right),subtree(subtree){}
  
  ~tree(){
    delete left;
    delete right;
    delete subtree;
  }
  
  void preorder_print() const{
    cout << d.deg << " ";
    if (left != nullptr) left->preorder_print();
    if (right != nullptr) right->preorder_print();
  }
  
  void print() const{
    cout << "holds size: " << d.size << endl;
    cout << "subtree: " << endl;
    subtree->preorder_print();
    cout << endl;
    
    if(left != nullptr) left->print();
    if(right != nullptr) right->print();
  }
  
  int find_best(long size, long deg){
    auto *stp = find_closest(true,size,this);
    auto *dtp = find_closest(false,deg,stp->subtree);
    
    const auto &timings = dtp->d.timings;
    long ind = 0;
    double min_time = timings[0];
    for (long i = 1; i < timings.size(); i++){
      if (min_time == -1){
        ind = i;
        min_time = timings[i];
      }
      if (timings[i] != -1 && min_time > timings[i]){
        ind = i;
        min_time = timings[i];
      }
    }
    cout << "size: " << size << ", deg: " << deg <<
    ", min_time: " << min_time << endl;
    return ind;
  }
  
  tree* find_closest(bool size, long n, tree* start){
    tree* at;
    long min = numeric_limits<long>::max();
    data::return_size = size;
    find_closest(start,at,min,n);
    return at;
  }
  
  long abs(long n){
    if (n < 0) return -n;
    return n;
  }
  
  void find_closest(tree* tp, tree* &at, long &min, const long val){
    if (tp == nullptr) return;
    if (abs(tp->d.val() - val) < min){
      at = tp;
      min = abs(tp->d.val() - val);
    }
    if (val < tp->d.val())
      find_closest(tp->left,at,min,val);
    else
      find_closest(tp->right,at,min,val);
  }
};

tree* create_tree_deg(const vector<data> &v){
  if (v.size() == 0) return nullptr;
  if (v.size() == 1){
    return new tree{v[0],nullptr,nullptr,nullptr};
  }
  
  auto n = v.size()/2;
  long deg = v[n].deg;
  vector<data> left,right;
  for (long i = 0; i < v.size(); i++){
    if (i !=n){
      if (v[i].deg < deg)
        left.emplace_back(v[i]);
      else
        right.emplace_back(v[i]);
    }
  }
  return new tree{v[n],
                  create_tree_deg(left),
                  create_tree_deg(right),
                  nullptr};
}

tree* create_tree_size(const vector<data> &v){
  if (v.size() == 0)
    return nullptr;
  if (v.size() == 1){
    return new tree{v[0],nullptr,nullptr,create_tree_deg(v)};
  }
  
  auto n = v.size()/2;
  auto &size = v[n].size;
  vector<data> left,right;
  vector<data> eq;
  for (long i = 0; i < v.size(); i++){
    if (i !=n){
      if (v[i].size < size)
        left.emplace_back(v[i]);
      else if (v[i].size > size)
        right.emplace_back(v[i]);
      else
        eq.emplace_back(v[i]);
    }
  }
  
  return new tree{v[n],
                  create_tree_size(left),
                  create_tree_size(right),
                  create_tree_deg(eq)};
}


vector<data> parse (string filename){
  vector<data> result;
  ifstream ifs{filename};
  string s;
  getline(ifs,s);
  while(getline(ifs,s)){
    istringstream iss{s};
    char c;
    long p,size,len;
    double t1,t2,t3,t4;
    iss >> p >> c >> size >> c >> len >> c
    >> t1 >> c >> t2 >> c >> t3 >> c >> t4;
    result.emplace_back(data{p,size,len,t1,t2,t3,t4});
  }  
  return result;
}

int main(){
  vector<data> v = parse("test.csv");
  vector<data> large_prime, small_prime, fft_prime;
  for (auto &i: v){
    if (i.p == 0) fft_prime.emplace_back(i);
    else if (i.p == 23068673) small_prime.emplace_back(i);
    else large_prime.emplace_back(i);
  }       
  
  auto *tp0 = create_tree_size(fft_prime);
  auto *tp1 = create_tree_size(small_prime);
  auto *tp2 = create_tree_size(large_prime);
  long p,size,deg;
  cout << "ENTER p(0,1,2), size, deg" << endl;
  while (cin >> p >> size >> deg){
    auto l = 0; 
    if (p == 0)
      l = tp0->find_best(size,deg);
    else if (p == 1)
      l = tp1->find_best(size,deg);
    else
      l = tp1->find_best(size,deg);
    cout << "index: " << l << endl << endl;
    cout << "ENTER p(0,1,2), size, deg" << endl;
  }
  delete tp0;
  delete tp1;
  delete tp2;
}




























