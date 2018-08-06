#include <iostream>
#include <vector>

using namespace std;
int main(int argc, char* argv[])
{
    long sizes[] = 
        {1,2,3,4,5,6,7,8,9,10,15,20,25,30,35,
         40,45,50,55,60,65,70,75,80,85,90,95,
         100,200,300,400};
    long ls = sizeof(sizes)/sizeof(long);
    long degrees[] = 
        {1,2,3,4,5,6,7,8,9,10,15,20,25,30,35,
         40,45,50,55,60,65,70,75,80,85,90,95,
         100,200,300,400,500};
    long ld = sizeof(degrees)/sizeof(long);
  
    vector<long> v1{sizes,sizes+ls};
    vector<long> v2{degrees,degrees+ld};
  
    cout << "#!/bin/bash" << endl;
    for (auto &i : v1)
        for (auto &j : v2){
            if (i^2*j < 200000000)
                cout << "./test_multiply " << i << " " << j << " >> " << argv[1] << endl;
        }
}
