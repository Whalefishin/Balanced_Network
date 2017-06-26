#include <string>
#include <iostream>
#include <utility>
#include <ctime>
#include <vector>

using namespace std;

int main(){

vector<double> test_v;

test_v.push_back(1);

for (int i=0;i<10;i++){
  test_v.push_back(test_v.back()+1);
}

for (int i=0;i<test_v.size();i++){
  cout << test_v[i]<< endl;
}
  return 0;
}
