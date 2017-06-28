#pragma once

#include <math.h>
#include <iostream>
#include <vector>

using namespace std;

inline double mean(vector<double> data){
  double sum=0;
  double size = data.size();
  for (int i=0;i<data.size();i++){
    sum += data[i];
  }
  return sum/size;
}

inline double standardDeviation(vector<double> data){
  double m = mean(data);
  double sum = 0;
  double size = data.size();

  for (int i=0;i<data.size();i++){
    sum += pow(data[i]-m,2);
  }
  return sqrt(sum/size);
};

inline double CoefficientVariation(vector<double> data){
  double m = mean(data);
  double sd = standardDeviation(data);
  return (sd/m);
}
