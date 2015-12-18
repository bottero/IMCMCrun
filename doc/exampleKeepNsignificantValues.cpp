#include<iostream>
#include<math.h>
#include<vector>
#include <iomanip>      // std::setprecision
#include <algorithm>    // std::sort

// To compile and run: g++ exampleKeepNsignificantValues.cpp -o file -lmpfr -lgmp; ./file

void findthresh(std::vector<double> vector, int N, double& t)
// At the end t is the Nth biggest value in vector
{
  std::sort(vector.begin(), vector.end(), std::greater<double>());
  t = vector.at(N-1);
}

void keepNsignificantValues(std::vector<double>* vector, int N)
// Keep the N biggest absolute values in vector. Put the others to 0
// !! Warning !! for example if vector = 10 21 12 12 7 8 and N = 2
// at the end vector = 0 21 12 12 0 0. We could fix that.
{
  int notZeros = 0;
  double thresh = 0.0;
  std::vector<double> temp_dwtoutput;
  for (unsigned int i = 0; i < vector->size();i++) {
    double temp = fabs((*vector)[i]);
    temp_dwtoutput.push_back(temp);
  }
  findthresh(temp_dwtoutput,N, thresh); // thresh is the Nth biggest value of temp_dwtoutput
  for (unsigned int i = 0; i < vector->size();i++) {
    double temp = fabs((*vector)[i]);
    if (temp < thresh or notZeros >= N)
      (*vector)[i] = 0.0;
    else
      notZeros++;
  }
}

int main (int argc, char **argv) {

  static const double arr[] ={-103.842,18636.9,19811.7,-103.842,-293.706,-293.706,-293.706,-293.706,-103.842,-103.842,-103.842,-103.842,-103.842,-103.842,-103.842,-103.842,-36.7139,-36.7139,-36.7139};    
  std::vector<double> vec (arr, arr + sizeof(arr) / sizeof(arr[0]) );

  std::cout << "vec:" << std::endl;
  for (std::vector<double>::const_iterator j = vec.begin(); j != vec.end(); ++j)
    std::cout << *j << ' ';
  std::cout << std::endl << std::endl;
  
  keepNsignificantValues(&vec,6);
  
  std::cout << "vec:" << std::endl;
  for (std::vector<double>::const_iterator j = vec.begin(); j != vec.end(); ++j)
    std::cout << *j << ' ';
  std::cout << std::endl << std::endl;
  
  return 0;
}
