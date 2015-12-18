/*
Example of the use of a vector of mpfr
*/

#include<iostream>
#include<sstream>
#include<mpfr.h>
#include<math.h>
#include<vector>
#include <stdlib.h>     // rand

#define LIST_OF_WAVELETS "haar","db1","db2","db3","db4","db5","db6","db7","db8","db9","db10","db11","db12","db13","db14", "db15","bior1.1","bio1.3","bior1.5","bior2.2","bior2.4","bior2.6","bior2.8","bior3.1","bior3.3","bior3.5","bior3.7","bior3.9","bior4.4","bior5.5","bior6.8","coif1","coif2","coif3","coif4","coif5"
// To compile and run: g++ exampleMpfr2.cpp -o file -lmpfr -lgmp; ./file

typedef struct{
  mpfr_t* dummy;
  mpfr_t* vectorOfmpfr;
  //mpfr_t* dummy,vectorOfmpfr; // DON'T WORK!!
}Chain;

typedef struct{
  std::vector<Chain*> chains;
}Run;

int main (int argc, char **argv) {

  mpfr_set_default_prec(256);
  int a=5;
  
  std::string listOfWavelet[] = {LIST_OF_WAVELETS};
  std::cout << "List of wavelets: ";
  for( unsigned int a = 0; a < sizeof(listOfWavelet)/sizeof(listOfWavelet[0]); a = a + 1 )
    std::cout << listOfWavelet[a] << " ";
  std::cout << std::endl;

  Run run;
  Chain* chain = new Chain();  // Create a pointer on a Chain
  chain->vectorOfmpfr = new mpfr_t[a];
  run.chains.push_back(chain); // Add this chain to the run
    
  mpfr_t temp;
  mpfr_init_set_d(temp,4.0,MPFR_RNDN);

  for(int iz=0;iz<a;iz++) {
    mpfr_init_set_d(run.chains[0]->vectorOfmpfr[iz],0.0,MPFR_RNDN);
    mpfr_add(run.chains[0]->vectorOfmpfr[iz],run.chains[0]->vectorOfmpfr[iz],temp,MPFR_RNDN);
  }

  for(int iz=0;iz<a;iz++) {
   std::cout << "iz:" << std::endl;
   mpfr_out_str(stdout,10,256,run.chains[0]->vectorOfmpfr[iz],MPFR_RNDN); // stdout, decimal base, precision, rounding option
   std::cout << std::endl;
   std::cout << std::endl;
  }
  
  delete[] run.chains[0]->vectorOfmpfr;
  
  return 0;
}
