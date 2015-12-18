/*
Example of how to compare exp(-BIGNUMBER) and exp(-OTHERBIGNUMBER) by using mpfr
*/
#include<iostream>
#include<sstream>
#include<mpfr.h>
#include <math.h>

// To compile and run: g++ exampleMpfr.cpp -o file -lmpfr -lgmp; ./file

int main (int argc, char **argv) {

    long double truc=28500;
    long double autreTruc=28500.0000000001;
    int prec=256;

    // == Convert doubles to str ==
    std::ostringstream strs;
    strs.precision(prec);
    strs << -truc;
    std::string strstr = strs.str();
    const char* str1 = strstr.c_str();
    std::ostringstream strs2;
    strs2.precision(prec);
    strs2 << -autreTruc;
    std::string strstr2 = strs2.str();
    const char* str2 = strstr2.c_str();
    std::cout << "str1:" << str1 << std::endl;
    std::cout << "str2:" << str2 << std::endl;
    std::cout << std::endl;
    // == end of convert to doubles to str ==

    mpfr_set_default_prec(prec);
    mpfr_t a,b,c,d; // Declare four numbers
    mpfr_inits2 (prec,a,b,c,d,(mpfr_ptr) 0); // and initialize them with precision 256
    mpfr_set_str(a,str1,10,MPFR_RNDN); // Variable where we store, value as a char*, decimal base, rounding option
    mpfr_set_str(b,str2,10,MPFR_RNDN); // Variable where we store, value as a char*, decimal base, rounding option
    // mpfr_set_str(a,"15", 10,MPFR_RNDN); // Variable where we store, value as a char*, decimal base, rounding option
 
     std::cout << "a:" << std::endl;
    mpfr_out_str(stdout,10,prec,a,MPFR_RNDN); // stdout, decimal base, precision, rounding option
    std::cout << std::endl;
    std::cout << "b:" << std::endl;
    mpfr_out_str(stdout,10,prec,b,MPFR_RNDN); // stdout, decimal base, precision, rounding option
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout.precision(prec);
 
    mpfr_exp(c,a,MPFR_RNDN);
    mpfr_exp(d,b,MPFR_RNDN);
    std::cout << "c:" << std::endl;
    mpfr_out_str(stdout,10,prec,c,MPFR_RNDN); // stdout, decimal base, precision, rounding option
    std::cout << std::endl;
    std::cout << "d:" << std::endl;
    mpfr_out_str(stdout,10,prec,d,MPFR_RNDN); // stdout, decimal base, precision, rounding option
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout.precision(prec);
    if (mpfr_cmp(c,d) == 1)
      std::cout << "exp(-" << truc << ") > exp(-" << autreTruc << ")" << std::endl;
    else
      std::cout << "exp(-" << truc << ") <= exp(-" << autreTruc << ")" << std::endl;
    //mpfr_out_str(stdout,10,256,b,MPFR_RNDN); // stdout, decimal base, precision, rounding option
    
    mpfr_clears(a,b,c,d,(mpfr_ptr) 0); // Clear variables

    return 0;
}
